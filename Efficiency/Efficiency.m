(* ::Package:: *)

(* Wolfram Language Raw Program *)
Get["FileReadingModules/FileReadingModules.m"]


cceFormula[z_?NumericQ, gamma_?NumericQ, tau_?NumericQ] := 
 1 - gamma*Exp[-z/tau]

gammaVals = {0.6, 0.9};
tauVals = {50, 100};(*tau vals are in nm*)

trajSorterModImsil[list_, entranceWindow_, binSize_: 100] := 
 Block[{depthBins, maxDepth, listSize, maxBin, binList, enLossTab, 
   posIndex, negativePos, flag = True},
  listSize = Length[list];
  maxDepth = Max[list[[All, 3]]];
  If[maxDepth > entranceWindow,
   maxBin = Ceiling[maxDepth - entranceWindow, binSize];
   binList = Table[i, {i, 0, maxBin, binSize}];
   enLossTab = Table[0, {i, Length[binList] - 1}];
   Do[If[list[[i, 3]] > entranceWindow,
     posIndex = 
      FirstPosition[binList, 
         x_ /; x >= list[[i, 3]] - entranceWindow][[1]] - 1;
     enLossTab[[posIndex]] = enLossTab[[posIndex]] + list[[i, 7]];
     ]
    , {i, listSize}];
   negativePos = Flatten[Position[enLossTab, x_ /; x < 0], 1];
   If[negativePos == {}, flag = False];
   While[flag == True,
    If[negativePos != {},
     Do[
       Which[i != 1, 
         enLossTab[[i - 1]] = enLossTab[[i - 1]] + enLossTab[[i]]; 
         enLossTab[[i]] = 0,
         i == 1, enLossTab[[i]] = 0];
       , {i, negativePos}];
     ];
    negativePos = Flatten[Position[enLossTab, x_ /; x < 0], 1]; 
    If[negativePos == {}, flag = False]];
   Return[{binList, enLossTab}],
   Return[{{0, binSize}, {0}}]
   ]]

cceModImsil[list_, gammaVal_: gammaVals[[2]], tauVal_: tauVals[[2]], 
  binSize_: 100] := Block[{cceTable, totalIons, len},
  len = Length[list[[2]]];
  cceTable = 
   Table[list[[2, i]]/
      3.6*(cceFormula[list[[1, i]]/10, gammaVal, tauVal] + 
        cceFormula[list[[1, i + 1]]/10, gammaVal, tauVal])/2, {i, len}];
  Return[{list[[1]], cceTable}];
  ]

avgEhpairDepthModImsil[list_, totalIons_] := 
 Block[{maxBin, binWidth, binsX, binsY, pos, plotTable},
  maxBin = Max[list[[All, 1]]];
  binWidth = list[[1, 1, 2]] - list[[1, 1, 1]];
  binsX = Table[i, {i, 0, maxBin, binWidth}];
  binsY = Table[0, {i, 1, Length[binsX] - 1}];
  Do[
   Do[binsY[[j]] = binsY[[j]] + list[[i, 2, j]]
    , {j, Length[list[[i, 2]]]}]
   , {i, Length[list]}];
  binsY = binsY/totalIons;
  plotTable = 
   Table[{(binsX[[i + 1]] + binsX[[i]])/2/10, binsY[[i]]}, {i, 
     Length[binsY]}];
  Return[{plotTable}]
  ]

histogramModImsil[list_, binSize_] := 
 Block[{histogramList, plotList, totalSignal},
  totalSignal = Table[Total[list[[i, 2]]], {i, Length[list]}];
  histogramList = HistogramList[totalSignal, {binSize}];
  plotList = 
   Table[{(histogramList[[1, i]] + histogramList[[1, i + 1]])/2, 
     histogramList[[2, i]]}, {i, Length[histogramList[[2]]]}];
  Return[{histogramList, plotList}];
  ]

effiencyModImsil[histoList_, noIons_: 50000] := 
 Block[{len, binsWoPedes, signalWoPedes, resBinX, resBinY, errBinY, 
   effRes, effResErr, undetectRes, errSignal, resBinYwErr, 
   undtctResWErr}, binsWoPedes = histoList[[2 ;; All, 1]];
  signalWoPedes = histoList[[2 ;; All, 2]];
  errSignal = 
   Table[Sqrt[signalWoPedes[[i]]], {i, Length[signalWoPedes]}];
  len = Length[signalWoPedes];
  resBinY = 
   Table[Sum[signalWoPedes[[j]], {j, i, len}]/noIons, {i, 1, len}];
  errBinY = 
   Sqrt[Table[Sum[errSignal[[j]]^2, {j, i, len}], {i, 1, len}]]/noIons;
  resBinX = binsWoPedes;
  effRes = Transpose[{resBinX, resBinY}];
  resBinYwErr = 
   Table[Around[resBinY[[i]], errBinY[[i]]], {i, Length[resBinY]}];
  effResErr = Transpose[{resBinX, resBinYwErr}];
  undetectRes = Transpose[{resBinX, 100.*(1 - resBinY)}];
  undtctResWErr = Transpose[{resBinX, 100.*(1 - resBinYwErr)}];
  Return[{effRes, undetectRes, effResErr, undtctResWErr, errBinY}];]
  
  imsilTrajWriterFunc[listName_, windowA_: 150] := 
 Block[{trajData, ionData, fileNameSplit, en, ang, runNo, 
   outputFile},
  $HistoryLength = 1;
  fileNameSplit = 
   StringSplit[StringTrim[StringSplit[listName, "/"][[-1]], ".tra"], 
    "_"];
  runNo = fileNameSplit[[2]];
  en = StringTrim[fileNameSplit[[3]], "en"];
  ang = StringTrim[fileNameSplit[[4]], "ang"];
  trajData = trajReaderMod[listName];
  ionData = 
   Table[trajSorterModImsil[trajData[[i]], windowA], {i, 
     Length[trajData]}];
  outputFile = 
   "ImsilPlgadDataEn" <> en <> "p" <> runNo <> "ang" <> ang <> 
    "IonList.mx";
  Export[FourTbDir <> 
    "MD_Simulations/IMSIL/For_Paper/pLGAD_Paper/"<>outputFile, ionData];
  Clear[trajData, ionData];
  Return[outputFile];
  ]
ClearAll[imsilTrajImportMod]
imsilTrajImportMod[enSrchStr_, 
  dir_: FourTbDir <> "MD_Simulations/IMSIL/For_Paper/pLGAD_Paper/"] :=
  Block[{locFileLocation, fileNames}, If[StringQ[enSrchStr],
   locFileLocation = FileNames["*En" <> enSrchStr <> "p*.mx", dir];
   
   Which[$OperatingSystem=="Windows",
   fileNames = StringSplit[#, {"/", ".mx"}][[-1]] & /@ locFileLocation;
   fileNames=StringSplit[#, {"\\", ".mx"}][[-1]] & /@ fileNames,
   $OperatingSystem=="Unix",
   fileNames = StringSplit[#, {"/", ".mx"}][[-1]] & /@ locFileLocation;
   ];
   fileNames = StringDelete[#, "."] & /@ fileNames;
   Table[If[! ListQ[ToExpression[fileNames[[i]]]],
     With[{\[FormalS] = Symbol[fileNames[[i]]]}, \[FormalS] = 
       Import[locFileLocation[[i]]]]], {i, Length[fileNames]}];
   Return[fileNames];
   , Print["Error: Enter a search string"]; Return[]]
  ]
  
  imsilTrajFileJoinerMod[fileNames_, ang_] := Block[{selFiles, data},
  selFiles = 
   Select[fileNames, 
    StringMatchQ[#, "*ang" <> ToString[ang] <> "*"] &];
  If[selFiles == {},
   selFiles = 
    Select[fileNames, 
     StringMatchQ[#, "*p" ~~ NumberString ~~ "Ion*"] &];
   If[selFiles == {},
    Return[
     Print["The list does not contain entries from this angle or has \
another naming format"]]
    ](*end of second if*)
   ];(*end of first if*)
  
  data = Flatten[DeleteDuplicates[ToExpression[#] & /@ selFiles], 1];
  Return[data];
  ]
  
  
  totalWStatErr[list_] := 
 Block[{var}, var = Total[list]; Return[{var, Sqrt[var]}]]

histogramModImsilErr[list_, binSize_] := 
 Block[{histogramList, plotList, totalSignal, hlErr, plotListErr}, 
  totalSignal = 
   Table[totalWStatErr[list[[i, 2]]], {i, Length[list]}];
  histogramList = HistogramList[totalSignal[[All, 1]], {binSize}];
  hlErr = 
   Table[Sqrt[histogramList[[2, i]]], {i, Length[histogramList[[2]]]}];
  plotList = 
   Table[{(histogramList[[1, i]] + histogramList[[1, i + 1]])/2, 
     histogramList[[2, i]]}, {i, Length[histogramList[[2]]]}];
  plotListErr = 
   Table[{(histogramList[[1, i]] + histogramList[[1, i + 1]])/2, 
     Around[histogramList[[2, i]], hlErr[[i]]]}, {i, 
     Length[histogramList[[2]]]}];
  Return[{histogramList, plotList, plotListErr, totalSignal, hlErr}];]
