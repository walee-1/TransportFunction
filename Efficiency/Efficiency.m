(* Wolfram Language Raw Program *)
Get["FileReadingModules/FileReadingModules.m"]
Get["General/PathSettings.m"]
Get["General/PlotStyles.m"]

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
      3.6*(cceFormula[list[[1, i]], gammaVal, tauVal] + 
        cceFormula[list[[1, i + 1]], gammaVal, tauVal])/2, {i, len}];
  Return[{list[[1]], cceTable}];
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

effiencyModImsil[histoList_, noIons_: 1000] := 
 Block[{len, binsWoPedes, signalWoPedes, resBinX, resBinY, effRes, 
   undetectRes},
  binsWoPedes = histoList[[2 ;; All, 1]];
  signalWoPedes = histoList[[2 ;; All, 2]];
  len = Length[signalWoPedes];
  resBinY = 
   Table[Sum[signalWoPedes[[j]], {j, i, len}]/noIons, {i, 1, len}];
  resBinX = binsWoPedes;
  effRes = Transpose[{Reverse[resBinX], resBinY}];
  undetectRes = Transpose[{resBinX, 100.*(1 - resBinY)}];
  Return[{effRes, undetectRes}];
  ]
