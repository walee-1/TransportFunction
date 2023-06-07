(* ::Package:: *)

Get["General/easeFunctions.m"];
Get["FileReadingModules/commonTables.m"];


depthCutOffMod[array_] := 
 Module[{energy = {}, depth = {}}, 
  Table[If[array[[i, 5]] == "-", 
    If[array[[i, 6]]*0.1 > 30, AppendTo[energy, array[[i, 4]]]], 
    If[array[[i, 5]]*0.1 > 30, AppendTo[energy, array[[i, 3]]]]], {i, 
    1, array // Length}]; 
  Table[If[array[[i, 5]] == "-", 
    If[array[[i, 6]]*0.1 > 30, AppendTo[depth, array[[i, 6]]*0.1]], 
    If[array[[i, 5]]*0.1 > 30, 
     AppendTo[depth, array[[i, 5]]*0.1]]], {i, 1, array // Length}]; 
  Return[{energy, depth}]]

backscatteringReaderMod[fileName_, numberIons_] := 
  Module[{dataFile, backscatPercentage, totalBackscat}, 
   dataFile = Import[fileName, "Table", HeaderLines -> 34]; 
   totalBackscat = 
    dataFile[[Position[dataFile, "Backscattered"][[1, 1]], 4]]; 
   backscatPercentage = totalBackscat/numberIons; 
   Return[{totalBackscat, backscatPercentage}]];
   
   rangeReadingMod[fileName_, headerLines_: 24, germanGuard_: False] := 
 Module[{dataFile, depth, len, decDelimiter}, 
  decDelimiter = If[germanGuard, ",", "."]; 
  dataFile = 
   Import[fileName, "Table", HeaderLines -> headerLines, 
    NumberPoint -> decDelimiter]; 
  If[NumericQ[Check[dataFile[[dataFile // Length, 2]], False]], 
   len = dataFile // Length, len = (dataFile // Length) - 1];
  depth = Table[dataFile[[i, 2]], {i, 1, len}];
  Return[{dataFile, depth}]]
  
  BSJoinerMod[fileName1_, fileName2_] := 
 Block[{joinedList, list1, list2, totalIons, totalBS, BSPercent},
  list1 = fileReadingModNew2[fileName1];
  list2 = fileReadingModNew2[fileName2];
  totalIons = list1[[5]] + list2[[5]];
  totalBS = (list1[[3]] + list2[[3]]);
  BSPercent = totalBS/totalIons;
  joinedList = {Join[list1[[1]], list2[[1]]], BSPercent, totalBS, 
    list1[[4]], totalIons};
  Return[joinedList]
  ](*deprecated MOD, use srimMultiRunImporterInstead*)
  
  
BSJoinerMod::usage="Deprecated Mod, please use srimMultiRunImporter instead as it can handle everything with far more ease"  

srimMultiRunImporter[loc_] := 
 Block[{runNo, runTable, totalIons, totalBs, BSPercent, finalList},
  runNo = StringSplit[#, {"/"}][[-1]] & /@ FileNames["Run*", loc];
  Which[Length[runNo] > 1,
   runTable = 
    Table[fileReadingModNew2[loc <> i <> "/BACKSCAT.txt"], {i, runNo}];
   totalIons = Total[runTable[[All, 5]]];
   totalBs = Total[runTable[[All, 3]]];
   BSPercent = totalBs/totalIons;
   finalList = {Flatten[runTable[[All, 1]], 1], BSPercent, 
      totalBs, runTable[[1, 4]], totalIons};,
   Length[runNo] == 1, 
   finalList = 
     fileReadingModNew2[loc <> runNo[[1]] <> "/BACKSCAT.txt"];,True,
   Return[Print["There are no Run folders at the location"]];
   ];
  Return[finalList];
  ]

srimMultiRunImporter::usage="srimMultiRunImporter[location] loads SRIM Backscattering files for all possible runs for a specific energy and angle and combines them.
	Parameters:
		Location is the location of the run folders. The folders have to start with Run*.
	Also works for single runs, so just overall better SRIM backscattering file reader function"


fileReadingModNew2[fileName_, numberIons_] := 
 Block[{dataFile, numberIons2, totalBackscat, backscatPercentage, len,
    col, decDelimiter, bLength, dataFileBefore, dataFileAfter, key, 
   germanGuard, char},
  If[StringContainsQ[fileName, "BACK"] || 
    StringContainsQ[fileName, "Back"] || 
    StringContainsQ[fileName, "back"], char = "B", char = "T"];
  dataFile = OpenRead[fileName];
  germanGuard = 
   StringContainsQ[Table[ReadLine[dataFile], {i, 13}][[13]], ","];
  Close[dataFile]; Clear[dataFile];
  decDelimiter = If[germanGuard, ",", "."];
  dataFile = 
   Import[fileName, "Table", HeaderLines -> 12, 
    NumberPoint -> decDelimiter];
  len = dataFile // Length;
  bLength = Length[Position[dataFile, char]];
  dataFileBefore = 
   Delete[dataFile[[1 ;; bLength]], Position[dataFile, char]];
  dataFileAfter = 
   Transpose@
    ReplacePart[Transpose[dataFile[[bLength + 1 ;; All]]], 
     1 -> ToExpression[
       StringReplace[char -> ""][
        dataFile[[bLength + 1 ;; All, 1]]]]];
  dataFile = Join[dataFileBefore, dataFileAfter];
  numberIons2 = Round[dataFile[[Length[dataFile]]][[1]], 10000];
  backscatPercentage = len/numberIons;
  totalBackscat = len;
  If[char == "B", 
   key = {{1, "ionNo"}, {2, "AtomNo"}, {3, 
      "energy at leaving (ev)"}, {4, " -sign"}, {5, 
      "xPos when decl. BS(A)"}, {6, "lateralY (A)"}, {7, 
      "lateralZ (A)"}, {8, "Cos(X)"}, {9, "Cos(Y)"}, {10, "Cos(Z)"}}, 
   key = {{1, "ionNo"}, {2, "AtomNo"}, {3, 
      "energy at leaving (ev)"}, {4, "xPos when decl. BS(A)"}, {5, 
      "lateralY (A)"}, {6, "lateralZ (A)"}, {7, "Cos(X)"}, {8, 
      "Cos(Y)"}, {9, "Cos(Z)"}}];
  Return[{dataFile, backscatPercentage, totalBackscat, key, 
    numberIons}]]
fileReadingModNew2[fileName_] := 
 Block[{dataFile, numberIons, totalBackscat, backscatPercentage, len, 
   rangeString, decDelimiter, bLength, dataFileBefore, dataFileAfter, 
   key, germanGuard, char},
  If[StringContainsQ[fileName, "BACK"] || 
    StringContainsQ[fileName, "Back"] || 
    StringContainsQ[fileName, "back"], char = "B", char = "T"];
  Which[char == "B", 
   rangeString = StringReplace[fileName, "BACKSCAT" -> "RANGE_3D"], 
   char == "T", 
   rangeString = StringReplace[fileName, "TRANSMIT" -> "RANGE_3D"], 
   True, Print["Wrong File given"]; Return[]];
  dataFile = OpenRead[fileName];
  germanGuard = 
   StringContainsQ[Table[ReadLine[dataFile], {i, 13}][[13]], ","];
  Close[dataFile]; Clear[dataFile];
  decDelimiter = If[germanGuard, ",", "."];
  dataFile = 
   Import[fileName, "Table", HeaderLines -> 12, 
    NumberPoint -> decDelimiter];
  len = dataFile // Length;
  bLength = Length[Position[dataFile, char]];
  dataFileBefore = 
   Delete[dataFile[[1 ;; bLength]], Position[dataFile, char]];
  dataFileAfter = 
   Transpose@
    ReplacePart[Transpose[dataFile[[bLength + 1 ;; All]]], 
     1 -> ToExpression[
       StringReplace[char -> ""][
        dataFile[[bLength + 1 ;; All, 1]]]]];
  dataFile = Join[dataFileBefore, dataFileAfter];
  numberIons = 
   ToExpression[
    StringSplit[Import[rangeString, {"Data", -Range[1]}]][[1]]];
  backscatPercentage = len/numberIons;
  totalBackscat = len;
  If[char == "B", 
   key = {{1, "ionNo"}, {2, "AtomNo"}, {3, 
      "energy at leaving (ev)"}, {4, " -sign"}, {5, 
      "xPos when decl. BS(A)"}, {6, "lateralY (A)"}, {7, 
      "lateralZ (A)"}, {8, "Cos(X)"}, {9, "Cos(Y)"}, {10, "Cos(Z)"}}, 
   key = {{1, "ionNo"}, {2, "AtomNo"}, {3, 
      "energy at leaving (ev)"}, {4, "xPos when decl. BS(A)"}, {5, 
      "lateralY (A)"}, {6, "lateralZ (A)"}, {7, "Cos(X)"}, {8, 
      "Cos(Y)"}, {9, "Cos(Z)"}}];
  Return[{dataFile, backscatPercentage, totalBackscat, key, 
    numberIons}]]
fileReadingModNew2::usage="fileReadingModNew2[filename, totalIons] reads SRIM backscattering or transmitting files and gives the output with {dataFile, BSPercent, totalBS, key, totalIons}
fileReadingModNew2[filename] reads SRIM backscattering or transmitting files and gives the output with {dataFile, BSPercent, totalBS, key, totalIons}.
fileReadingModNew2[filename] is the better version but it required the existances of RANGE3D file, otherwise it can't work.
fileReadingModNew2[filename, totalIons] rounds the ions in BACKSCAT.txt to the nearest 10000, so be careful if the total number of simulated ions were in a different range."


transmissionFileReadingMod[fileName_, numberIons_] := 
 Module[{dataFile, energy, depth, transmitPercentage, len, 
   totalTransmit}, 
  dataFile = Import[fileName, "Table", HeaderLines -> 12];
  If[NumericQ[Check[dataFile[[dataFile // Length, 4]], False]], 
   len = dataFile // Length, len = (dataFile // Length) - 1];
  energy = 
   Table[If[dataFile[[i, 3]] == 1, dataFile[[i, 4]], 
     dataFile[[i, 3]]], {i, 1, len}];
  depth = 
   Table[If[dataFile[[i, 3]] == 1, dataFile[[i, 5]], 
     dataFile[[i, 4]]], {i, 1, len}];
  transmitPercentage = (dataFile // Length)/numberIons*100;
  totalTransmit = dataFile // Length;
  Return[{dataFile, energy, depth, transmitPercentage, totalTransmit, numberIons}]]
  


exyzReadingMod2[fileName_, germanGuard_: False] := 
 Module[{dataFile, data, fileStream}, 
  fileStream = OpenRead[fileName];
  Table[ReadLine[fileStream], {i, 15}];
  dataFile = 
   ReadList[
    fileStream, {Number, Number, Number, Number, Number, Number, 
     Number}];
  data = Gather[dataFile, First[#1] === First[#2] &];
  Return[data]]

exyzReadingMod[fileName_, germanGuard_: False] := 
 Module[{dataFile, file, data}, 
  If[germanGuard, 
   dataFile = 
    Import[fileName, "Table", HeaderLines -> 15, 
     "NumberPoint" -> ","], 
   dataFile = Import[fileName, "Table", HeaderLines -> 15]];
  file = Table[{dataFile[[i, 1]], dataFile[[i, 2]], dataFile [[i, 3]],
      dataFile[[i, 4]], dataFile[[i, 5]], dataFile[[i, 6]], 
     dataFile[[i, 7]]}, {i, dataFile // Length}];
  data = Gather[file, First[#1] === First[#2] &];
  Return[data]]
  
exyzReadingMod::usage = "Better to use exyzReadingMod2[fileName, germanGuard] instead as it uses lower level file handling and as a result less memory intensive"
exyzReadingMod2::usage = "exyzReadingMod2[fileName, germanGuard] reads the EXYZ file from SRIM and gives it in the form of a list.
Default value of germanGuard is False and should only be set to True if the delimiter for the numbers is a \",\""


mdRangeModule[fileName_] := Block[{list, listNormed},
  list = Import[fileName, "Table"];
  listNormed = 
   Transpose[{list[[All, 1]]/10, 
     list[[All, 2]]/Total[list[[All, 2]]]}];
  Return[{listNormed, list}];]
  
  mdRangeEnergyModule[fileName_] := 
 Block[{list, ionEn, nuclEn, listNormed},
  list = Import[fileName, "Table"];
  ionEn = Transpose[{list[[All, 1]], list[[All, 3]]}];
  nuclEn = Transpose[{list[[All, 1]], list[[All, 2]]}];
  Return[{ionEn, nuclEn, list}];
  ]
  
  cceFormula[z_?NumericQ, gamma_?NumericQ, tau_?NumericQ] := 
 1 - gamma*Exp[-z/tau]
 
 gammaVals = {0.6, 0.9};
tauVals = {50, 100};(*tau vals are in nm*)

imsilHistogramReaderMod[fileName_] := 
 Block[{importedFile, indexInfo = {1}, i, j = 1, k = 0, mainTable},
  importedFile = Import[fileName, "Table"];
  i = importedFile[[1, 2]] + 2;
  While[j < Length[importedFile], j = j + i + k; 
   If[j < Length[importedFile], i = importedFile[[j, 2]]]; k = 2; 
   AppendTo[indexInfo, j]];
  mainTable = 
   Table[importedFile[[
     indexInfo[[i]] + 2 ;; indexInfo[[i + 1]] - 1]], {i, 
     Length[indexInfo] - 1}];
  Return[mainTable]
  ]
  
  imsilBackReader[fileName_] := 
 Block[{backscatteredNo, backscatteredErr, file}, 
  file = OpenRead[fileName];
  Find[file, "Backscattered atoms per ion"];
  backscatteredNo = 
   ToExpression[StringSplit[Table[ReadLine[file], {4}][[4]]][[5]]];
  backscatteredErr = 
   ToExpression[StringSplit[Table[ReadLine[file], {5}][[5]]][[5]]];
  Close[file];
  Clear[file];
  Return[{backscatteredNo, backscatteredErr}];
  ]
  
imsilHistogramReaderMod2[fileName_, simIons_: 500000] := 
 Block[{importedFile, bsIons, bsInfo, indexInfo = {1}, i, j = 1, 
   k = 0, mainTable}, importedFile = Import[fileName, "Table"]; 
  i = importedFile[[1, 2]] + 2; 
  While[j < Length[importedFile], j = j + i + k; 
   If[j < Length[importedFile], i = importedFile[[j, 2]]]; k = 2; 
   AppendTo[indexInfo, j]]; 
  mainTable = 
   Table[importedFile[[
     indexInfo[[i]] + 2 ;; indexInfo[[i + 1]] - 1]], {i, 
     Length[indexInfo] - 1}];
  bsInfo = 
   imsilBackReader[
     StringSplit[fileName, ".his"][[1]] <> ".out"][[1]];
  bsIons = bsInfo*simIons;
  mainTable = Append[mainTable, bsIons];
  Return[mainTable]]

histoFixerMod[list_, bsIons_, simIons_: 500000] := 
 Block[{val, calcIons, retHisto, noIons}, 
  val = Table[{list[[i - 1, 1]], 
     list[[i, 2]] (list[[i, 1]] - list[[i - 1, 1]])}, {i, 3, 
     Length[list] - 1, 2}];
  noIons = Floor[Total[val[[All, 2]]*simIons]];
  If[noIons != bsIons, Print["ERROR"]; Return[0]];
  retHisto = 
   Transpose[{val[[All, 1]], val[[All, 2]]*simIons/noIons}];
  Return[retHisto];
  ]
  
  
  
  imsilRangeReader[fileName_] := 
 Block[{RangeNo, RangeErr, file}, file = OpenRead[fileName];
  Find[file, "Vertical moments (A)"];
  RangeNo = 
   ToExpression[
    StringSplit[Table[ReadLine[file], {4}][[4]], "|"][[3]]];
  RangeErr = 
   ToExpression[
    StringSplit[Table[ReadLine[file], {5}][[5]], "|"][[3]]];
  Close[file];
  Clear[file];
  Return[{RangeNo, RangeErr}];
  ]
  
  trajReaderMod[fileName_, debugFlag_: False] := 
 Block[{importFile, fileStream, noIons, ionStartList, returnList}, 
  If[debugFlag == True,
   fileStream = OpenRead[fileName];
   importFile = 
    ReadList[
     fileStream, {Number, Number, Number, Number, Number, Number, 
      Number, Number, Number, Number}];
   ionStartList = 
    Append[Flatten[Position[importFile[[All, 8; 10]], _?(# == 1 &)]], 
     Length[importFile] + 1];
   Print[ionStartList];
   noIons = Length[ionStartList] - 1;
   Print[noIons];
   returnList = 
    Table[importFile[[
      ionStartList[[i]] ;; ionStartList[[i + 1]] - 1]], {i, noIons}],
   fileStream = OpenRead[fileName];
   importFile = 
    ReadList[
     fileStream, {Number, Number, Number, Number, Number, Number, 
      Number, Number, Number, Number}];
   ionStartList = 
    Append[Flatten[Position[importFile[[All, 8; 10]], _?(# == 1 &)]], 
     Length[importFile] + 1];
   noIons = Length[ionStartList] - 1;
   returnList = 
    Table[
     importFile[[ionStartList[[i]] ;; ionStartList[[i + 1]] - 1]], {i,
       noIons}];
   ];
  Clear[importFile, ionStartList];
  Close[fileStream];
  Return[returnList]
  ]
  
  binningModIMSIL[filename_] := 
 Block[{datax, datay, return, file, binLength, evTable2, evTable}, 
  file = imsilHistogramReaderMod[filename];
  binLength = file[[1, 5, 1]] - file[[1, 4, 1]];
  datax = Table[file[[1, i, 1]], {i, 4, Length[file[[1]]], 2}];
  datay = Table[file[[1, i, 2]], {i, 4, Length[file[[1]]] - 1, 2}];
  return = 
   Table[{(datax[[i]] + datax[[i + 1]])/2, datay[[i]]}, {i, 
     Length[datay]}];
  evTable2 = {datax/10, datay};
  evTable = 
   Table[{(datax[[i]] + datax[[i + 1]])/2, datay[[i]]*binLength}, {i, 
     Length[datay]}];
  Return[{return, evTable, evTable2}]
  ]
  
  ehModIMSIL[filename_, window_, gammaVal_: gammaVals[[2]], 
  tauVal_: tauVals[[2]]] := 
 Block[{datax, datay, return, file, binLength, evTable, evHoles}, 
  file = imsilHistogramReaderMod[filename];
  binLength = file[[1, 5, 1]] - file[[1, 4, 1]];
  datax = Table[file[[1, i, 1]], {i, 4, Length[file[[1]]], 2}];
  datay = Table[file[[1, i, 2]], {i, 4, Length[file[[1]]] - 1, 2}];
  return = 
   Table[{(datax[[i]] + datax[[i + 1]])/2, datay[[i]]}, {i, 
     Length[datay]}];
  evTable = 
   Table[{(datax[[i]] + datax[[i + 1]])/2, datay[[i]]*binLength}, {i, 
     Length[datay]}];
  evHoles = 
   Table[{(datax[[i]] + datax[[i + 1]])/
      2, (evTable[[i, 2]]/
        3.6)*(cceFormula[datax[[i]]/10 - window, gammaVal, tauVal] + 
         cceFormula[datax[[i + 1]]/10 - window, gammaVal, tauVal])/
       2}, {i, FirstPosition[datax, x_ /; x > window*10][[1]], 
     Length[evTable]}];
  Return[{return, evTable, evHoles, binLength}]
  ]
  
  ehModSRIM[list_, window_, gammaVal_: gammaVals[[2]], 
  tauVal_: tauVals[[2]]] := Block[{binLength, evTable, evHoles},
  binLength = list[[2, 1]] - list[[1, 1]];
  evTable = 
   Table[{list[[i, 1]], list[[i, 2]]*binLength}, {i, Length[list]}];
  evHoles = 
   Table[{list[[i, 
      1]], (evTable[[i, 2]]/
        3.6)*(cceFormula[list[[i, 1]]/10 - window, gammaVal, tauVal] +
          cceFormula[list[[i + 1, 1]]/10 - window, gammaVal, tauVal])/
       2}, {i, FirstPosition[list[[All, 1]], x_ /; x > window*10][[
      1]], Length[evTable] - 1}];
  Return[{evTable, evHoles}]
  ]
  
  ehModSRIM[list_, window_] := Block[{binLength, evTable, evHoles},
  binLength = list[[2, 1]] - list[[1, 1]];
  evTable = 
   Table[{list[[i, 1]], list[[i, 2]]*binLength}, {i, Length[list]}];
  evHoles = 
   Table[{list[[i, 
      1]], (evTable[[i, 2]]/
        3.6)*(cceFormula[list[[i, 1]]/10 - window, gammaVals[[2]], 
          tauVals[[2]]] + 
         cceFormula[list[[i + 1, 1]]/10 - window, gammaVals[[2]], 
          tauVals[[2]]])/2}, {i, 
     FirstPosition[list[[All, 1]], x_ /; x > window*10][[1]], 
     Length[evTable] - 1}];
  Return[{evTable, evHoles}]
  ]
  
  diffListwErrorMod[list1_, list2_, angleList_: anglesTable] := 
  Block[{diffList}, 
   If[Dimensions[list1] == Dimensions[list2], 
    diffList = 
     Table[{angleList[[th]], list1[[en, th, 2]] - list2[[en, th, 2]], 
       Sqrt[list1[[en, th, 3]]^2 + list2[[en, th, 3]]^2]}, {en, 
       Length[list1]}, {th, Length[list1[[en]]]}]; Return[diffList], 
    Print["Dimensions of Lists not equal"]; Return[0]]];
    
    diffListwErrorMod2[list1_, list2_, angleList_: anglesTable] := 
  Block[{diffList}, 
   If[Dimensions[list1] == Dimensions[list2], 
    diffList = 
     Table[{angleList[[th]], 
       Abs[(list1[[en, th, 2]] - list2[[en, th, 2]])]/
        list2[[en, th, 2]], 
       Sqrt[(list2[[en, th, 3]]/
             list2[[en, th, 2]])^2 + (list1[[en, th, 3]]/
             list1[[en, th, 2]])^2]*
        list1[[en, th, 2]]/list2[[en, th, 2]]}, {en, 
       Length[list1]}, {th, Length[list1[[en]]]}]; Return[diffList], 
    Print["Dimensions of Lists not equal"]; Return[0]]];
    
    diffListwErrorMod12[list1_, list2_, angleList_: anglesTable] := 
  Block[{diffList}, 
   If[Dimensions[list1] == Dimensions[list2], 
    diffList = 
     Table[{angleList[[th]], 
       Abs[NumberForm[
         1 - list1[[en, th, 2]]/list2[[en, th, 2]], {2, 2}]]}, {en, 
       Length[list1]}, {th, Length[list1[[en]]]}]; Return[diffList], 
    Print["Dimensions of Lists are not equal"]; Return[0]]];
    
    
    cceWeightingModulePerIon[list_, binSize_: 100, 
  gammaVal_: gammaVals[[2]], tauVal_: tauVals[[2]]] := 
 Block[{len = Length[list], ehPairRaw, ehPairWeighted, totalehPair, 
   initDepth, finalDepth, dummy}, 
  If[list[[1, 2]] != 0, ehPairRaw = (1/3.6*#) & /@ list[[All, 2]];
    dummy = list[[All, 3]][[All]];
    initDepth = Table[dummy[[i]], {i, 1, Length[dummy]}];
    finalDepth = 
     Table[initDepth[[i + 1]], {i, 1, Length[initDepth] - 1}];
    If[Length[finalDepth] != 0, 
     AppendTo[finalDepth, finalDepth[[Length[finalDepth]]] + binSize];
      PrependTo[finalDepth, 0], AppendTo[finalDepth, 0 + binSize]; 
     PrependTo[finalDepth, 0]];
    ehPairWeighted = 
     Table[ehPairRaw[[i]]*
       cceWeightFunc[initDepth[[i]], finalDepth[[i]], gammaVal, 
        tauVal], {i, len}];
    totalehPair = Total[ehPairWeighted];
    Return[{list[[1, 1]], totalehPair, ehPairWeighted, finalDepth}], 
    Return[{list[[1, 1]], 0, {0}, {0}}]];
  ]

vecCalc[x1_, y1_, z1_, x2_, y2_, z2_] := {x2 - x1, y2 - y1, z2 - z1}

firstPointCalc[x1_, y1_, z1_, x2_, y2_, z2_] := 
 Block[{res1, x3 = 100, y3, z3, XDummy, percentX}, 
  res1 = vecCalc[x1, y1, z1, x2, y2, z2]; 
  percentX = (x2 - x3)/res1[[1]]; XDummy = percentX*res1; 
  Return[XDummy]]

(*THIS IS THE BETTER MODULE AND FASTER*)
depthModAll2[list_, window_, binSize_] := 
 Block[{vectorMag, vector, depth, ionizEn, totalIonizEn, maxDepth, 
   resList = {}, dummyBins, start = 0, maxList, i},
  maxList = Ceiling[Max[list[[All, 3]]] - window];
  totalIonizEn = Table[0, {i, 0, maxList, binSize}];
  maxDepth = 
   Table[i, {i, 0, Ceiling[Max[list[[All, 3]]] - window, binSize], 
     binSize}];
  Do[If[list[[j, 3]] >= window,
    If[j == 2, 
     vector = 
      firstPointCalc[list[[j - 1, 3]], list[[j - 1, 4]], 
       list[[j - 1, 5]], list[[j, 3]], list[[j, 4]], list[[j, 5]]],
     vector = 
      vecCalc[list[[j - 1, 3]], list[[j - 1, 4]], list[[j - 1, 5]], 
       list[[j, 3]], list[[j, 4]], list[[j, 5]]]];
    vectorMag = Norm[vector];
    depth = list[[j, 3]] - window;
    ionizEn = vectorMag*list[[j, 6]];
    i = Ceiling[depth/binSize];
    totalIonizEn[[i]] += ionizEn;
    ],
   {j, 2, Length[list]}];
  resList = 
   Table[{list[[1, 1]], totalIonizEn[[j]], maxDepth[[j]]}, {j, 
     Length[totalIonizEn]}];
  Return[resList];
  ]

mainMod[enList_, ionNo_, window_, binSize_: 100] := 
 Block[{start, list, temp, tempList2}, 
  tempList2 = enList[[ionNo]]; 
  start = FirstPosition[tempList2[[All, 3]], _?(# >= window &)][[
    1]](*finds the first position of where depth>dead area*); 
  If[start==1,start=2];
  If[start == "NotFound", Return[{{tempList2[[1, 1]], 0, 0, 0}}], 
   list = Table[
     tempList2[[i]], {i, (start - 1), 
      Length[tempList2]}];(*starts from active silicon - 1 depth*)
   temp = depthModAll2[list, window, binSize]; Return[temp];
   ]
  ]
  
  
  (*THERE IS AN ERROR IN THE CODE BELOW! RUN WITH 1 bin size otherwise \
electron hole pairs are a factor of binlength lower!*)

histogramMod2[list_, binSize_, flag_: 0] := 
 Block[{max, len, bins, result, result2}, 
  max = Round[Max[list[[All, 2]]], binSize];
  bins = Table[i, {i, 1, max, binSize}];
  len = Length[bins];
  If[bins[[len]] < max, AppendTo[bins, bins[[len]] + binSize]];
  PrependTo[bins, 0];
  result = HistogramList[list[[All, 2]], {bins}];
  result2 = 
   Table[{(result[[1, i]] + result[[1, i + 1]])/2, 
     result[[2, i]]}, {i, Length[result[[2]]]}];
  Return[result2];
  ]
  
  efficiencyMod[histoList_, noIons_: 50000] := 
 Block[{len, binsWoPedes, signalWoPedes, resBinX, resBinY, effRes, 
   undetectRes}, binsWoPedes = histoList[[2 ;; All, 1]];
  signalWoPedes = histoList[[2 ;; All, 2]];
  len = Length[signalWoPedes];
  resBinY = 
   Table[Sum[signalWoPedes[[j]], {j, i, len}]/noIons, {i, 1, len}];
  resBinX = binsWoPedes;
  effRes = Transpose[{resBinX, resBinY}];
  undetectRes = Transpose[{resBinX, 100.*(1 - resBinY)}];
  Return[{effRes, undetectRes}];]
  
  cceWeightFunc[initDepth_, finalDepth_, gammaVal_, tauVal_] := 
 Block[{initDepthnm, finalDepthnm, depthnm, normVal, normValPos, 
   weight}, initDepthnm = initDepth/10;
  finalDepthnm = finalDepth/10;
  If[initDepthnm > 400, weight = 1,
   weight = (cceFormula[initDepthnm, gammaVal, tauVal] + 
       cceFormula[finalDepthnm, gammaVal, tauVal])/2];
  Return[weight];]
  
  
  
  