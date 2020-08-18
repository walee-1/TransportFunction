(* Wolfram Language package *)
Get["General/PathSettings.m"]
Get["General/Colors.m"]

saveFile[] := NotebookSave@EvaluationNotebook[];

  rebinningMod[list_, binsRebin_] := 
 Block[{binLength, retList, remainderBins, remainderList}, 
  binLength = list[[2, 1]] - list[[1, 1]];
  remainderBins = Mod[Length[list], binsRebin];
  If[remainderBins == 0,
   retList = 
    Table[{Sum[list[[i + j, 1]], {j, 0, binsRebin - 1}]/binsRebin, 
      Sum[list[[i + j, 2]], {j, 0, binsRebin - 1}]}, {i, 1, 
      Length[list], binsRebin}], 
   retList = 
    Table[{Sum[list[[i + j, 1]], {j, 0, binsRebin - 1}]/binsRebin, 
      Sum[list[[i + j, 2]], {j, 0, binsReabin - 1}]}, {i, 1, 
      Length[list] - remainderBins, binsRebin}]; 
   remainderList = 
    Join[list[[-remainderBins ;;]], 
     Table[{list[[Length[list], 1]] + i*binLength, 0}, {i, 1, 
       remainderBins}]]; 
   retList = 
    AppendTo[
     retList, {Total[remainderList[[All, 1]]]/Length[remainderList], 
      Total[remainderList[[All, 2]]]}];];
  Return[retList];
  ]
  
    PercCalc[val_, per_, round_: 0.001] := 
 Block[{err, return}, err = val*per; 
  return = Sort[{val + err, val - err}]; Return[Round[return, round]]]
(* Wolfram Language Raw Program *)
saveFile[] := NotebookSave@EvaluationNotebook[];

    statErrorCalc[percent_, TotalNumber_] := 
 percent/Sqrt[(percent*TotalNumber)]
 
 relErrCalc[valEXP_, valOrig_] := (valEXP - valOrig)/valOrig
  
SetAttributes[imageWrite, HoldFirst];
imageWrite[var_, file_: 0 path_: imgDir] := 
  Block[{varName, fileName}, varName = SymbolName[Unevaluated@var];
   If[Evaluate[file] == 0, fileName = varName,fileName=file];
   Export[path <> fileName <> ".png", ToExpression[varName]];
   Export[path <> fileName <> ".svg", ToExpression[varName]];];


   
   texWrite[var_, path_: writeDir] := 
  Export[path <> var <> ".tex", ToExpression[var]];

SetOptions[ListPlot, ImageSize -> 800, Frame -> True, 
  GridLines -> Automatic, 
  LabelStyle -> Directive[FontSize -> 20, Bold, Black], 
  FrameStyle -> Thick, FrameTicksStyle -> Thick];

SetOptions[Plot, ImageSize -> 800, Frame -> True, 
  GridLines -> Automatic, 
  LabelStyle -> Directive[FontSize -> 20, Bold, Black], 
  FrameStyle -> Thick, FrameTicksStyle -> Thick];

SetOptions[ListLogPlot, ImageSize -> 800, Frame -> True, 
  GridLines -> Automatic, 
  LabelStyle -> Directive[FontSize -> 20, Bold, Black], 
  FrameStyle -> Thick, FrameTicksStyle -> Thick];
  
  errorBarMod[prob_, totalN_] := 
 100/totalN*
  Sqrt[totalN*
    prob*(1 - 
      prob)] (*because we convert the mean to a percentage by \
dividing by totalN and multiplying by 100, 
  so similar thing should be done for the errorbars... basic \
propogation of errors is a thing dumbass*)

scTicks[start_, end_, step_, power_] := 
 Table[{i*10^power, 
   Style[ToString[i] <> "\[Times]\!\(\*SuperscriptBox[\(10\), \(" <> 
     ToString[power] <> "\)]\)"]}, {i, start, end, step}]

texWrite[var_, path_: writeDir] := 
  Export[path <> var <> ".tex", ToExpression[var]];
  
  SetOptions[ListLogLogPlot, ImageSize -> 800, Frame -> True, 
  GridLines -> Automatic, 
  LabelStyle -> Directive[FontSize -> 20, Bold, Black], 
  FrameStyle -> Thick, FrameTicksStyle -> Thick];