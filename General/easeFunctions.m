(* ::Package:: *)

(* Wolfram Language package *)
Get["General/PathSettings.m"]
Get["General/Colors.m"]
Get["General/PlottingOptions.m"]
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
  
saveFile[] := NotebookSave@EvaluationNotebook[];


SetAttributes[imageWrite, HoldFirst];
imageWrite[var_, file_: 0 ,path_: imgDir] := 
  Block[{varName, fileName}, varName = SymbolName[Unevaluated@var];
   If[Evaluate[file] == 0, fileName = varName,fileName=file];
   Export[path <> fileName <> ".png", ToExpression[varName]];
   Export[path <> fileName <> ".svg", ToExpression[varName]];];
imageWrite::usage = "Writes images to the imgDir folder in svg and png formats
Arguments:
	var = Variable name (without strings)
	file = File name, by default it is the variable name
	path = if you want a different path than imgDir"

   
   texWrite[var_, path_: writeDir] := 
  Export[path <> var <> ".tex", ToExpression[var]];

LogarithmicScaling[x_, min_, max_] := Log[x/min]/Log[max/min]

texWrite[var_, path_: writeDir] := 
  Export[path <> var <> ".tex", ToExpression[var]];
  
  SetOptions[ListLogLogPlot, ImageSize -> 800, Frame -> True, 
  GridLines -> Automatic, 
  LabelStyle -> Directive[FontSize -> 20, Bold, Black], 
  FrameStyle -> Thick, FrameTicksStyle -> Thick];
