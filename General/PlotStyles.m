(* Wolfram Language package *)
imageWrite[var_, path_: imgDir] := 
  Block[{}, Export[path <> var <> ".png", ToExpression[var]]; 
   Export[path <> var <> ".svg", ToExpression[var]];];

texWrite[var_, path_: writeDir] := 
  Export[path <> var <> ".tex", ToExpression[var]];
  
  SetOptions[ListLogLogPlot, ImageSize -> 800, Frame -> True, 
  GridLines -> Automatic, 
  LabelStyle -> Directive[FontSize -> 20, Bold, Black], 
  FrameStyle -> Thick, FrameTicksStyle -> Thick];

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
