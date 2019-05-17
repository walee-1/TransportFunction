(* Wolfram Language package *)



XDataCreation[{xmin_?NumericQ,xmax_?NumericQ},BinN_]:=Table[x,{x,xmin,xmax,(xmax-xmin)/(BinN-1)}]