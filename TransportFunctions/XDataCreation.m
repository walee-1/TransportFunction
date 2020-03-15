(* Wolfram Language package *)



XDataCreation[{xmin_?NumericQ,xmax_?NumericQ,xOff_},BinN_]:=Table[x+xOff,{x,xmin,xmax,(xmax-xmin)/(BinN-1)}]

XYDataCreation::usage = "Creates XY data with min and max values and an Offset for each dimension.";
XYDataCreation[{xmin_?NumericQ,xmax_?NumericQ,xOff_,xbins_},{ymin_?NumericQ,ymax_?NumericQ,yOff_,ybins_}]:=
Flatten[Table[
	{x+xOff,y+yOff},
		{x,xmin,xmax,(xmax-xmin)/(xbins-1)},
		{y,ymin,ymax,(ymax-ymin)/(ybins-1)}],1]