(* Wolfram Language package *)


XYBinIntBoundaries[{xmin_,xmax_,xOff_,xbins_},{ymin_,ymax_,yOff_,ybins_}]:=
Module[
	{xlength,ylength,xbinlength,ybinlength},
	xlength=xmax-xmin;
	ylength=ymax-ymin;
	xbinlength=xlength/xbins;
	ybinlength=ylength/ybins;
	Flatten[Table[
		{
			xOff+xmin+{i*xbinlength,(i+1)*xbinlength},
			yOff+ymin+{j*ybinlength,(j+1)*ybinlength}
		},
		{i,0,xbins-1},
		{j,0,ybins-1}
	],1]
		
]