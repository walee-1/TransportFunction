(* ::Package:: *)

(* Wolfram Language package *)

Get["Spectra/ProtonSpectrumNachtmann.m"];
Get["Bfield+Drift/Bfields.m"]
Get["Bfield+Drift/Drifts.m"];
Get["Common/CommonFunctions.m"];
Get["NeutronBeam/NeutronBeamDefinition.m"];
Get["Aperture/BooleApert.m"];
Get["Spectra/Backscattering.m"];


(*SetSystemOptions["CacheOptions" -> "Numeric" -> "Cache" -> False];
SetSystemOptions["CacheOptions" -> "Symbolic" -> "Cache" -> False];
$HistoryLength=0;*)

charge= 1;   

BackScatterBoole = True;
TofP[p_] := p^2/(2*mp);
pacc[p_, U_] := Sqrt[2*mp*(p^2/(2*mp) - U)];
thDetAcc[p_, th0_, rD_, U_] := 
 ArcSin[Sqrt[rD*p^2*Sin[th0]^2/(2*mp*(p^2/(2*mp) - U))]];
fitFunc[en_, 
  th_] := (1 - 0.859228 E^(-0.097232 th^2)) (2.53741 - 0.0800777 en)
   
xGCAShift[xDV_, phiDV_, p_, th0_, BRxB_, rRxB_, rA_, XAShift_] := (xDV - rG[p, th0, BRxB/rRxB]*Cos[phiDV])*Sqrt[1/rA] + XAShift
yGCAShift[yDV_, phiDV_, p_, th0_, BRxB_, rRxB_, rA_, YAShift_] := (yDV - rG[p, th0, BRxB/rRxB]*Sin[phiDV])*Sqrt[1/rA] + YAShift
yRxBGCShift[yDV_, phiDV_, p_, th0_, BRxB_, rRxB_, YRxBShift_] := (yDV - rG[p, th0, BRxB/rRxB]*Sin[phiDV]) *Sqrt[1/rRxB] + YRxBShift  
  
DeltayDVShift[phiDV_, yD_, p_,alpha_, th0_, BRxB_, rRxB_, rD_, phiDet_, YDetShift_] := 
(	
	(yD - rG[p, th0, rD*BRxB/rRxB]*Sin[phiDet] - YDetShift)*Sqrt[rD]
) + rG[p, th0, BRxB/rRxB]*Sin[phiDV];

DeltaxDVShift[phiDV_, yDV_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rD_, phiDet_, G1_, G2_, YRxBShift_, XDetShift_] := 

(
	(xD - rG[p, th0, rD*BRxB/rRxB]*Cos[phiDet] - XDetShift)*Sqrt[rD/rRxB](*error here*) -
	charge*D1stBPolyGradwoR[p, alpha, th0, BRxB, rRxB, yRxBGCShift[yDV, phiDV, p, th0, BRxB, rRxB, YRxBShift], G1, G2,YRxBShift]
   	(*DxBPolyGradRx[p, alpha, th0, BRxB, rRxB, yRxBGCShift[yDV, phiDV, p, th0, BRxB, rRxB, YRxBShift], R, G1, G2, charge, YRxBShift]*)
)*Sqrt[rRxB] + rG[p, th0, BRxB/rRxB]*Cos[phiDV];
   
   
DeltaxDV2Shift[phiDV_, yD_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rD_, phiDet_, G1_, G2_, {YRxBShift_, XDetShift_, YDetShift_}] := 
 	DeltaxDVShift[phiDV, 
  		DeltayDVShift[phiDV, yD, p, alpha, th0, BRxB, rRxB, rD, phiDet, YDetShift],
  		xD, p, th0, alpha, BRxB, rRxB, rD, phiDet, G1, G2, YRxBShift, XDetShift
  	]

xGCAShift2[yD_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, G1_, G2_, {XAShift_, YRxBShift_, XDetShift_, YDetShift_}] :=
	xGCAShift[
 		DeltaxDV2Shift[phiDV, yD, xD, p, th0, alpha, BRxB, rRxB, rD, phiDet, G1, G2, {YRxBShift, XDetShift, YDetShift}], 
 		phiDV, p, th0, BRxB, rRxB, rA, XAShift
 	]
 	
yGCAShift2[yD_, p_, alpha_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, YAShift_, YDetShift_] :=
	yGCAShift[
 		DeltayDVShift[phiDV, yD, p, alpha, th0, BRxB, rRxB, rD, phiDet, YDetShift], 
 		phiDV, p, th0, BRxB, rRxB, rA, YAShift
 	]
   
xAShift[phiA_, yD_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, G1_, G2_, {XAShift_, YRxBShift_, XDetShift_, YDetShift_}] := 
 	xGCAShift[
 		DeltaxDV2Shift[phiDV, yD, xD, p, th0, alpha, BRxB, rRxB, rD, phiDet, G1, G2, {YRxBShift, XDetShift, YDetShift}], 
 		phiDV, p, th0, BRxB, rRxB, rA, XAShift
 	] + rG[p, th0, rA*BRxB/rRxB]*Cos[phiA]
  
yAShift[phiA_, yD_, p_, alpha_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, YAShift_, YDetShift_] := 
 	yGCAShift[
 		DeltayDVShift[phiDV, yD, p, alpha, th0, BRxB, rRxB, rD, phiDet, YDetShift], 
 		phiDV, p, th0, BRxB, rRxB, rA, YAShift
 	] + rG[p, th0, rA*BRxB/rRxB]*Sin[phiA]


(*Compiled Function*)
DeltaxDV2ShiftSet[phiDV_, yD_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rD_, phiDet_, G1_, G2_, {YRxBShift_, XDetShift_, YDetShift_}] = 
 	DeltaxDVShift[phiDV, 
  		DeltayDVShift[phiDV, yD, p, alpha, th0, BRxB, rRxB, rD, phiDet, YDetShift],
  		xD, p, th0, alpha, BRxB, rRxB, rD, phiDet, G1, G2, YRxBShift, XDetShift
  	]

DeltaxDV2ShiftCompiled =
Compile[
	{
		{phiDV,_Real}, {yD,_Real}, {xD,_Real}, {p,_Real}, {th0,_Real}, {alpha,_Real}, {BRxB,_Real}, {rRxB,_Real}, {rD,_Real},
		{phiDet,_Real}, {G1,_Real}, {G2,_Real}, {YRxBShift,_Real}, {XDetShift,_Real}, {YDetShift,_Real}
	},
	DeltaxDV2ShiftSet[phiDV, yD, xD, p, th0, alpha, BRxB, rRxB, rD, phiDet, G1, G2, {YRxBShift, XDetShift, YDetShift}],
  CompilationTarget -> "C"(*, RuntimeAttributes -> {Listable}*), RuntimeOptions -> "Speed"
];
   
   
   

(*Compiled Integrand*)

If[BackScatterBoole,
	
	Integrand2DwNBeamCompiledShiftRxBS[
	a_, yD_, xD_, {phiDV_?NumericQ, phiA_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, {alpha_, BRxB_, rRxB_, rA_,rD_, G1_, G2_},
	{XAShift_, YAShift_, YRxBShift_, XDetShift_, YDetShift_}, 
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] :=
 
 pmomNormed[p,a]*Sin[th0]*(1-fitFunc[TofP[pacc[p, -15000]]/1000, thDetAcc[p, th0, rD, -15000]*180/Pi]/100)*(*Change the backscattering description correctly*)
  
  ApertBooleCompiled[
   xAShift[phiA, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, {XAShift, YRxBShift, XDetShift, YDetShift}], 
   yAShift[phiA, yD, p, alpha, th0, BRxB, rRxB, rA, rD, phiDet, YAShift, YDetShift], 
   	xAA, yAA, xOff, yOff],
   	
   	
   	Integrand2DwNBeamCompiledShiftRxBS[
	a_, yD_, xD_, {phiDV_?NumericQ, phiA_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, {alpha_, BRxB_, rRxB_, rA_,rD_, G1_, G2_},
	{XAShift_, YAShift_, YRxBShift_, XDetShift_, YDetShift_}, 
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] :=
 
 pmomNormed[p,a]*Sin[th0]*
    
  ApertBooleCompiled[
   xAShift[phiA, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, {XAShift, YRxBShift, XDetShift, YDetShift}], 
   yAShift[phiA, yD, p, alpha, th0, BRxB, rRxB, rA, rD, phiDet, YAShift, YDetShift], 
   	xAA, yAA, xOff, yOff]
]  



(* ::Section:: *)
(* Integrand with phiA limits and compiled integrand *)


phiAminApertShift[yD_, xD_, plocal_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, G1_, G2_, xOff_, xAA_,{XAShift_, YRxBShift_, XDetShift_, YDetShift_}] := 
	ArcCos[(xOff + xAA/2 - 
		xGCAShift2[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, {XAShift, YRxBShift, XDetShift, YDetShift}]
		)/rG[plocal, th0, rA*BRxB/rRxB]]
  (*Solve[xOff + xAA/2 == xAShift[phiAlocal, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, {XAShift, YRxBShift, XDetShift, YDetShift}], phiAlocal][[All, 1, 2]];*)
phiAmaxApertShift[yD_, xD_, plocal_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, G1_, G2_, xOff_, xAA_,{XAShift_, YRxBShift_, XDetShift_, YDetShift_}] := 
   ArcCos[(
   	xOff - xAA/2 - 
   	xGCAShift2[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, {XAShift, YRxBShift, XDetShift, YDetShift}]
   	)/rG[plocal, th0, rA*BRxB/rRxB]]
   (*Solve[xOff - xAA/2 == xAShift[phiAlocal, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, {XAShift, YRxBShift, XDetShift, YDetShift}], phiAlocal][[All, 1, 2]];*)


(*phiACases[yD_?NumericQ, xD_?NumericQ, plocal_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_] :=
 {phiA,
 	Sequence@@Re[{
  phiAmaxApert[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[1]],
  phiAminApert[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[1]],
  phiAminApert[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[2]],
  phiAmaxApert[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[2]]
  }]
 }*)



(*Integrand with manual summation of phiA integration*)
(*ManualphiAIntegrand[b_, yD_, xD_, {phiDV_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, 
	{alpha_, BRxB_, rRxB_, rA_, rD_, R_, G1_, G2_}, {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] := 
	Module[
  		{phimin2, phimax2, phimean, int5Dvalue},
  		phimin2 = Re[phiAminApert[yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA]];
 		phimax2 = Re[phiAmaxApert[yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA]];
  		phimean = (phimax2[[1]] + phimin2[[1]])/2;
  		int5Dvalue = Integrand2DwNBeamCompiled[b, yD, xD, {phiDV, phimean, phiDet, p, th0}, {alpha, BRxB, rRxB, rA, rD, 
     R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, 
     k3y}, {xAA, yAA, xOff, yOff}];
  		(phimin2[[1]] - phimax2[[1]])* int5Dvalue + (phimax2[[2]] - phimin2[[2]])*int5Dvalue
  ]*)
  

(*phiA limits from Y Aperture*)

phiAminApertfromYShift[yD_, plocal_, alpha_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_, YAShift_, YDetShift_] := 
	ArcSin[(yOff + yAA/2 - yGCAShift2[yD, plocal, alpha, th0, BRxB, rRxB, rA, rD, phiDet, YAShift, YDetShift])/rG[plocal, th0, rA*BRxB/rRxB]]
 (*Solve[yOff + yAA/2 == yAShift[phiAlocal, yD, plocal, alpha, th0, BRxB, rRxB, rA, rD, phiDet, YAShift, YDetShift, R], phiAlocal][[All, 1, 2]]*)
 
phiAmaxApertfromYShift[yD_, plocal_, alpha_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_, YAShift_, YDetShift_] := 
 ArcSin[(yOff - yAA/2 - yGCAShift2[yD, plocal, alpha, th0, BRxB, rRxB, rA, rD, phiDet, YAShift, YDetShift])/rG[plocal, th0, rA*BRxB/rRxB]]
 (*Solve[yOff - yAA/2 == yAShift[phiAlocal, yD, plocal, alpha, th0, BRxB, rRxB, rA, rD, phiDet, YAShift, YDetShift, R], phiAlocal][[All, 1, 2]]*)

  
  
  
  
(*new manual phiA integration with all phi limits*)
ManualphiAIntegrandAllLimitsShift[
	a_, yD_, xD_, {phiDV_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, 
	{alpha_, BRxB_, rRxB_, rA_, rD_, G1_, G2_}, {XAShift_, YAShift_, YRxBShift_, XDetShift_, YDetShift_}, 
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] := 
	Block[
  		{
  			phimin = Re[phiAminApertShift[yD, xD, If[p==0.,0.00000001,p], th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}]],
  			phimax = Re[phiAmaxApertShift[yD, xD, If[p==0.,0.00000001,p], th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}]], 
  			phiminY = Re[phiAminApertfromYShift[yD, If[p==0.,0.00000001,p], alpha, th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA, YAShift, YDetShift]], 
  			phimaxY = Re[phiAmaxApertfromYShift[yD, If[p==0.,0.00000001,p], alpha, th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA, YAShift, YDetShift]],
			IntValues, philimitlistSorted, Deltalimits, phiCenters, Integrated
		},
  		
  		philimitlistSorted = (*Flatten[*){-Pi//N, phimin, phimax,-phimin, -phimax(*also minus because ArcCos and negative both true*), phiminY, phimaxY, Pi//N}(*]*);
		philimitlistSorted = If[# > Pi, # - 2 Pi, #] & /@ philimitlistSorted;
		philimitlistSorted = Union[philimitlistSorted];
		Deltalimits = philimitlistSorted[[2;;-1]] - philimitlistSorted[[1;;-2]];
		phiCenters = philimitlistSorted[[1;;-2]] + Deltalimits/2;
		IntValues = Integrand2DwNBeamCompiledShiftRxBS[a, yD, xD, {phiDV, #, phiDet, p, th0}, {alpha, BRxB, rRxB, rA, rD, G1,G2}, 
			{XAShift, YAShift, YRxBShift, XDetShift, YDetShift},
			{twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA, yAA, xOff, yOff}] & /@ phiCenters;
		Integrated = Total[Deltalimits*IntValues];
		(*ClearAll[phimin,phimax,phiminY,phimaxY,IntValues,philimitlistSorted,Deltalimits,phiCenters];*)
  		Integrated
  		
  ]  




(* ::Section:: *)
(* Integrand with p limits from aperture *)


pminApertShift[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, G1_, G2_, xOff_, xAA_,
	{XAShift_, YRxBShift_, XDetShift_, YDetShift_}] :=
 	NSolve[xOff + xAA/2 == 
 		xAShift[phiA, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1,G2, {XAShift, YRxBShift, XDetShift, YDetShift}],
 		plocal,Reals][[All, 1, 2]]
      
      
pminApertCasesShift[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, G1_, G2_, xOff_, xAA_,
	{XAShift_, YRxBShift_, XDetShift_, YDetShift_}] := 
Block[
  {
  	pmin = Cases[
  	pminApertShift[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}],
  	 _?( -2*ppmax <= # <= 2*ppmax &)]
  	 },
  Piecewise[
   {
    {Print["pmin >2"], Length[pmin] > 2},
    {Which[#<0,0.,#>ppmax,ppmax,True,#]&/@pmin, Length[pmin] <= 2}(*,
    {
    	Print[
    		"pmin: Length > 2",
    		Length[pmin],
    		pmin,phiA," ",th0," ",phiDet," ",
    		pminApertShift[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}]
    	], 
    	Length[pmin]>2(*Length[pmin]!=1 && Length[pmin]!=0*) 
    }*)
    },
    Print["pmin: other Case",pminApertShift[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}]]
   ]
]
  
  
pmaxApertShift[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, G1_, G2_, xOff_, xAA_,
	{XAShift_,YRxBShift_, XDetShift_, YDetShift_}] := 
 		NSolve[xOff - xAA/2 == xAShift[phiA, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1,G2, {XAShift, YRxBShift, XDetShift, YDetShift}], plocal,
 			Reals][[All, 1, 2]]


pmaxApertCasesShift[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, G1_, G2_, xOff_, xAA_,
	{XAShift_, YRxBShift_, XDetShift_, YDetShift_}] := 
Block[
  {
  	ppmaxx =Cases[
  		pmaxApertShift[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}],
  		 _?(-2*ppmax <= # <= 2*ppmax &)]},
  
  Piecewise[
   {
   	{Print["ppmax >2"], Length[ppmaxx] > 2},
    {Which[#<0,0.,#>ppmax,ppmax,True,#]&/@ppmaxx, Length[ppmaxx] <= 2}(*,
    {	
    	Print["ppmax: Length > 2",Length[ppmaxx],ppmaxx,phiA," ",th0," ",phiDet," ",
    		ppmaxApertShift[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}]
    	],
    	Length[ppmaxx]>2(*Length[ppmaxx]!=1 && Length[ppmaxx]!=0*) 
    }*)
    },
    Print["ppmax: other case"]
   ]

  ]
  

pLimitsApertXListShift[yD_?NumericQ, xD_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, G1_, G2_, xOff_, xAA_,
	{XAShift_, YRxBShift_, XDetShift_, YDetShift_}] := 
(*pLimitsApertXList[yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA] =*)
Block[
  {
   pminXmPi = pminApertCasesShift[-Pi, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}],
   pminX0 = pminApertCasesShift[0., yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}],
   ppmaxXmPi = pmaxApertCasesShift[-Pi, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}],
   ppmaxX0 = pmaxApertCasesShift[0., yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}],
   plimitlist
   },
	plimitlist=Union[Flatten[{pminXmPi, pminX0, ppmaxXmPi, ppmaxX0}]];
	Which[
		Length[plimitlist]==1,Join[plimitlist,plimitlist],
		True,plimitlist
	]
   
]
  


(*p limits from apertY*)

(*pminApertY[phiA_, yD_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_] =
 Module[{plocal}, Solve[yOff + yAA/2 == yA[phiA, yD, plocal, th0, BRxB, rRxB, rA, rD, phiDet], plocal, Reals][[1, 1, 2, 1]]]
 
ppmaxApertY[phiA_, yD_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_] =
Module[{plocal}, Solve[yOff - yAA/2 == yA[phiA, yD, plocal, th0, BRxB, rRxB, rA, rD, phiDet], plocal][[1, 1, 2]]]

pminApertYCases[yD_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_] := Module[
  {phiAlocal = If[yD > yAA/2 + yOff, -Pi/2, Pi/2], pminlocal},
  pminlocal = pminApertY[phiAlocal, yD, th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA];
  Piecewise[{{pminlocal, 0. < pminlocal < ppmax}}, 0.]
  ]

ppmaxApertYCases[yD_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_] := Module[
  {phiAlocal = If[yD < -yAA/2 + yOff, Pi/2, -Pi/2], ppmaxlocal},
  ppmaxlocal = ppmaxApertY[phiAlocal, yD, th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA];
  Piecewise[{{ppmaxlocal, 0. < ppmaxlocal < ppmax}}, 0.]
  ]

pLimitsApertAllList[yD_?NumericQ, xD_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, {xAA_, xOff_,yAA_, yOff_}] := 
Module[
  {
   pminXmPi = pminApertCases[-Pi, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA],
   pminX0 = pminApertCases[0., yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA],
   ppmaxXmPi = ppmaxApertCases[-Pi, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA],
   ppmaxX0 = ppmaxApertCases[0., yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA],
   pminY = pminApertYCases[yD, th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA],
   ppmaxY = ppmaxApertYCases[yD, th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA]
   },
 	Sort[{pminXmPi, pminX0, ppmaxXmPi, ppmaxX0,pminY,ppmaxY}]
  ]*)



(* ::Section:: *)
(* 2 step-Integration *)


pXDomainShift[yD_?NumericQ, xD_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, G1_, G2_, xOff_, xAA_,
	{XAShift_, YRxBShift_, XDetShift_, YDetShift_}]:=
{p,Sequence@@pLimitsApertXListShift[yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}]}


Step1IntShift[a_, yD_, xD_, {phiDet_?NumericQ, th0_?NumericQ}, 
	{alpha_, BRxB_, rRxB_, rA_, rD_, G1_, G2_},{XAShift_, YAShift_, YRxBShift_, XDetShift_, YDetShift_}, 
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},
	method1_,PrecGoal1_,AccGoal1_,MinRec1_,MaxRec1_]:=
	NIntegrate[
		ManualphiAIntegrandAllLimitsShift[a, yD, xD, {phiDV, phiDet, p, th0}, 
			{alpha, BRxB, rRxB, rA, rD, G1, G2}, {XAShift, YAShift, YRxBShift, XDetShift, YDetShift}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA, yAA, xOff, yOff}],
		Evaluate[pXDomainShift[yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}]],
		{phiDV,-Pi,Pi},
		Method->method1,PrecisionGoal->PrecGoal1,AccuracyGoal->AccGoal1, MinRecursion->MinRec1, MaxRecursion->MaxRec1	
]


Step2IntShift[a_, xD_?NumericQ,yD_?NumericQ, {alpha_, BRxB_, rF_, rRxB_, rA_, rD_, G1_, G2_},{XAShift_, YAShift_, YRxBShift_, XDetShift_, YDetShift_}, 
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_},
	{xAA_, yAA_, xOff_, yOff_}, method1_,PrecGoal1_,AccGoal1_,MinRec1_,MaxRec1_, method2_,PrecGoal2_,AccGoal2_,MinRec2_,MaxRec2_]:=

		NIntegrate[
		Step1IntShift[a, yD, xD, {phiDet, th0}, 
			{alpha, BRxB, rRxB, rA, rD, G1, G2}, {XAShift, YAShift, YRxBShift, XDetShift, YDetShift}, 
			{twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA, yAA, xOff, yOff},
			method1,PrecGoal1,AccGoal1,MinRec1,MaxRec1],
			{th0,0.,thetamax[rF]},{phiDet,-Pi,Pi},Method->method2,PrecisionGoal->PrecGoal2,AccuracyGoal->AccGoal2,MinRecursion->MinRec2, MaxRecursion->MaxRec2		
		]
		
		


(* ::Section:: *)
(* Add xD and yD integration *)


BinIntShift[a_, OneBinList_, {alpha_, BRxB_, rF_, rRxB_, rA_, rD_, G1_, G2_},{XAShift_, YAShift_, YRxBShift_, XDetShift_, YDetShift_}, 
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_},
	{xAA_, yAA_, xOff_, yOff_}, 
	{method1_,PrecGoal1_,AccGoal1_,MinRec1_,MaxRec1_}, 
	{method2_,PrecGoal2_,AccGoal2_,MinRec2_,MaxRec2_}, 
	{method3_,PrecGoal3_,AccGoal3_,MinRec3_,MaxRec3_}]:=
	NIntegrate[
		
		Step2IntShift[a, xD, yD, {alpha, BRxB, rF, rRxB, rA, rD, G1, G2}, {XAShift, YAShift, YRxBShift, XDetShift, YDetShift}, 
			{twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA, yAA, xOff, yOff},
			method1,PrecGoal1,AccGoal1,MinRec1,MaxRec1,method2,PrecGoal2,AccGoal2,MinRec2,MaxRec2],
			
		{xD, OneBinList[[1,1]], OneBinList[[1,2]]}, {yD, OneBinList[[2,1]], OneBinList[[2,2]]},
		PrecisionGoal->PrecGoal3,Method->method3,MinRecursion->MinRec3,MaxRecursion->MaxRec3,AccuracyGoal->AccGoal3
	]


bin2DGen[xStart_, xEnd_, yStart_, yEnd_, binNoX_, binNoY_] := 
 Block[{binLengthX, binLengthY, binList, binsX, binsY},
  binLengthX = Round[(xEnd - xStart)/binNoX,10^-20.];
  binLengthY = Round[(yEnd - yStart)/binNoY,10^-20.];
  binsX = 
   Table[{xStart + (i - 1)*binLengthX, xStart + i*binLengthX}, {i, 
     binNoX}];
  binsY = 
   Table[{yStart + (i - 1)*binLengthY, yStart + i*binLengthY}, {i, 
     binNoY}];
  binList = 
   Flatten[Table[{binsX[[i]], binsY[[j]]}, {i, Length[binsX]}, {j, 
      Length[binsY]}], 1];
  Return[binList];
  ]
  
bin2DGen2[xStart_, xEnd_, yStart_, yEnd_, binLengthX_, binNoY_] := 
 Block[{binNoX, binLengthY, binList, binsX, binsY}, 
  binNoX = Round[(xEnd - xStart)/binLengthX];
  binLengthY = Round[(yEnd - yStart)/binNoY, 10^-20.];
  binsX = 
   Table[{xStart + (i - 1)*binLengthX, xStart + i*binLengthX}, {i, 
     binNoX}];
  binsY = 
   Table[{yStart + (i - 1)*binLengthY, yStart + i*binLengthY}, {i, 
     binNoY}];
  binList = 
   Flatten[Table[{binsX[[i]], binsY[[j]]}, {i, Length[binsX]}, {j, 
      Length[binsY]}], 1];
  Return[binList];]
