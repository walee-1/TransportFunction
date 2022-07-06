(* ::Package:: *)

(* Wolfram Language package *)

Get["Spectra/ProtonSpectrumNachtmann.m"];
Get["Bfield+Drift/Bfields.m"]
Get["Bfield+Drift/Drifts.m"];
Get["Common/CommonFunctions.m"];
Get["Aperture/BooleApert.m"];
Get["Spectra/Backscattering.m"];


(*SetSystemOptions["CacheOptions" -> "Numeric" -> "Cache" -> False];
SetSystemOptions["CacheOptions" -> "Symbolic" -> "Cache" -> False];
$HistoryLength=0;*)

charge= 1;   

BackScatterBoole = False;
   


yRxBGCShift[yDet_, phiDet_, p_, th0_, BRxB_, rRxB_, rD_, yDetShift_] := (yDet - rG[p, th0, rD*BRxB/rRxB]*Sin[phiDet]-yDetShift) *Sqrt[rD/rRxB] 
  
DeltayAShift[phiA_, yD_, p_,alpha_, th0_, BRxB_, rRxB_, rD_, rA_,phiDet_, YDetShift_] := 
(	
	(yD - rG[p, th0, rD*BRxB/rRxB]*Sin[phiDet] - YDetShift)*Sqrt[rD/rA](*-
	DYSimple[p,alpha,BRxB,th0,rRxB,R]*)
)(**Sqrt[rRxB]*) + rG[p, th0, rA*BRxB/rRxB]*Sin[phiA];

DeltaxAShift[phiA_, yA_, xD_, yD_, p_, th0_, alpha_, BRxB_, rRxB_, rD_, rA_, phiDet_, G1_, G2_, yDetShift_, YRxBShift_,XDetShift_] := 

(
	(xD - rG[p, th0, rD*BRxB/rRxB]*Cos[phiDet] - XDetShift)*Sqrt[rD/rRxB] -
	charge*D1stBPolyGradwoR[p, alpha, th0, BRxB, rRxB, yRxBGCShift[yD, phiDet, p, th0, BRxB, rRxB, rD,yDetShift], G1, G2,YRxBShift]
   	(*DxBPolyGradRx[p, alpha, th0, BRxB, rRxB, yRxBGCShift[yDV, phiDV, p, th0, BRxB, rRxB, YRxBShift], R, G1, G2, charge, YRxBShift]*)
)*Sqrt[rRxB/rA] + rG[p, th0, rA*BRxB/rRxB]*Cos[phiA];
   
   

   
xAShift[phiA_, yA_, xD_, yD_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, G1_, G2_, {YRxBShift_, XDetShift_, YDetShift_}] := 
 	
 		DeltaxAShift[phiA, yA, xD, yD, p, th0, alpha, BRxB, rRxB, rD, rA, phiDet, G1, G2, YDetShift, YRxBShift, XDetShift]
  
yAShift[phiA_, yD_, p_, alpha_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, YAShift_, YDetShift_] := 
 		DeltayAShift[phiA, yD, p, alpha, th0, BRxB, rRxB, rD, rA, phiDet, YDetShift]


   
   
   

(*Compiled Integrand*)

If[BackScatterBoole,
	
	Integrand2DwNBeamCompiledShiftRxBS[
	a_, yD_, xD_, {phiA_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, {alpha_, BRxB_, rRxB_, rA_,rD_, G1_, G2_},
	{XAShift_, YAShift_, YRxBShift_, XDetShift_, YDetShift_}, {xAA_, yAA_, xOff_, yOff_}] :=
 
 pmomNormed[p,a]*Sin[th0]*(1-bsFuncElectron2Compiled[p, th0, rD])*(*Change the backscattering description correctly*)
  
  ApertBooleCompiled[
   xAShift[phiA, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, {XAShift, YRxBShift, XDetShift, YDetShift}], 
   yAShift[phiA, yD, p, alpha, th0, BRxB, rRxB, rA, rD, phiDet, YAShift, YDetShift], 
   	xAA, yAA, xOff, yOff],
   	
   	
   	Integrand2DwNBeamCompiledShiftRxBS[
	a_, yD_, xD_, {phiA_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, {alpha_, BRxB_, rRxB_, rA_,rD_, G1_, G2_},
	{XAShift_, YAShift_, YRxBShift_, XDetShift_, YDetShift_}, {xAA_, yAA_, xOff_, yOff_}] :=
 
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
		(xD - rG[plocal, th0, rD*BRxB/rRxB]*Cos[phiDet] - XDetShift)*Sqrt[rD/rRxB] -
	charge*D1stBPolyGradwoR[p, alpha, th0, BRxB, rRxB, yRxBGCShift[yD, phiDet, plocal, th0, BRxB, rRxB, rD,YDetShift], G1, G2,YRxBShift]
   	(*DxBPolyGradRx[p, alpha, th0, BRxB, rRxB, yRxBGCShift[yDV, phiDV, p, th0, BRxB, rRxB, YRxBShift], R, G1, G2, charge, YRxBShift]*)
)*Sqrt[rRxB/rA]/rG[plocal, th0, rA*BRxB/rRxB]]
  (*Solve[xOff + xAA/2 == xAShift[phiAlocal, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, {XAShift, YRxBShift, XDetShift, YDetShift}], phiAlocal][[All, 1, 2]];*)
phiAmaxApertShift[yD_, xD_, plocal_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, G1_, G2_, xOff_, xAA_,{XAShift_, YRxBShift_, XDetShift_, YDetShift_}] := 
   ArcCos[(
   	(xD - rG[plocal, th0, rD*BRxB/rRxB]*Cos[phiDet] - XDetShift)*Sqrt[rD/rRxB] -
	charge*D1stBPolyGradwoR[p, alpha, th0, BRxB, rRxB, yRxBGCShift[yD, phiDet, plocal, th0, BRxB, rRxB, rD,YDetShift], G1, G2,YRxBShift]
   	(*DxBPolyGradRx[p, alpha, th0, BRxB, rRxB, yRxBGCShift[yDV, phiDV, p, th0, BRxB, rRxB, YRxBShift], R, G1, G2, charge, YRxBShift]*)
)*Sqrt[rRxB/rA]/rG[plocal, th0, rA*BRxB/rRxB]]
   (*Solve[xOff - xAA/2 == xAShift[phiAlocal, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, {XAShift, YRxBShift, XDetShift, YDetShift}], phiAlocal][[All, 1, 2]];*)


(*phiA limits from Y Aperture*)

phiAminApertfromYShift[yD_, plocal_, alpha_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_, YAShift_, YDetShift_] := 
	ArcSin[(yOff + yAA/2 - (yD - rG[plocal, th0, rD*BRxB/rRxB]*Sin[phiDet] - YDetShift)*Sqrt[rD/rA])/rG[plocal, th0, rA*BRxB/rRxB]]
 (*Solve[yOff + yAA/2 == yAShift[phiAlocal, yD, plocal, alpha, th0, BRxB, rRxB, rA, rD, phiDet, YAShift, YDetShift, R], phiAlocal][[All, 1, 2]]*)
 
phiAmaxApertfromYShift[yD_, plocal_, alpha_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_, YAShift_, YDetShift_] := 
 ArcSin[(yOff - yAA/2 - (yD - rG[plocal, th0, rD*BRxB/rRxB]*Sin[phiDet] - YDetShift)*Sqrt[rD/rA])/rG[plocal, th0, rA*BRxB/rRxB]]
 (*Solve[yOff - yAA/2 == yAShift[phiAlocal, yD, plocal, alpha, th0, BRxB, rRxB, rA, rD, phiDet, YAShift, YDetShift, R], phiAlocal][[All, 1, 2]]*)

  
  
  
  
(*new manual phiA integration with all phi limits*)
ManualphiAIntegrandAllLimitsShift[
	a_, yD_, xD_, {phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, 
	{alpha_, BRxB_, rRxB_, rA_, rD_, G1_, G2_}, {XAShift_, YAShift_, YRxBShift_, XDetShift_, YDetShift_}, {xAA_, yAA_, xOff_, yOff_}] := 
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
		IntValues = Integrand2DwNBeamCompiledShiftRxBS[a, yD, xD, { #, phiDet, p, th0}, {alpha, BRxB, rRxB, rA, rD, G1,G2}, 
			{XAShift, YAShift, YRxBShift, XDetShift, YDetShift},{xAA, yAA, xOff, yOff}] & /@ phiCenters;
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
  



(* ::Section:: *)
(* 2 step-Integration *)


pXDomainShift[yD_?NumericQ, xD_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, G1_, G2_, xOff_, xAA_,
	{XAShift_, YRxBShift_, XDetShift_, YDetShift_}]:=
{p,Sequence@@pLimitsApertXListShift[yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}]}


Step1IntShift[a_, yD_, xD_, {phiDet_?NumericQ, th0_?NumericQ}, 
	{alpha_, BRxB_, rRxB_, rA_, rD_, G1_, G2_},{XAShift_, YAShift_, YRxBShift_, XDetShift_, YDetShift_}, {xAA_, yAA_, xOff_, yOff_},
	method1_,PrecGoal1_,AccGoal1_,MinRec1_,MaxRec1_]:=
		ManualphiAIntegrandAllLimitsShift[a, yD, xD, {phiDet, p, th0}, 
			{alpha, BRxB, rRxB, rA, rD, G1, G2}, {XAShift, YAShift, YRxBShift, XDetShift, YDetShift}, {xAA, yAA, xOff, yOff}],
		Evaluate[pXDomainShift[yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}]]


Step2IntShift[a_, xD_?NumericQ,yD_?NumericQ, {alpha_, BRxB_, rF_, rRxB_, rA_, rD_, G1_, G2_},{XAShift_, YAShift_, YRxBShift_, XDetShift_, YDetShift_},
	{xAA_, yAA_, xOff_, yOff_}, method1_,PrecGoal1_,AccGoal1_,MinRec1_,MaxRec1_, method2_,PrecGoal2_,AccGoal2_,MinRec2_,MaxRec2_]:=

		NIntegrate[
		Step1IntShift[a, yD, xD, {phiDet, th0}, 
			{alpha, BRxB, rRxB, rA, rD, G1, G2}, {XAShift, YAShift, YRxBShift, XDetShift, YDetShift}, {xAA, yAA, xOff, yOff},
			method1,PrecGoal1,AccGoal1,MinRec1,MaxRec1],
			{th0,0.,thetamax[rF]},{phiDet,-Pi,Pi},Method->method2,PrecisionGoal->PrecGoal2,AccuracyGoal->AccGoal2,MinRecursion->MinRec2, MaxRecursion->MaxRec2		
		]
		
		


(* ::Section:: *)
(* Add xD and yD integration *)


BinIntShift[a_, OneBinList_, {alpha_, BRxB_, rF_, rRxB_, rA_, rD_, G1_, G2_},{XAShift_, YAShift_, YRxBShift_, XDetShift_, YDetShift_},
	{xAA_, yAA_, xOff_, yOff_}, 
	{method1_,PrecGoal1_,AccGoal1_,MinRec1_,MaxRec1_}, 
	{method2_,PrecGoal2_,AccGoal2_,MinRec2_,MaxRec2_}, 
	{method3_,PrecGoal3_,AccGoal3_,MinRec3_,MaxRec3_}]:=
	NIntegrate[
		
		Step2IntShift[a, xD, yD, {alpha, BRxB, rF, rRxB, rA, rD, G1, G2}, {XAShift, YAShift, YRxBShift, XDetShift, YDetShift}, {xAA, yAA, xOff, yOff},
			method1,PrecGoal1,AccGoal1,MinRec1,MaxRec1,method2,PrecGoal2,AccGoal2,MinRec2,MaxRec2],
			
		{xD, OneBinList[[1,1]], OneBinList[[1,2]]}, {yD, OneBinList[[2,1]], OneBinList[[2,2]]},
		PrecisionGoal->PrecGoal3,Method->method3,MinRecursion->MinRec3,MaxRecursion->MaxRec3,AccuracyGoal->AccGoal3
	]


bin2DGen[xStart_, xEnd_, yStart_, yEnd_, binNoX_, binNoY_] := 
 Block[{binLengthX, binLengthY, binList, binsX, binsY},
  binLengthX = Round[(xEnd - xStart)/binNoX,10^-6.];
  binLengthY = Round[(yEnd - yStart)/binNoY,10^-6.];
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
     
