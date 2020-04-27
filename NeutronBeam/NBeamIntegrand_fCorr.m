(* Wolfram Language package *)
Get["Spectra/ElectronSpectrum.m"];
Get["Bfield+Drift/Bfields.m"]
Get["Bfield+Drift/Drifts.m"];
Get["Common/CommonFunctions.m"];
Get["NeutronBeam/NeutronBeamDefinition.m"];
Get["Aperture/BooleApert.m"];
   
   
fCorr[p_,th0_,th1max_, rRxB_]:= 1 + p^2*theta2[th0,rRxB]^3*th1max^2/pmax^2/thetamax[2.]^3/1.5^2*0.00029  
   
xGCA[xDV_, phiDV_, p_, th0_, BRxB_, rRxB_, 
  rA_] := (xDV - rG[p, th0, BRxB/rRxB]*Cos[phiDV])*Sqrt[1/rA]
yGCA[yDV_, phiDV_, p_, th0_, BRxB_, rRxB_, 
  rA_] := (yDV - rG[p, th0, BRxB/rRxB]*Sin[phiDV])*Sqrt[1/rA]
  
  
DeltayDV[phiDV_, yD_, p_, th0_, BRxB_, rRxB_, rA_, rD_, 
  phiDet_] := (yD - rG[p, th0, rD*BRxB/rRxB]*Sin[phiDet])*
   Sqrt[rD/rRxB]*Sqrt[rA] + rG[p, th0, BRxB/rRxB]*Sin[phiDV]
DeltaxDVfCorr[phiDV_, yDV_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, 
  phiDet_, R_, G1_, 
  G2_,fCorrf_] := ((xD - rG[p, th0, rD*BRxB/rRxB]*Cos[phiDet])*Sqrt[rD/rRxB] -
      D1stBPolyGrad[p, alpha, th0, BRxB, rRxB, 
      yGCA[yDV, phiDV, p, th0, BRxB, rRxB, rA], R, G1, G2]*fCorrf)*Sqrt[rA] +
   rG[p, th0, BRxB/rRxB]*Cos[phiDV]
   
   
DeltaxDV2fCorr[phiDV_, yD_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, 
  phiDet_, R_, G1_, G2_,fCorrf_] := 
 DeltaxDVfCorr[phiDV, 
  DeltayDV[phiDV, yD, p, th0, BRxB, rRxB, rA, rD, phiDet], xD, p, th0,
   alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2,fCorrf]
   
   
xAfCorr[phiA_, yD_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_,
   R_, G1_, G2_,fCorrf_] := 
 xGCA[DeltaxDV2fCorr[phiDV, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, 
    phiDet, R, G1, G2,fCorrf], phiDV, p, th0, BRxB, rRxB, rA] + 
  rG[p, th0, rA*BRxB/rRxB]*Cos[phiA]
  
yA[phiA_, yD_, p_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_] := 
 yGCA[DeltayDV[phiDV, yD, p, th0, BRxB, rRxB, rA, rD, phiDet], phiDV, 
   p, th0, BRxB, rRxB, rA] + rG[p, th0, rA*BRxB/rRxB]*Sin[phiA]
   

(*Compiled Integrand*)

Integrand2DwNBeamCompiledfCorr[
	b_, yD_, xD_, {phiDV_?NumericQ, phiA_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, {alpha_, BRxB_, rRxB_, rA_,rD_, R_, G1_, G2_,fCorrf_}, 
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] :=
 
 wmomNormedWb[\[Lambda]0, \[Kappa]0, b, p]*Sin[th0]*
  
  TrapezNBeamCompiledNormed[
   DeltaxDV2fCorr[phiDV, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2,fCorrf]-xOff, {twx, plx, k1x, k2x, k3x}]*
  TrapezNBeamCompiledNormed[
   DeltayDV[phiDV, yD, p, th0, BRxB, rRxB, rA, rD, phiDet]-yOff, {twy, ply,k1y, k2y, k3y}]*
  
  ApertBooleCompiled[
   xAfCorr[phiA, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2,fCorrf], yA[phiA, yD, p, th0, BRxB, rRxB, rA, rD, phiDet], 
   	xAA, yAA, xOff, yOff]  



(* ::Section:: *)
(* Integrand with phiA limits and compiled integrand *)




phiAminApertfCorr[yD_, xD_, plocal_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_,fCorrf_] = 
  Solve[xOff + xAA/2 == xAfCorr[phiAlocal, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2,fCorrf], phiAlocal][[All, 1, 2, 1, 1]];
phiAmaxApertfCorr[yD_, xD_, plocal_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_,fCorrf_] = 
   Solve[xOff - xAA/2 == xAfCorr[phiAlocal, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2,fCorrf], phiAlocal][[All, 1, 2, 1, 1]];


phiACasesfCorr[yD_?NumericQ, xD_?NumericQ, plocal_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_,fCorrf_] :=
 {phiA,
 	Sequence@@Re[{
  phiAmaxApertfCorr[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA,fCorrf][[1]],
  phiAminApertfCorr[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA,fCorrf][[1]],
  phiAminApertfCorr[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA,fCorrf][[2]],
  phiAmaxApertfCorr[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA,fCorrf][[2]]
  }]
 }



(*Integrand with manual summation of phiA integration*)
ManualphiAIntegrandfCorr[b_, yD_, xD_, {phiDV_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, 
	{alpha_, BRxB_, rRxB_, rA_, rD_, R_, G1_, G2_,fCorrf_}, {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] := 
	Module[
  		{phimin2, phimax2, phimean, int5Dvalue},
  		phimin2 = Re[phiAminApertfCorr[yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA,fCorrf]];
 		phimax2 = Re[phiAmaxApertfCorr[yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA,fCorrf]];
  		phimean = (phimax2[[1]] + phimin2[[1]])/2;
  		int5Dvalue = Integrand2DwNBeamCompiledfCorr[b, yD, xD, {phiDV, phimean, phiDet, p, th0}, {alpha, BRxB, rRxB, rA, rD, 
     R, G1, G2,th1max}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, 
     k3y}, {xAA, yAA, xOff, yOff}];
  		(phimin2[[1]] - phimax2[[1]])* int5Dvalue + (phimax2[[2]] - phimin2[[2]])*int5Dvalue
  ]
  

(*phiA limits from Y Aperture*)

phiAminApertfromY[yD_, plocal_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_] = 
 Solve[yOff + yAA/2 == yA[phiAlocal, yD, plocal, th0, BRxB, rRxB, rA, rD, phiDet], phiAlocal][[All, 1, 2, 1, 1 ;; -2]]
 
phiAmaxApertfromY[yD_, plocal_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_] = 
 Solve[yOff - yAA/2 == yA[phiAlocal, yD, plocal, th0, BRxB, rRxB, rA, rD, phiDet], phiAlocal][[All, 1, 2, 1, 1 ;; -2]]

  
(*new manual phiA integration with all phi limits*)
ManualphiAIntegrandAllLimitsfCorr[
	b_, yD_, xD_, {phiDV_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, 
	{alpha_, BRxB_, rRxB_, rA_, rD_, R_, G1_, G2_,fCorrf_}, {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] := 
	Module[
  		{
  			phimin = Re[phiAminApertfCorr[yD, xD, If[p==0.,0.00000001,p], th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA,fCorrf]],
  			phimax = Re[phiAmaxApertfCorr[yD, xD, If[p==0.,0.00000001,p], th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA,fCorrf]], 
  			phiminY = Re[phiAminApertfromY[yD, If[p==0.,0.00000001,p], th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA]], 
  			phimaxY = Re[phiAmaxApertfromY[yD, If[p==0.,0.00000001,p], th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA]],
			IntValues, philimitlist, philimitlistSorted, Deltalimits, phiCenters, Integrated
		},
  		
  		philimitlist = Flatten[{-Pi//N, phimin, phimax, phiminY, phimaxY, Pi//N}];
		philimitlist = If[# > Pi, # - 2 Pi, #] & /@ philimitlist;
		philimitlistSorted = Union[philimitlist];
		Deltalimits = Table[philimitlistSorted[[i + 1]] - philimitlistSorted[[i]], {i, 1, Length[philimitlistSorted] - 1}];
		phiCenters = Table[philimitlistSorted[[i]] + Deltalimits[[i]]/2, {i, 1, Length[philimitlistSorted] - 1}];
		IntValues = Integrand2DwNBeamCompiledfCorr[b, yD, xD, {phiDV, #, phiDet, p, th0}, {alpha, BRxB, rRxB, rA, rD, R, G1,G2,fCorrf},
			{twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA, yAA, xOff, yOff}] & /@ phiCenters;
		Integrated = Total[Deltalimits*IntValues];
  		Integrated
  		
  ]  




(* ::Section:: *)
(* Integrand with p limits from aperture *)


pminApertfCorr[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_,th1max_] :=
Module[
	{plocal}, 
 	NSolve[xOff + xAA/2 == xAfCorr[phiA, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1,G2,fCorr[plocal,th0,th1max,rRxB]], plocal][[All, 1, 2]]
]
      
      
pminApertCasesfCorr[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_,th1max_] := 
Module[
  {pmin = Cases[pminApertfCorr[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA,th1max], _?(# \[Element] Reals && 0. <= # <= 2*pmax &)]},
  Piecewise[
   {
    {0, Length[pmin] == 0},
    {pmax, Length[pmin] == 1 && pmin[[1]] > pmax},
    {pmin[[1]], Length[pmin] == 1},
    {
    	Print["pmin: Length =",Length[pmin],pmin,phiA," ",th0," ",phiDet," ",pminApertfCorr[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA,th1max]], 
    	Length[pmin]>1(*Length[pmin]!=1 && Length[pmin]!=0*) 
    }
    },
    Print["pmin: other Case"]
   ]
 
]
  
pmaxApertfCorr[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_,th1max_] :=
	Module[
		{plocal}, 
 		NSolve[xOff - xAA/2 == xAfCorr[phiA, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1,G2,fCorr[plocal,th0,th1max,rRxB]], plocal][[All, 1, 2]]
	]


pmaxApertCasesfCorr[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_,th1max_] := 
Module[
  {pmaxx =Cases[pmaxApertfCorr[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA,th1max], _?(# \[Element] Reals && 0. <= # <= 2*pmax &)]},
  
  Piecewise[
   {
    {pmax, Length[pmaxx] == 0},
    {pmax, Length[pmaxx] == 1 && pmaxx[[1]] > pmax},
    {pmaxx[[1]], Length[pmaxx] == 1},
    {	
    	Print["pmax: Length =",Length[pmaxx],pmaxx,phiA," ",th0," ",phiDet," ",pmaxApertfCorr[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA,th1max]],
    	Length[pmaxx]>1(*Length[pmaxx]!=1 && Length[pmaxx]!=0*) 
    }
    },
    Print["pmax: other case"]
   ]

  ]
  
  
pLimitsApertXListfCorr[yD_?NumericQ, xD_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_,th1max_] := 
(*pLimitsApertXList[yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA] =*)
Module[
  {
   pminXmPi = pminApertCasesfCorr[-Pi, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA,th1max],
   pminX0 = pminApertCasesfCorr[0., yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA,th1max],
   pmaxXmPi = pmaxApertCasesfCorr[-Pi, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA,th1max],
   pmaxX0 = pmaxApertCasesfCorr[0., yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA,th1max],
   plimitlist
   },
	plimitlist=Union[{pminXmPi, pminX0, pmaxXmPi, pmaxX0}];
	Which[
		Length[plimitlist]==1,Join[plimitlist,plimitlist],
		True,plimitlist
	]
]
  


(*p limits from apertY*)

pminApertY[phiA_, yD_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_] =
 Module[{plocal}, Solve[yOff + yAA/2 == yA[phiA, yD, plocal, th0, BRxB, rRxB, rA, rD, phiDet], plocal, Reals][[1, 1, 2, 1]]]
 
pmaxApertY[phiA_, yD_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_] =
Module[{plocal}, Solve[yOff - yAA/2 == yA[phiA, yD, plocal, th0, BRxB, rRxB, rA, rD, phiDet], plocal][[1, 1, 2]]]

pminApertYCases[yD_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_] := Module[
  {phiAlocal = If[yD > yAA/2 + yOff, -Pi/2, Pi/2], pminlocal},
  pminlocal = pminApertY[phiAlocal, yD, th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA];
  Piecewise[{{pminlocal, 0. < pminlocal < pmax}}, 0.]
  ]

pmaxApertYCases[yD_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_] := Module[
  {phiAlocal = If[yD < -yAA/2 + yOff, Pi/2, -Pi/2], pmaxlocal},
  pmaxlocal = pmaxApertY[phiAlocal, yD, th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA];
  Piecewise[{{pmaxlocal, 0. < pmaxlocal < pmax}}, 0.]
  ]

pLimitsApertAllListfCorr[yD_?NumericQ, xD_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_,fCorrf_, {xAA_, xOff_,yAA_, yOff_}] := 
Module[
  {
   pminXmPi = pminApertCasesfCorr[-Pi, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA,fCorrf],
   pminX0 = pminApertCasesfCorr[0., yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA,fCorrf],
   pmaxXmPi = pmaxApertCasesfCorr[-Pi, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA,fCorrf],
   pmaxX0 = pmaxApertCasesfCorr[0., yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA,fCorrf],
   pminY = pminApertYCases[yD, th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA],
   pmaxY = pmaxApertYCases[yD, th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA]
   },
 	Sort[{pminXmPi, pminX0, pmaxXmPi, pmaxX0,pminY,pmaxY}]
  ]



(* ::Section:: *)
(* 2 step-Integration *)

pXDomainfCorr[yD_?NumericQ, xD_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_,th1max_]:=
{p,Sequence@@pLimitsApertXListfCorr[yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA,th1max]}


Step1IntfCorr[b_, yD_, xD_, {phiDet_?NumericQ, th0_?NumericQ}, 
	{alpha_, BRxB_, rRxB_, rA_, rD_, R_, G1_, G2_,th1max_}, {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},
	method1_,PrecGoal1_,AccGoal1_,MinRec1_,MaxRec1_]:=
	NIntegrate[
		ManualphiAIntegrandAllLimitsfCorr[b, yD, xD, {phiDV, phiDet, p, th0}, 
			{alpha, BRxB, rRxB, rA, rD, R, G1, G2,fCorr[p,th0,th1max,rRxB]}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA, yAA, xOff, yOff}],
		Evaluate[pXDomainfCorr[yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA,th1max]],
		{phiDV,-Pi,Pi},
		Method->method1,PrecisionGoal->PrecGoal1,AccuracyGoal->AccGoal1, MinRecursion->MinRec1, MaxRecursion->MaxRec1	
]


Step2IntfCorr[b_, xD_?NumericQ,yD_?NumericQ, {alpha_, BRxB_, rF_, rRxB_, rA_, rD_, R_, G1_, G2_,th1max_}, {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_},
	{xAA_, yAA_, xOff_, yOff_}, method1_,PrecGoal1_,AccGoal1_,MinRec1_,MaxRec1_, method2_,PrecGoal2_,AccGoal2_,MinRec2_,MaxRec2_]:=

		NIntegrate[
		Step1IntfCorr[b, yD, xD, {phiDet, th0}, 
			{alpha, BRxB, rRxB, rA, rD, R, G1, G2,th1max}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA, yAA, xOff, yOff},
			method1,PrecGoal1,AccGoal1,MinRec1,MaxRec1],
			{th0,0.,thetamax[rF]},{phiDet,-Pi,Pi},Method->method2,PrecisionGoal->PrecGoal2,AccuracyGoal->AccGoal2,MinRecursion->MinRec2, MaxRecursion->MaxRec2		
		]
		
		
(* ::Section:: *)
(* Add xD and yD integration *)

BinIntfCorr[b_, OneBinList_, {alpha_, BRxB_, rF_, rRxB_, rA_, rD_, R_, G1_, G2_,th1max_}, {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_},
	{xAA_, yAA_, xOff_, yOff_}, {method1_,PrecGoal1_,AccGoal1_,MinRec1_,MaxRec1_}, {method2_,PrecGoal2_,AccGoal2_,MinRec2_,MaxRec2_}, {method3_,PrecGoal3_,AccGoal3_}]:=
	NIntegrate[
		
		Step2IntfCorr[b, xD, yD, {alpha, BRxB, rF, rRxB, rA, rD, R, G1, G2,th1max}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA, yAA, xOff, yOff},
			method1,PrecGoal1,AccGoal1,MinRec1,MaxRec1,method2,PrecGoal2,AccGoal2,MinRec2,MaxRec2],
			
		{xD, OneBinList[[1,1]], OneBinList[[1,2]]}, {yD, OneBinList[[2,1]], OneBinList[[2,2]]},
		PrecisionGoal->PrecGoal3,Method->method3,MinRecursion->0,MaxRecursion->1,AccuracyGoal->AccGoal3
	]


     
   