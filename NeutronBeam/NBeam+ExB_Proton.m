(* Wolfram Language package *)

Get["Spectra/GlueckProtonSpectrum.m"];
Get["Bfield+Drift/Drifts.m"];
Get["Common/CommonFunctions.m"];
Get["NeutronBeam/NeutronBeamDefinition.m"];
Get["Aperture/BooleApert.m"];


(*Acceleration formulas*)
pacc[p_,U_]:=Sqrt[2*mp*(p^2/(2*mp)-U)]

thDetAcc[p_,th0_,rD_,U_]:=ArcSin[Sqrt[ rD*p^2*Sin[th0]^2/(2*mp*(p^2/(2*mp)-U))]]

(*ExB drift formulas*)
kExB[p_, scale_] := scale*(20 - 6.5*p/ppmax)/50

(* ::Section:: *)
(*forward transfer definitions*)
xGCAF[xDV_, phiDV_, p_, th0_, BRxB_, rRxB_, rA_] := 
	(xDV - rG[p, th0, BRxB/rRxB]*Cos[phiDV])*Sqrt[1/rA]
yGCAF[yDV_, phiDV_, p_, th0_, BRxB_, rRxB_, rA_] := 
	(yDV - rG[p, th0, BRxB/rRxB]*Sin[phiDV])*Sqrt[1/rA]
  

xRxBGCF[phiDV_,xDV_, yDV_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, R_, G1_, G2_]:=
  xGCAF[xDV,phiDV,p,th0,BRxB,rRxB,rA] + D1stBPolyGrad[p, alpha, th0, BRxB, rRxB, yGCAF[yDV, phiDV, p, th0, BRxB, rRxB, rA], R, G1, G2]
  
(*ExB drift applied as well as Det Spread*)
xDetGCF[phiDV_,xDV_, yDV_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, R_, G1_, G2_, rD_, ExBscale_]:=
	(
	xRxBGCF[phiDV,xDV, yDV, p, th0, alpha, BRxB, rRxB, rA, R, G1, G2] 
	- kExB[p, ExBscale]*yGCAF[yDV, phiDV, p, th0, BRxB, rRxB, rA]
	)*Sqrt[rRxB/rD]
	
yDetGCF[phiDV_,xDV_, yDV_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, R_, G1_, G2_, rD_, ExBscale_]:=
	(
	yGCAF[yDV, phiDV, p, th0, BRxB, rRxB, rA] 
	+ kExB[p, ExBscale]*xRxBGCF[phiDV,xDV, yDV, p, th0, alpha, BRxB, rRxB, rA, R, G1, G2]
	)*Sqrt[rRxB/rD]
	
xDetF[phiDV_,xDV_, yDV_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, R_, G1_, G2_, rD_, ExBscale_, U_, phiDet_]:=
	xDetGCF[phiDV,xDV, yDV, p, th0, alpha, BRxB, rRxB, rA, R, G1, G2, rD, ExBscale] + rG[ pacc[p,U], thDetAcc[p,th0,rD,U], rD*BRxB/rRxB] * Cos[phiDet]
	
yDetF[phiDV_,xDV_, yDV_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, R_, G1_, G2_, rD_, ExBscale_, U_, phiDet_]:=
	yDetGCF[phiDV,xDV, yDV, p, th0, alpha, BRxB, rRxB, rA, R, G1, G2, rD, ExBscale] + rG[ pacc[p,U], thDetAcc[p,th0,rD,U], rD*BRxB/rRxB] * Sin[phiDet]
	
(* ::Section:: *)	
(* backwards transfer definitions *)
xDetGCB[xDet_, p_, th0_, phiDet_, BRxB_, rRxB_, rD_, U_]:= xDet - rG[ pacc[p,U], thDetAcc[p,th0,rD,U], rD*BRxB/rRxB ] * Cos[phiDet]
yDetGCB[yDet_, p_, th0_, phiDet_, BRxB_, rRxB_, rD_, U_]:= yDet - rG[ pacc[p,U], thDetAcc[p,th0,rD,U], rD*BRxB/rRxB ] * Sin[phiDet]

xExBGCB[xDet_, p_, th0_, phiDet_, BRxB_, rRxB_, rD_, U_]:= xDetGCB[xDet, p, th0, phiDet, BRxB, rRxB, rD, U]/Sqrt[rRxB/rD]
yExBGCB[yDet_, p_, th0_, phiDet_, BRxB_, rRxB_, rD_, U_]:= yDetGCB[yDet, p, th0, phiDet, BRxB, rRxB, rD, U]/Sqrt[rRxB/rD]

(*xRxBGCB[]:= xExBGCB[] + kExB[p, ExBscale] * yAGCB*)

(* ::Section:: *)
(*hand calculated xDV yDV*)

yDVHandSolve[phiDV_, xDet_, yDet_, p_, th0_, BRxB_, rRxB_, rA_, rD_, ExBscale_, U_, phiDet_] := 
(yExBGCB[yDet, p, th0, phiDet, BRxB, rRxB, rD, U] - 
      kExB[p, ExBscale]*
       xExBGCB[xDet, p, th0, phiDet, BRxB, rRxB, rD, U])/(1 + 
      kExB[p, ExBscale]^2)/Sqrt[1/rA] + 
  rG[p, th0, BRxB/rRxB]*Sin[phiDV]
  
xDVHandSolve[phiDV_, xDet_, yDet_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, R_, G1_, G2_, ExBscale_, U_, phiDet_] :=
 (
    xExBGCB[xDet, p, th0, phiDet, BRxB, rRxB, rD, U]
     + kExB[p, ExBscale]/Sqrt[1/rA]*
      (
       yDVHandSolve[phiDV, xDet, yDet, p, th0, BRxB, rRxB, rA, rD, 
         ExBscale, U, phiDet] - rG[p, th0, BRxB/rRxB]*Sin[phiDV]
       )
     - D1stBPolyGrad[p, alpha, th0, BRxB, rRxB, 
      yGCAF[yDVHandSolve[phiDV, xDet, yDet, p, th0, BRxB, rRxB, rA, 
        rD, ExBscale, U, phiDet], phiDV, p, th0, BRxB, rRxB, rA], R, 
      G1, G2]
    )/Sqrt[1/rA] +
  rG[p, th0, BRxB/rRxB]*Cos[phiDV]

 (*xA for Aperture function*)
 
 
xAAccel[xGCVar_, phiA_, p_, th0_, BRxB_, rRxB_, rA_] := 
	xGCVar + 
	rG[p, th0, rA*BRxB/rRxB]*Cos[phiA]
	
yAAccel[yGCVar_, phiA_, p_, th0_, BRxB_, rRxB_, rA_] := 
	yGCVar + 
	rG[p, th0, rA*BRxB/rRxB]*Sin[phiA]
  
 
  
  
(*Compiled Integrand*)

IntegrandProton2DwNBeamAccelCompiled[
	lambda_, a_, b_, xDet_, yDet_, {phiDV_?NumericQ, phiA_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, 
	{alpha_, BRxB_, rRxB_, rA_,rD_, R_, G1_, G2_, ExBscale_,U_}, 
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] :=
 
 pmomNormedGlueck[p, lambda, a, b]*Sin[th0]*
  
  TrapezNBeamCompiledNormed[
   xDVHandSolve[phiDV, xDet, yDet, p, th0, alpha, BRxB, rRxB, rA, rD, R, G1, G2, ExBscale, U, phiDet]-xOff, {twx, plx, k1x, k2x, k3x}]*
  TrapezNBeamCompiledNormed[
   yDVHandSolve[phiDV, xDet, yDet, p, th0, BRxB, rRxB, rA, rD, ExBscale, U, phiDet]-yOff, {twy, ply,k1y, k2y, k3y}]*
  
  ApertBooleCompiled[
   xAAccel[
   	xGCAF[xDVHandSolve[phiDV, xDet, yDet, p, th0, alpha, BRxB, rRxB, rA, rD, R, G1, G2, ExBscale, U, phiDet],
   	phiDV, p, th0, BRxB, rRxB, rA] , phiA, p, th0, BRxB, rRxB, rA
   ], 
   yAAccel[
   	yGCAF[yDVHandSolve[phiDV, xDet, yDet, p, th0, BRxB, rRxB, rA, rD, ExBscale, U, phiDet], 
   	phiDV, p, th0, BRxB, rRxB, rA], phiA, p, th0, BRxB, rRxB, rA
   ], 
   xAA, yAA, xOff, yOff] 
  
  
(* ::Section::*)(*get limits of integration*)
phiAminApertAccel[xGCVar_, p_, th0_, BRxB_, rRxB_, rA_, xAA_, xOff_] = 
  Solve[xOff + xAA/2 == xAAccel[xGCVar, phiAlocal, p, th0, BRxB, rRxB, rA], phiAlocal][[All, 1, 2, 1, 1]];
phiAmaxApertAccel[xGCVar_, p_, th0_, BRxB_, rRxB_, rA_, xAA_, xOff_] = 
  Solve[xOff - xAA/2 == xAAccel[xGCVar, phiAlocal, p, th0, BRxB, rRxB, rA], phiAlocal][[All, 1, 2, 1, 1]];


(*phiACasesAccel[yD_?NumericQ, xD_?NumericQ, plocal_?NumericQ, 
  th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_,
   G1_, G2_, xOff_, xAA_, ExBscale_, U_] := {phiA, 
  Sequence @@ 
   Re[{phiAmaxApertAccel[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, 
       rD, phiDet, R, G1, G2, xOff, xAA, ExBscale, U][[1]], 
     phiAminApertAccel[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD,
        phiDet, R, G1, G2, xOff, xAA, ExBscale, U][[1]], 
     phiAminApertAccel[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD,
        phiDet, R, G1, G2, xOff, xAA, ExBscale, U][[2]], 
     phiAmaxApertAccel[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD,
        phiDet, R, G1, G2, xOff, xAA, ExBscale, U][[2]]}]}*)

(*phiA limits from Y Aperture*)

phiAminApertfromYAccel[yGCVar_, p_, th0_, BRxB_, rRxB_, rA_, yAA_, yOff_] = Solve[
   yOff + yAA/2 == yAAccel[yGCVar, phiAlocal, p, th0, BRxB, rRxB, rA], phiAlocal][[All, 1, 2, 1, 1 ;; -2]]

phiAmaxApertfromYAccel[yGCVar_, p_, th0_, BRxB_, rRxB_, rA_, yAA_, yOff_] = Solve[
   yOff - yAA/2 == yAAccel[yGCVar, phiAlocal, p, th0, BRxB, rRxB, rA], phiAlocal][[All, 1, 2, 1, 1 ;; -2]]


(*new manual phiA integration with all phi limits*)

ManualphiAIntegrandAllLimitsAccel[lambda_, a_, b_, xDet_, yDet_, {phiDV_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, 
	{alpha_, BRxB_, rRxB_, rA_, rD_, R_, G1_, G2_, ExBscale_, U_}, 
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] := 
 	Module[
 		{
 		phimin = Re[phiAminApertAccel[xGCAF[xDVHandSolve[phiDV, xDet, yDet, If[p == 0., 0.00000001, p], th0, alpha, BRxB, rRxB, rA, rD, R, G1, G2, ExBscale, U, phiDet],
   			phiDV, p, th0, BRxB, rRxB, rA], If[p == 0., 0.00000001, p], th0, BRxB, rRxB, rA, xAA, xOff]], 
   		phimax = Re[phiAmaxApertAccel[xGCAF[xDVHandSolve[phiDV, xDet, yDet, If[p == 0., 0.00000001, p], th0, alpha, BRxB, rRxB, rA, rD, R, G1, G2, ExBscale, U, phiDet],
   			phiDV, p, th0, BRxB, rRxB, rA], If[p == 0., 0.00000001, p], th0, BRxB, rRxB, rA, xAA, xOff]], 
   		phiminY = Re[phiAminApertfromYAccel[yGCAF[yDVHandSolve[phiDV, xDet, yDet, If[p == 0., 0.00000001, p], th0, BRxB, rRxB, rA, rD, ExBscale, U, phiDet], 
   			phiDV, p, th0, BRxB, rRxB, rA], If[p == 0., 0.00000001, p], th0, BRxB, rRxB, rA, yAA, yOff]], 
   		phimaxY = Re[phiAmaxApertfromYAccel[yGCAF[yDVHandSolve[phiDV, xDet, yDet, If[p == 0., 0.00000001, p], th0, BRxB, rRxB, rA, rD, ExBscale, U, phiDet], 
   			phiDV, p, th0, BRxB, rRxB, rA], If[p == 0., 0.00000001, p], th0, BRxB, rRxB, rA, yAA, yOff]], 
   		IntValues, philimitlist, philimitlistSorted, Deltalimits, phiCenters, Integrated
   		}, 
  		philimitlist = Flatten[{-Pi // N, phimin, phimax, phiminY, phimaxY, Pi // N}];
		philimitlist = If[# > Pi, # - 2 Pi, #] & /@ philimitlist;
  		philimitlistSorted = Union[philimitlist];
  		Deltalimits = Table[philimitlistSorted[[i + 1]] - philimitlistSorted[[i]], {i, 1, Length[philimitlistSorted] - 1}];
  		phiCenters = Table[philimitlistSorted[[i]] + Deltalimits[[i]]/2, {i, 1, Length[philimitlistSorted] - 1}];
  		IntValues = 
   			IntegrandProton2DwNBeamAccelCompiled[lambda, a, b, xDet, yDet, {phiDV, #, phiDet, p, th0}, {alpha, BRxB, rRxB, rA, rD, R, G1, G2, ExBscale, U}, 
   				{twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA, yAA, xOff, yOff}] & /@ phiCenters;
  		Integrated = Total[Deltalimits*IntValues];
  		Integrated
  	]




(* ::Section:: *)
(* Integrand with p limits from aperture *)


pminApertAccel[(*phiDV_, *)phiA_?NumericQ, xDet_?NumericQ, yDet_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_, ExBscale_, U_] :=
Module[
	{plocal}, 
 	NSolve[xOff + xAA/2 == 
 		xAAccel[
   			xGCAF[xDVHandSolve[phiDV, xDet, yDet, plocal, th0, alpha, BRxB, rRxB, rA, rD, R, G1, G2, ExBscale, U, phiDet],
   				phiDV, plocal, th0, BRxB, rRxB, rA] , phiA, plocal, th0, BRxB, rRxB, rA],
   		plocal][[All, 1, 2]]
]

pmaxApertAccel[(*phiDV_, *)phiA_?NumericQ, xDet_?NumericQ, yDet_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_, ExBscale_, U_] :=
Module[
	{plocal}, 
 	NSolve[xOff - xAA/2 == 
 		xAAccel[
   			xGCAF[xDVHandSolve[phiDV, xDet, yDet, plocal, th0, alpha, BRxB, rRxB, rA, rD, R, G1, G2, ExBscale, U, phiDet],
   				phiDV, plocal, th0, BRxB, rRxB, rA] , phiA, plocal, th0, BRxB, rRxB, rA],
   		plocal][[All, 1, 2]]
]
      
      
pminApertCasesAccel[(*phiDV_, *)phiA_?NumericQ, xDet_?NumericQ, yDet_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_, ExBscale_, U_] := 
Module[
  {pmin = Cases[pminApertAccel[(*phiDV, *)phiA, xDet, yDet, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, ExBscale, U], _?(# \[Element] Reals && 0. <= # <= 2*ppmax &)]},
  Piecewise[
   {
    {0, Length[pmin] == 0},
    {ppmax, Length[pmin] == 1 && pmin[[1]] > ppmax},
    {pmin[[1]], Length[pmin] == 1},
    {
    	Print["pmin: Length =",Length[pmin],pmin,phiA," ",th0," ",phiDet," ",pminApertAccel[(*phiDV, *)phiA, xDet, yDet, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, ExBscale, U]], 
    	Length[pmin]>1(*Length[pmin]!=1 && Length[pmin]!=0*) 
    }
    },
    Print["pmin: other Case"]
   ]
 
]

pmaxApertCasesAccel[(*phiDV_, *)phiA_?NumericQ, xDet_?NumericQ, yDet_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_, ExBscale_, U_] := 
Module[
  {pmaxx =Cases[pmaxApertAccel[(*phiDV, *)phiA, xDet, yDet, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, ExBscale, U], _?(# \[Element] Reals && 0. <= # <= 2*ppmax &)]},
  
  Piecewise[
   {
    {ppmax, Length[pmaxx] == 0},
    {ppmax, Length[pmaxx] == 1 && pmaxx[[1]] > ppmax},
    {pmaxx[[1]], Length[pmaxx] == 1},
    {	
    	Print["pmax: Length =",Length[pmaxx],pmaxx,phiA," ",th0," ",phiDet," ",pmaxApertAccel[(*phiDV, *)phiA, xDet, yDet, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, ExBscale, U]],
    	Length[pmaxx]>1(*Length[pmaxx]!=1 && Length[pmaxx]!=0*) 
    }
    },
    Print["pmax: other case"]
   ]

  ]
  
  
pLimitsApertXListAccel[(*phiDV_, *)xDet_?NumericQ, yDet_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_, ExBscale_, U_] := 
(*pLimitsApertXList[yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA] =*)
Module[
  {
   pminXmPi = pminApertCasesAccel[ -Pi, xDet, yDet, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, ExBscale, U],
   pminX0 = pminApertCasesAccel[ 0., xDet, yDet, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, ExBscale, U],
   pmaxXmPi = pmaxApertCasesAccel[ -Pi, xDet, yDet, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, ExBscale, U],
   pmaxX0 = pmaxApertCasesAccel[ 0., xDet, yDet, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, ExBscale, U],
   plimitlist
   },
	plimitlist=Union[{pminXmPi, pminX0, pmaxXmPi, pmaxX0}];
	Which[
		Length[plimitlist]==1,Join[plimitlist,plimitlist],
		True,plimitlist
	]
]


(* ::Section:: *)
(* 2 step-Integration *)

pXDomainAccel[xDet_?NumericQ, yDet_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_, ExBscale_, U_]:=
{p,Sequence@@pLimitsApertXListAccel[xDet, yDet, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, ExBscale, U]}


Step1IntAccel[lambda_, a_, b_, xDet_, yDet_, {phiDet_?NumericQ, th0_?NumericQ}, 
	{alpha_, BRxB_, rRxB_, rA_, rD_, R_, G1_, G2_, ExBscale_, U_}, {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},
	method1_,PrecGoal1_,AccGoal1_,MinRec1_,MaxRec1_]:=
	NIntegrate[
		ManualphiAIntegrandAllLimitsAccel[lambda, a, b, xDet, yDet, {phiDV, phiDet, p, th0}, 
			{alpha, BRxB, rRxB, rA, rD, R, G1, G2, ExBscale, U}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA, yAA, xOff, yOff}],
		Evaluate[pXDomainAccel[xDet, yDet, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, ExBscale, U]],
		{phiDV,-Pi,Pi},
		Method->method1,PrecisionGoal->PrecGoal1,AccuracyGoal->AccGoal1, MinRecursion->MinRec1, MaxRecursion->MaxRec1	
]


Step2IntAccel[lambda_, a_, b_, xDet_?NumericQ,yDet_?NumericQ, 
	{alpha_, BRxB_, rF_, rRxB_, rA_, rD_, R_, G1_, G2_, ExBscale_, U_}, {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_},
	{xAA_, yAA_, xOff_, yOff_}, method1_,PrecGoal1_,AccGoal1_,MinRec1_,MaxRec1_, method2_,PrecGoal2_,AccGoal2_,MinRec2_,MaxRec2_]:=

		NIntegrate[
		Step1IntAccel[lambda, a, b, xDet, yDet, {phiDet, th0}, 
			{alpha, BRxB, rRxB, rA, rD, R, G1, G2, ExBscale, U}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA, yAA, xOff, yOff},
			method1,PrecGoal1,AccGoal1,MinRec1,MaxRec1],
			{th0,0.,thetamax[rF]},{phiDet,-Pi,Pi},Method->method2,PrecisionGoal->PrecGoal2,AccuracyGoal->AccGoal2,MinRecursion->MinRec2, MaxRecursion->MaxRec2		
		]
		
		
(* ::Section:: *)
(* Add xD and yD integration *)

BinIntAccel[lambda_, a_, b_, OneBinList_, {alpha_, BRxB_, rF_, rRxB_, rA_, rD_, R_, G1_, G2_, ExBscale_, U_}, {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_},
	{xAA_, yAA_, xOff_, yOff_}, {method1_,PrecGoal1_,AccGoal1_,MinRec1_,MaxRec1_}, {method2_,PrecGoal2_,AccGoal2_,MinRec2_,MaxRec2_}, {method3_,PrecGoal3_,AccGoal3_}]:=
	NIntegrate[
		
		Step2IntAccel[lambda, a, b, xDet, yDet, {alpha, BRxB, rF, rRxB, rA, rD, R, G1, G2, ExBscale, U}, 
			{twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA, yAA, xOff, yOff},
			method1,PrecGoal1,AccGoal1,MinRec1,MaxRec1,method2,PrecGoal2,AccGoal2,MinRec2,MaxRec2],
			
		{xDet, OneBinList[[1,1]], OneBinList[[1,2]]}, {yDet, OneBinList[[2,1]], OneBinList[[2,2]]},
		PrecisionGoal->PrecGoal3,Method->method3,MinRecursion->0,MaxRecursion->1,AccuracyGoal->AccGoal3
	]
  
  