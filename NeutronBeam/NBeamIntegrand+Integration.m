(* Wolfram Language package *)
Get["Spectra/ElectronSpectrum.m"];
Get["Bfield+Drift/Bfields.m"]
Get["Bfield+Drift/Drifts.m"];
Get["Common/CommonFunctions.m"];
Get["NeutronBeam/NeutronBeamDefinition.m"];


ApertBoole[xtraj_, ytraj_, {xA_, yA_, xOff_, yOff_}] := 
 Boole[xtraj < xOff + xA/2 && xtraj > xOff - xA/2 && 
   ytraj < yOff + yA/2 && ytraj > yOff - yA/2]
   
ApertBooleCompiled = Compile[
  {{xtraj, _Real}, {ytraj, _Real}, {xA, _Real}, {yA, _Real}, {xOff, _Real}, {yOff, _Real}},
  Boole[xtraj < xOff + xA/2 && xtraj > xOff - xA/2 && ytraj < yOff + yA/2 && ytraj > yOff - yA/2], 
  CompilationTarget -> "C", RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"
  ];   
   
   
   
   
   
xGCA[xDV_, phiDV_, p_, th0_, BRxB_, rRxB_, 
  rA_] := (xDV - rG[p, th0, BRxB/rRxB]*Cos[phiDV])*Sqrt[1/rA]
yGCA[yDV_, phiDV_, p_, th0_, BRxB_, rRxB_, 
  rA_] := (yDV - rG[p, th0, BRxB/rRxB]*Sin[phiDV])*Sqrt[1/rA]
  
  
DeltayDV[phiDV_, yD_, p_, th0_, BRxB_, rRxB_, rA_, rD_, 
  phiDet_] := (yD - rG[p, th0, rD*BRxB/rRxB]*Sin[phiDet])*
   Sqrt[rD/rRxB]*Sqrt[rA] + rG[p, th0, BRxB/rRxB]*Sin[phiDV]
DeltaxDV[phiDV_, yDV_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, 
  phiDet_, R_, G1_, 
  G2_] := ((xD - rG[p, th0, rD*BRxB/rRxB]*Cos[phiDet])*Sqrt[rD/rRxB] -
      D1stBPolyGrad[p, alpha, th0, BRxB, rRxB, 
      yGCA[yDV, phiDV, p, th0, BRxB, rRxB, rA], R, G1, G2])*Sqrt[rA] +
   rG[p, th0, BRxB/rRxB]*Cos[phiDV]
   
   
DeltaxDV2[phiDV_, yD_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, 
  phiDet_, R_, G1_, G2_] := 
 DeltaxDV[phiDV, 
  DeltayDV[phiDV, yD, p, th0, BRxB, rRxB, rA, rD, phiDet], xD, p, th0,
   alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2]
   
   
xA[phiA_, yD_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_,
   R_, G1_, G2_] := 
 xGCA[DeltaxDV2[phiDV, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, 
    phiDet, R, G1, G2], phiDV, p, th0, BRxB, rRxB, rA] + 
  rG[p, th0, rA*BRxB/rRxB]*Cos[phiA]
  
yA[phiA_, yD_, p_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_] := 
 yGCA[DeltayDV[phiDV, yD, p, th0, BRxB, rRxB, rA, rD, phiDet], phiDV, 
   p, th0, BRxB, rRxB, rA] + rG[p, th0, rA*BRxB/rRxB]*Sin[phiA]
   
   
   
   
Integrand2DwNBeamv31[b_, yD_, 
  xD_, {phiDV_?NumericQ, phiA_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, {alpha_, BRxB_, rRxB_, rA_,
    rD_, R_, G1_, G2_}, {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, 
   k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] :=
 
 wmomNormedWb[\[Lambda]0, \[Kappa]0, b, p]*Sin[th0]*
  
  TrapezNBeamNewNormed[
   DeltaxDV2[phiDV, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,
     R, G1, G2]-xOff, {twx, plx, k1x, k2x, k3x}]*
  TrapezNBeamNewNormed[
   DeltayDV[phiDV, yD, p, th0, BRxB, rRxB, rA, rD, phiDet]-yOff, {twy, ply,
     k1y, k2y, k3y}]*
  
  ApertBoole[
   xA[phiA, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, 
    G2], yA[phiA, yD, p, th0, BRxB, rRxB, rA, rD, phiDet], {xAA, yAA, 
    xOff, yOff}]   
   
   
Integration2DwNBeamc31[b_, XYData_List, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	 {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_]:=

ParallelTable[
  NIntegrate[
   Integrand2DwNBeamv31[b, XYData[[bin,2]], XYData[[bin,1]],
   	 {phiDV, phiA, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}
   	 ],
     
     {th0, 0, thetamax[rF]}, {phiDet, 0, 2*Pi}, {phiA, 0, 2 Pi}, {phiDV, 0, 2*Pi},{p, 0, pmax},
   PrecisionGoal -> IntPrec, Method -> "AdaptiveMonteCarlo", AccuracyGoal -> 5, MaxPoints -> 2500000]
	, {bin,1,Length[XYData]}, Method -> "FinestGrained"]
  
  
  
Integration2DwNBeamc31WOTheta[b_, XYData_List, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	 {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_]:=

ParallelTable[
  NIntegrate[
   Integrand2DwNBeamv31[b, XYData[[bin,2]], XYData[[bin,1]],
   	 {phiDV, phiA, phiDet, p, Pi/4},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}
   	 ],
     
     {phiDet, 0, 2*Pi}, {phiA, 0, 2 Pi}, {phiDV, 0, 2*Pi},{p, 0, pmax},
   PrecisionGoal -> IntPrec, Method -> "AdaptiveMonteCarlo", AccuracyGoal -> 5, MaxPoints -> 2500000]
	, {bin,1,Length[XYData]}, Method -> "FinestGrained"]


(*Compiled Integrand*)

Integrand2DwNBeamCompiled[
	b_, yD_, xD_, {phiDV_?NumericQ, phiA_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, {alpha_, BRxB_, rRxB_, rA_,rD_, R_, G1_, G2_}, 
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] :=
 
 wmomNormedWb[\[Lambda]0, \[Kappa]0, b, p]*Sin[th0]*
  
  TrapezNBeamCompiledNormed[
   DeltaxDV2[phiDV, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2]-xOff, {twx, plx, k1x, k2x, k3x}]*
  TrapezNBeamCompiledNormed[
   DeltayDV[phiDV, yD, p, th0, BRxB, rRxB, rA, rD, phiDet]-yOff, {twy, ply,k1y, k2y, k3y}]*
  
  ApertBooleCompiled[
   xA[phiA, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2], yA[phiA, yD, p, th0, BRxB, rRxB, rA, rD, phiDet], 
   	xAA, yAA, xOff, yOff]  

Integration2DwNBeamCompiled[b_, XYData_List, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	 {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_]:=

ParallelTable[
  NIntegrate[
   Integrand2DwNBeamCompiled[b, XYData[[bin,2]], XYData[[bin,1]], {phiDV, phiA, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}
   	 ],
     
     {th0, 0, thetamax[rF]}, {phiDet, 0, 2*Pi}, {phiA, 0, 2 Pi}, {phiDV, 0, 2*Pi},{p, 0, pmax},
   PrecisionGoal -> IntPrec, Method -> "AdaptiveMonteCarlo", AccuracyGoal -> 5, MaxPoints -> 10000000]
	, {bin,1,Length[XYData]}, Method -> "FinestGrained"]



(* ::Section:: *)
(* Integrand with phiA limits and compiled integrand *)




phiAminApert[yD_, xD_, plocal_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_] = 
  Solve[xOff + xAA/2 == xA[phiAlocal, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2], phiAlocal][[All, 1, 2, 1, 1]];
phiAmaxApert[yD_, xD_, plocal_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_] = 
   Solve[xOff - xAA/2 == xA[phiAlocal, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2], phiAlocal][[All, 1, 2, 1, 1]];


phiACases[yD_?NumericQ, xD_?NumericQ, plocal_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_] :=
 {phiA,
 	Sequence@@Re[{
  phiAmaxApert[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[1]],
  phiAminApert[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[1]],
  phiAminApert[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[2]],
  phiAmaxApert[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[2]]
  }]
 }

(*we take the same integrand as above*)
Int2DwNBeamCompiledphiALimits[b_, XYData_List, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_]:=

ParallelTable[
  NIntegrate[
   Integrand2DwNBeamCompiled[b, XYData[[bin,2]], XYData[[bin,1]], {phiDV, phiA, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}
   	 ],
     
     {th0, 0, thetamax[rF]}, {phiDet, 0, 2*Pi}, {p, 0, pmax}, {phiDV, 0, 2*Pi}, 
     {
     	phiA,
     	Re[phiAmaxApert[XYData[[bin,2]], XYData[[bin,1]], p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[1]]],
     	Re[phiAminApert[XYData[[bin,2]], XYData[[bin,1]], p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[1]]],
     	Re[phiAminApert[XYData[[bin,2]], XYData[[bin,1]], p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[2]]],
     	Re[phiAmaxApert[XYData[[bin,2]], XYData[[bin,1]], p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[2]]]
     },
   PrecisionGoal -> IntPrec, Method -> "AdaptiveMonteCarlo", AccuracyGoal -> 5, MaxPoints -> 10000000], 
   	
{bin,1,Length[XYData]}, Method -> "FinestGrained"]


(*GLOBAL ADAPTIVE METHOD*)
Int2DwNBeamCompiledphiALimitsGLOBADAP[b_, XYData_List, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_]:=

ParallelTable[
  NIntegrate[
   Integrand2DwNBeamCompiled[b, XYData[[bin,2]], XYData[[bin,1]], {phiDV, phiA, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}
   	 ],
     
     {th0, 0, thetamax[rF]}, {phiDet, 0, 2*Pi}, {p, 0, pmax}, {phiDV, 0, 2*Pi}, 
     {
     	phiA,
     	Re[phiAmaxApert[XYData[[bin,2]], XYData[[bin,1]], p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[1]]],
     	Re[phiAminApert[XYData[[bin,2]], XYData[[bin,1]], p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[1]]],
     	Re[phiAminApert[XYData[[bin,2]], XYData[[bin,1]], p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[2]]],
     	Re[phiAmaxApert[XYData[[bin,2]], XYData[[bin,1]], p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[2]]]
     },
   PrecisionGoal -> IntPrec, Method -> {"GlobalAdaptive", "MaxErrorIncreases" -> 5000}, MinRecursion-> 5, MaxRecursion->15, AccuracyGoal -> 5], 
   	
{bin,1,Length[XYData]}, Method -> "FinestGrained"]




(*same thing with separated integrals for the separated phiA regions*)
Int2DwNBeamCompiledphiALimitsSplit[b_, XYData_List, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_,MaxPointN_Integer]:=

ParallelTable[
  NIntegrate[
   Integrand2DwNBeamCompiled[b, XYData[[bin,2]], XYData[[bin,1]], {phiDV, phiA, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}
   	 ],
     
     {th0, 0., thetamax[rF]}, {phiDet, -Pi, Pi}, {p, 0., pmax}, {phiDV, -Pi, Pi}, 
     {
     	phiA,
     	Re[phiAmaxApert[XYData[[bin,2]], XYData[[bin,1]], p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[1]]],
     	Re[phiAminApert[XYData[[bin,2]], XYData[[bin,1]], p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[1]]]
     },
   PrecisionGoal -> IntPrec, Method -> "AdaptiveMonteCarlo", AccuracyGoal -> 5, MaxPoints -> MaxPointN/2]
   +
   NIntegrate[
   Integrand2DwNBeamCompiled[b, XYData[[bin,2]], XYData[[bin,1]], {phiDV, phiA, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}
   	 ],
     
     {th0, 0., thetamax[rF]}, {phiDet, -Pi, Pi}, {p, 0., pmax}, {phiDV, -Pi, Pi}, 
     {
     	phiA,
     	Re[phiAmaxApert[XYData[[bin,2]], XYData[[bin,1]], p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[1]]],
     	Re[phiAminApert[XYData[[bin,2]], XYData[[bin,1]], p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[1]]],
     	Re[phiAminApert[XYData[[bin,2]], XYData[[bin,1]], p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[2]]],
     	Re[phiAmaxApert[XYData[[bin,2]], XYData[[bin,1]], p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[2]]]
     },
   PrecisionGoal -> IntPrec, Method -> "AdaptiveMonteCarlo", AccuracyGoal -> 5, MaxPoints -> MaxPointN/2], 
   	
{bin,1,Length[XYData]}, Method -> "FinestGrained"]



(*Integrand with manual summation of phiA integration*)
ManualphiAIntegrand[b_, yD_, xD_, {phiDV_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, 
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
  ]

Int2DwNBeamCompiledManualphiA[b_, XYData_List, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_]:=

ParallelTable[
  NIntegrate[
   ManualphiAIntegrand[b, XYData[[bin,2]], XYData[[bin,1]], {phiDV, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}],
     
     {th0, 0, thetamax[rF]}, {phiDet, -Pi, Pi}, {p, 0, pmax}, {phiDV, -Pi, Pi}, 
   PrecisionGoal -> IntPrec, Method -> "AdaptiveMonteCarlo", AccuracyGoal -> 5, MaxPoints -> 10000000], 
   	
{bin,1,Length[XYData]}, Method -> "FinestGrained"]




(* ::Section:: *)
(* Integrand with p limits from aperture *)


pminApert[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_] :=
Module[
	{plocal}, 
 	NSolve[xOff + xAA/2 == xA[phiA, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1,G2], plocal][[All, 1, 2]]
]
      
      
pminApertCases[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_] := 
Module[
  {pmin = Cases[pminApert[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA], _?(# \[Element] Reals && 0. <= # <= 2*pmax &)]},
  Piecewise[
   {
    {0, Length[pmin] == 0},
    (*{pmax, Length[pmin] == 1 && pmin[[1]] > pmax},*)
    {pmin[[1]], Length[pmin] == 1},
    {
    	Print["pmin: Length =",Length[pmin],pmin,phiA," ",th0," ",phiDet," ",pminApert[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA]], 
    	Length[pmin]>1(*Length[pmin]!=1 && Length[pmin]!=0*) 
    }
    },
    Print["mup"]
   ]
 
]
  
  
  
pmaxApert[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_] :=
	Module[
		{plocal}, 
 		NSolve[xOff - xAA/2 == xA[phiA, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1,G2], plocal][[All, 1, 2]]
	]


pmaxApertCases[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_] := 
Module[
  {pmaxx =Cases[pmaxApert[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA], _?(# \[Element] Reals && 0. <= # <= 2*pmax &)]},
  
  Piecewise[
   {
    {pmax, Length[pmaxx] == 0},
    {pmax, Length[pmaxx] == 1 && pmaxx[[1]] > pmax},
    {pmaxx[[1]], Length[pmaxx] == 1},
    {	
    	Print["pmax: Length =",Length[pmaxx],pmaxx,phiA," ",th0," ",phiDet," ",pmaxApert[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA]],
    	Length[pmaxx]>1(*Length[pmaxx]!=1 && Length[pmaxx]!=0*) 
    }
    },
    Print["mep"]
   ]

  ]


  
Integration2DwNBeamc31wpLimits[b_, XYData_List, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	 {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_]:=

ParallelTable[
  NIntegrate[
   Integrand2DwNBeamv31[b, XYData[[bin,2]], XYData[[bin,1]],
   	 {phiDV, phiA, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}
   	 ],
     
     {th0, 0, thetamax[rF]}, {phiDet, 0, 2*Pi}, {phiA, 0, 2 Pi}, {phiDV, 0, 2*Pi},
     {
     	p,
     	pminApertCases[phiA, XYData[[bin,2]], XYData[[bin,1]], th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA], 
     	pmaxApertCases[phiA, XYData[[bin,2]], XYData[[bin,1]], th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA]
     },
   PrecisionGoal -> IntPrec, Method -> "AdaptiveMonteCarlo", AccuracyGoal -> 5, MaxPoints -> 2500000]
	, {bin,1,Length[XYData]}, Method -> "FinestGrained"]


Integration2DwNBeamc31wpLimitsWOTheta[b_, XYData_List, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	 {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_]:=

ParallelTable[
  NIntegrate[
   Integrand2DwNBeamv31[b, XYData[[bin,2]], XYData[[bin,1]],
   	 {phiDV, phiA, phiDet, p, Pi/4},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}
   	 ],
     
     {phiDet, 0, 2*Pi}, {phiA, 0, 2 Pi}, {phiDV, 0, 2*Pi},
     {
     	p,
     	pminApertCases[phiA, XYData[[bin,2]], XYData[[bin,1]], Pi/4, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA], 
     	pmaxApertCases[phiA, XYData[[bin,2]], XYData[[bin,1]], Pi/4, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA]
     },
   PrecisionGoal -> IntPrec, Method -> "AdaptiveMonteCarlo", AccuracyGoal -> 5, MaxPoints -> 2500000]
	, {bin,1,Length[XYData]}, Method -> "FinestGrained"]
  
     
   