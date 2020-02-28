(* Wolfram Language package *)

Get["NeutronBeam/NBeamIntegrand.m"];




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
	
	
	
(*with phiA limits*)

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



(*Manual phiA integration*)
Int2DwNBeamCompiledManualphiA[b_, XYData_List, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_]:=

ParallelTable[
  NIntegrate[
   ManualphiAIntegrand[b, XYData[[bin,2]], XYData[[bin,1]], {phiDV, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}],
     
     {th0, 0, thetamax[rF]}, {phiDet, -Pi, Pi}, {p, 0, pmax}, {phiDV, -Pi, Pi}, 
   PrecisionGoal -> IntPrec, Method -> "AdaptiveMonteCarlo", AccuracyGoal -> 5, MaxPoints -> 10000000], 
   	
{bin,1,Length[XYData]}, Method -> "FinestGrained"]


(*Manual phiA integration V2 with all phiA limits*)
Int2DwNBeamCompiledManualphiAAllLimits[b_, XYData_List, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_]:=

ParallelTable[
  NIntegrate[
   ManualphiAIntegrandAllLimits[b, XYData[[bin,2]], XYData[[bin,1]], {phiDV, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}],
     
     {th0, 0., thetamax[rF]}, {phiDet, -Pi//N, Pi//N}, {p, 0., pmax}, {phiDV, -Pi//N, Pi//N}, 
   PrecisionGoal -> IntPrec, Method -> "LocalAdaptive", AccuracyGoal -> 5,MinRecursion->3,MaxRecursion->10], 
   	
{bin,1,Length[XYData]}, Method -> "FinestGrained"]


(*Manual phiA integration V2 with all phiA limits + p limits from APerture X*)
Int2DwNBeamCompiledManualphiAAllLimitspLimitsX[b_, XYData_List, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_]:=

ParallelTable[
	Module[
		{xbin=XYData[[bin,1]],ybin=XYData[[bin,2]]},
  NIntegrate[
   ManualphiAIntegrandAllLimits[b, ybin, xbin, {phiDV, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}],
     
     {th0, 0., thetamax[rF]}, {phiDet, -Pi//N, Pi//N},
     {
     	p, 
     	0.,
     	Sequence@@pLimitsApertXList[ybin, xbin, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA],
     	pmax
     }, 
     {phiDV, -Pi//N, Pi//N}, 
   PrecisionGoal -> IntPrec, Method -> "LocalAdaptive", AccuracyGoal -> 5,MinRecursion->3,MaxRecursion->10]
	], 
   	
{bin,1,Length[XYData]}, Method -> "CoarsestGrained"]


(*Manual phiA integration V2 with all phiA limits + p limits from APerture X*)
(*approach 2 where the p intermediate limits are given explicitly*)
Int2DwNBeamCompiledManualphiAAllLimitspLimitsXManualInt[b_, XYData_List,{phiDet_?NumericQ,th0_?NumericQ}, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_]:=

ParallelTable[
	Module[
		{xbin=XYData[[bin,1]],ybin=XYData[[bin,2]]},
  NIntegrate[
   ManualphiAIntegrandAllLimits[b, ybin, xbin, {phiDV, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff},False],
     
     (*{th0, 0., thetamax[rF]}, {phiDet, -Pi//N, Pi//N},*)
     {
     	p, 
     	0.,
     	Sequence@@pLimitsApertXList[ybin, xbin, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA](*[[1]]*),
     	(*pLimitsApertXList[ybin, xbin, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA][[2]],
     	pLimitsApertXList[ybin, xbin, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA][[3]],
     	pLimitsApertXList[ybin, xbin, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA][[4]],*)
     	pmax
     }, 
     {phiDV, -Pi//N, Pi//N}, 
   PrecisionGoal -> IntPrec, Method -> "LocalAdaptive", AccuracyGoal -> 5,MinRecursion->3,MaxRecursion->10]
	], 
   	
{bin,1,Length[XYData]}, Method -> "CoarsestGrained"]





(*Integration with p limits*)

  
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
  
	
	
	