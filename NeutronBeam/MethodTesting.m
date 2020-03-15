(* Wolfram Language package *)

Get["NeutronBeam/NBeamIntegrand.m"];



(* ::Section:: *)
(*2D*)
(*Manual phiA integration V2 with all phiA limits, only phiDV and p Int*)
(*no Table, specific yD and xD*)
(*Added Integration Options*)

Int2DwNBeamCompiledManualphiAAllLimitsPhiDVAndMomentumInt[b_, yD_, xD_,{phiDet_?NumericQ,th0_?NumericQ}, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_,Intopts_:OptionsPattern[NIntegrate]]:=

  Reap[NIntegrate[
   ManualphiAIntegrandAllLimits[b, yD, xD, {phiDV, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}],
     
     (*{th0, 0., thetamax[rF]}, {phiDet, -Pi//N, Pi//N},*) {p, 0., pmax}, {phiDV, -Pi//N, Pi//N}, 
   PrecisionGoal -> IntPrec, AccuracyGoal -> 5, Sequence@Intopts,
   EvaluationMonitor:>Sow[{p,phiDV}]
  ]
]
   
   
(*now with the p boundaries*)
Int2DwNBeamCompiledManualphiAAllLimitsPhiDVAndMomentumIntpLimits[b_, yD_, xD_,{phiDet_?NumericQ,th0_?NumericQ}, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_,Intopts_:OptionsPattern[NIntegrate]]:=

  Reap[NIntegrate[
   ManualphiAIntegrandAllLimits[b, yD, xD, {phiDV, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}],
     
     (*{th0, 0., thetamax[rF]}, {phiDet, -Pi//N, Pi//N},*) 
     {
     	p, 
     	0.,
     	Sequence@@pLimitsApertXList[yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA], 
     	pmax
     },
     {phiDV, -Pi//N, Pi//N}, 
   PrecisionGoal -> IntPrec, AccuracyGoal -> 5, Sequence@Intopts,
   EvaluationMonitor:>Sow[{p,phiDV}]
  ]
   
]

(* ::Section:: *)
(* 3D *)
(*with phiDet integration*)
Int2DwNBeamCompiledManualphiAAllLimitsPhiDVMomentumPhiDetInt[b_, yD_, xD_,{th0_?NumericQ}, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_,Intopts_:OptionsPattern[NIntegrate]]:=

 NIntegrate[
   ManualphiAIntegrandAllLimits[b, yD, xD, {phiDV, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}],
     
     (*{th0, 0., thetamax[rF]},*) {phiDet, -Pi//N, Pi//N}, {p, 0., pmax}, {phiDV, -Pi//N, Pi//N}, 
   PrecisionGoal -> IntPrec, AccuracyGoal -> 5, Sequence@Intopts
  ]
   
   
(*now with the p boundaries from ApertX*)
Int2DwNBeamCompiledManualphiAAllLimitsPhiDVMomentumPhiDetIntpLimits[b_, yD_, xD_,{th0_?NumericQ}, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_,Intopts_:OptionsPattern[NIntegrate]]:=

  NIntegrate[
   ManualphiAIntegrandAllLimits[b, yD, xD, {phiDV, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}],
     
     (*{th0, 0., thetamax[rF]},*) {phiDet, -Pi//N, Pi//N}, 
     {
     	p, 
     	0.,
     	Sequence@@pLimitsApertXList[yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA], 
     	pmax
     },
     {phiDV, -Pi//N, Pi//N}, 
   PrecisionGoal -> IntPrec, AccuracyGoal -> 5, Sequence@Intopts
]
  

  
(*now with the ALL p boundaries from Apert*)
Int2DwNBeamCompiledManualphiAAllLimitsPhiDVMomentumPhiDetIntpLimitsAll[b_, yD_, xD_,{th0_?NumericQ}, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_,Intopts_:OptionsPattern[NIntegrate]]:=

  NIntegrate[
   ManualphiAIntegrandAllLimits[b, yD, xD, {phiDV, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}],
     
     (*{th0, 0., thetamax[rF]},*) {phiDet, -Pi//N, Pi//N}, 
     {
     	p, 
     	0.,
     	Sequence@@pLimitsApertAllList[yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, {xAA, xOff, yAA, yOff}], 
     	pmax
     },
     {phiDV, -Pi//N, Pi//N}, 
   PrecisionGoal -> IntPrec, AccuracyGoal -> 5, Sequence@Intopts
]


(* ::Section:: *)
(* 4D integrals *)
(*now with the p boundaries from ApertX*)
Int2DwNBeamCompiledManualphiAAllIntpXLimits[b_, yD_, xD_, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_,Intopts_:OptionsPattern[NIntegrate]]:=

  NIntegrate[
   ManualphiAIntegrandAllLimits[b, yD, xD, {phiDV, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}],
     
     {th0, 0., thetamax[rF]}, {phiDet, -Pi//N, Pi//N}, 
     {
     	p,
     	0., 
      	Sequence@@pLimitsApertXList[yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA],
      	pmax
     },
     {phiDV, -Pi//N, Pi//N}, 
   PrecisionGoal -> IntPrec, AccuracyGoal -> 5, Sequence@Intopts
]



(*now without the p boundaries from ApertX*)
Int2DwNBeamCompiledManualphiAAllInt[b_, yD_, xD_, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_,Intopts_:OptionsPattern[NIntegrate]]:=

  NIntegrate[
   ManualphiAIntegrandAllLimits[b, yD, xD, {phiDV, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}],
     
     {th0, 0., thetamax[rF]}, {phiDet, -Pi//N, Pi//N}, 
     {
     	p,
     	0., 
      	pmax
     },
     {phiDV, -Pi//N, Pi//N}, 
   PrecisionGoal -> IntPrec, AccuracyGoal -> 5, Sequence@Intopts
]
