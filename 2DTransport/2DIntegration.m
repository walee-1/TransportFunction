(* Wolfram Language package *)

Get["2DTransport/2DIntegrand.m"];


Integration2D[b_?NumericQ, {alpha_, BRxB_, rRxB_, rD_,rF_}, {xA_, yA_, xOff_, yOff_, rA_,prec_},BinN_,XYList_,IntPrec_]:=
ParallelTable[
	NIntegrate[
	   Integrand2D[b, {XYList[[bin,1]],XYList[[bin,2]],phi, p, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff, yOff, rA,prec}],
	   {th0,0,thetamax[rF]},
	   {phi,0,2*Pi},
	   {p,
	   	pmin2DCases[{XYList[[bin,1]], phi, th0, BRxB, alpha, rA, rRxB,rD}, {xA,xOff}],
	   	pmax2DCases[{XYList[[bin,1]], phi, th0, BRxB, alpha, rA, rRxB,rD}, {xA,xOff}]},
	   PrecisionGoal->IntPrec
	],
	{bin,1,BinN},Method->"FinestGrained"
]

Integration2DApertMC[b_, {alpha_, BRxB_, rRxB_, rD_,rF_}, {xA_, yA_, xOff_, yOff_, rA_,steps_},BinN_,XYList_,IntPrec_]:=
ParallelTable[
	NIntegrate[
	   Integrand2DApertMC[b, {XYList[[bin,1]],XYList[[bin,2]],phi, p, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff, yOff, rA,steps}],
	   {th0,0,thetamax[rF]},
	   {phi,0,2*Pi},
	   {
	   	p,
	   	pmin2DCases[{XYList[[bin,1]], phi, th0, BRxB, alpha, rA, rRxB,rD}, {xA,xOff}],
	   	pmax2DCases[{XYList[[bin,1]], phi, th0, BRxB, alpha, rA, rRxB,rD}, {xA,xOff}]
	   	},
	   PrecisionGoal->IntPrec
	],
	{bin,1,BinN},Method->"FinestGrained"
]  


 Integration2DArith[b_?NumericQ, {alpha_, BRxB_, rRxB_, rD_,rF_}, {xA_, yA_, xOff_, yOff_, rA_},BinN_,XYList_,IntPrec_]:=
ParallelTable[
	NIntegrate[
	  
	   Integrand2DArith[b, {XYList[[bin,1]],XYList[[bin,2]],phi, p, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff, yOff, rA}],
	   
	   	{th0,0,thetamax[rF]},
	   {phi,0,2*Pi},
	   {p,
	   	pmin2DCases[{XYList[[bin,1]], phi, th0, BRxB, alpha, rA, rRxB,rD}, {xA,xOff}],
	   	pmax2DCases[{XYList[[bin,1]], phi, th0, BRxB, alpha, rA, rRxB,rD}, {xA,xOff}]
	   	},
	   PrecisionGoal->IntPrec
	],
	{bin,1,BinN},Method->"FinestGrained"
]

Integration2DArithIntMC[b_?NumericQ, {alpha_, BRxB_, rRxB_, rD_,rF_}, {xA_, yA_, xOff_, yOff_, rA_},BinN_,XYList_,IntPrec_]:=
ParallelTable[
	NIntegrate[
	  
	   Integrand2DArith[b, {XYList[[bin,1]],XYList[[bin,2]],phi, p, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff, yOff, rA}],
	   
	   	{th0,0,thetamax[rF]},
	   {phi,0,2*Pi},
	   {p,
	   	pmin2DCases[{XYList[[bin,1]], phi, th0, BRxB, alpha, rA, rRxB,rD}, {xA,xOff}],
	   	pmax2DCases[{XYList[[bin,1]], phi, th0, BRxB, alpha, rA, rRxB,rD}, {xA,xOff}]
	   	},
	   PrecisionGoal->IntPrec, Method -> {"AdaptiveMonteCarlo", "MaxPoints" -> 500000, 
  								Method -> {"MonteCarloRule", "AxisSelector" -> {"MinVariance", "SubsampleFraction" -> 1/2}}, 
  								"SymbolicProcessing" -> False}
	],
	{bin,1,BinN},Method->"FinestGrained"
]



(* ::Section:: *)
(* with Bgradient*)



IntTable2DArithBGrad[b_?NumericQ, {alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, rD_?NumericQ,rF_?NumericQ},
	 {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, yOff_?NumericQ, rA_?NumericQ},{R0_?NumericQ,G1_?NumericQ,G2_?NumericQ},
	 XYList_List,IntPrec_,kernels_]:=
Block[{intresult},
	LaunchKernels[kernels];
	intresult=ParallelTable[
	NIntegrate[
	  Re[
	   Integrand2DArithBGrad[b, {bin[[1]],bin[[2]],phi, p, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff, yOff, rA},{R0,G1,G2}]
	  ],
	   
	   {th0,0,thetamax[rF]},
	   {phi,0,2*Pi},
	   {p,
	   	pmin2DBGradCases[{bin[[1]],bin[[2]], phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff},{R0,G1,G2}],
	   	pmax2DBGradCases[{bin[[1]],bin[[2]], phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff},{R0,G1,G2}]
	   	},
	   PrecisionGoal->IntPrec,AccuracyGoal->5,MaxRecursion->5,MaxPoints->300000
	],
	{bin,XYList},Method->"FinestGrained"
	];
	CloseKernels[];
	intresult
]



(* ::Section:: *)
(* with Bgradient and Det spread *)

IntTable2DArithBGradDetSpread[b_?NumericQ, {alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, rD_?NumericQ,rF_?NumericQ},
	 {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, yOff_?NumericQ, rA_?NumericQ},{R0_?NumericQ,G1_?NumericQ,G2_?NumericQ},
	 XYList_List,IntPrec_,kernels_]:=
Block[{intresult},
	LaunchKernels[kernels];
	intresult=ParallelTable[
	Quiet[NIntegrate[
	  Re[
	   Integrand2DArithBGradDetSpread[b, {bin[[1]],bin[[2]],phi, p, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff, yOff, rA},{R0,G1,G2}]
	   ],
	   
	   {th0,0,thetamax[rF]},
	   {phi,0,2*Pi},
	   {p,
	   	pmin2DBGradDetSpreadCases[{bin[[1]],bin[[2]], phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff},{R0,G1,G2}],
	   	pmax2DBGradDetSpreadCases[{bin[[1]],bin[[2]], phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff},{R0,G1,G2}]
	   	},
	   PrecisionGoal->IntPrec,AccuracyGoal->4,MaxRecursion->5,MaxPoints->80000,Method -> {"AdaptiveMonteCarlo", "MaxPoints" -> 500000, 
  								Method -> {"MonteCarloRule", "AxisSelector" -> {"MinVariance", "SubsampleFraction" -> 1/2}}, 
  								"SymbolicProcessing" -> False}
	],{NIntegrate::slwcon}],
	{bin,XYList},Method->"FinestGrained"
	];
	CloseKernels[];
	intresult
]



(* ::Section:: *)
(* with Bgradient and Det spread inc. phi restrictions *)


IntTable2DArithBGradDetSpreadPhiLimits[b_?NumericQ, {alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, rD_?NumericQ,rF_?NumericQ},
	 {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, yOff_?NumericQ, rA_?NumericQ},{R0_?NumericQ,G1_?NumericQ,G2_?NumericQ},
	 XYList_List,IntPrec_,kernels_]:=
Block[{intresult},
	LaunchKernels[kernels];
	intresult=ParallelTable[
	NIntegrate[
	  Re[
	   Integrand2DArithBGradDetSpread[b, {bin[[1]],bin[[2]],phi, p, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff, yOff, rA},{R0,G1,G2}]
	   ],
	   
	   {th0,0,thetamax[rF]},
	   {phi, PhiLimits[bin[[2]], th0, BRxB, rA, rRxB, yA, yOff, rD][[1]],Sequence@@PhiLimits[bin[[2]], th0, BRxB, rA, rRxB, yA, yOff, rD][[2;;-1]]},
	   {p,
	   	pmin2DBGradDetSpreadCases[{bin[[1]],bin[[2]], phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff},{R0,G1,G2}],
	   	pmax2DBGradDetSpreadCases[{bin[[1]],bin[[2]], phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff},{R0,G1,G2}]
	   	},
	   PrecisionGoal->IntPrec,AccuracyGoal->4,MaxRecursion->5,MaxPoints->80000
	],
	{bin,XYList},Method->"FinestGrained"
	];
	CloseKernels[];
	intresult
]







