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
