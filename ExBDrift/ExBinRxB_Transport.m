(* ::Package:: *)

(* Wolfram Language package *)

Get["Common/Constants.m"];
Get["Spectra/ElectronSpectrum.m"];
Get["Aperture/ApertureDefinition.m"];
Get["Bfield+Drift/Bfields.m"]
Get["Bfield+Drift/DetectorGyrationWeighting.m"];
Get["Bfield+Drift/Drifts.m"];
Get["Common/CommonFunctions.m"];
Get["Aperture/ArithApert.m"];
Get["AddedAngle_Alpha/Definitions.m"];
Get["2DTransport/2DIntegrand.m"];



DExBinRxB[p_,alpha_,R_,BRxB_,ExBinRxBScale_]:=ExBinRxBScale*(0.5+1.94*p/pmax)*alpha*R*mp/p/BRxB/c


Integrand2DArithExBinRxB[
  b_?NumericQ, {xD_,yD_,phi_?NumericQ, p_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rD_}, {xA_, yA_, xOff_, yOff_, rA_},
  {R_,ExBinRxBScale_}] :=
 	
   wmomNormedWb[\[Lambda]0,\[Kappa]0,b, p]*
    Sin[th0]*
    ArithApert[
    	xA, yA, xOff, yOff, 
    	xD + rG[p,theta2[th0,rD],rD*BRxB/rRxB]*Sin[phi] + D1stSimple[p, alpha, BRxB, th0, rRxB],
    	yD + rG[p,theta2[th0,rD],rD*BRxB/rRxB]*Cos[phi] + DExBinRxB[p,alpha,R,BRxB,ExBinRxBScale],
    	rG[p, theta2[th0, rA], BRxB]
]


(*the same pmin and pmax still work, because the argument for x0 is the same, only y0 changes, which would result in different additional conditions*)

 IntTable2DArithExBinRxB[b_?NumericQ, {alpha_, BRxB_, rRxB_, rD_,rF_}, {xA_, yA_, xOff_, yOff_, rA_},{R_,ExBinRxBScale_},BinN_,XYList_,IntPrec_]:=
ParallelTable[
	NIntegrate[
	  
	   Integrand2DArithExBinRxB[b, {XYList[[bin,1]],XYList[[bin,2]],phi, p, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff, yOff, rA},{R,ExBinRxBScale}],
	   
	   {th0,0,thetamax[rF]},
	   {phi,0,2*Pi},
	   {p,
	   	pmin2DCases[{XYList[[bin,1]], phi, th0, BRxB, alpha, rA, rRxB,rD}, {xA,xOff}],
	   	pmax2DCases[{XYList[[bin,1]], phi, th0, BRxB, alpha, rA, rRxB,rD}, {xA,xOff}]
	   	},
	   PrecisionGoal->IntPrec,AccuracyGoal->5,MaxRecursion->5,MaxPoints->200000
	],
	{bin,1,BinN},Method->"FinestGrained"
]

