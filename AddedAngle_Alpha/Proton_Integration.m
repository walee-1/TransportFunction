(* Wolfram Language package *)

Get["AddedAngle_Alpha/Integrand_Transport.m"];
Get["Spectra/ProtonSpectrumNachtmann.m"];



ProtonIntegrandArithAddedAngleAlphaPrime[
  a_?NumericQ, {xi_?NumericQ, p_?NumericQ, th0_?NumericQ, alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, rD_?NumericQ},
   {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, yOff_?NumericQ, x0_?NumericQ, y0_?NumericQ, rA_?NumericQ},
   	{R_,AlphaTrans_,zt_,AddedScale_}] :=
    
    ArithApert[xA, yA, xOff, yOff, x0, y0, rG[ p,theta2[th0, rA], rA*BRxB/rRxB]]*
     pmomNormed[ p,a]*
     Sin[th0]*
     WeightrG[xi - x0 + DSumAddedAngleAlpha[p,th0,alpha,BRxB,rRxB,{R,AlphaTrans,zt,AddedScale}], p, theta2[th0, rD], rD*BRxB/rRxB]/
     (xA + 2*rG[pmax, theta2[th0,rA], rA*BRxB/rRxB])
     
     
     
ProtonIntTableArithAddedAngleAlphaPrime[a_?NumericQ, {alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, rD_?NumericQ,rF_}, {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, 
   yOff_?NumericQ, y0_?NumericQ, rA_?NumericQ},{R_,AlphaTrans_,zt_,AddedScale_},XList_,IntPrec_] := ParallelTable[
   NIntegrate[
   Re[
    ProtonIntegrandArithAddedAngleAlphaPrime[a, {xD, p, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff, yOff, x0, y0, rA},{R,AlphaTrans,zt,AddedScale}]
   ],
    
    	{th0, 0, thetamax[rF]},
    
    	{x0, -xA/2 - rG[pmax, theta2[th0,rA], rA*BRxB/rRxB]+ xOff, xA/2 + rG[pmax, theta2[th0,rA], rA*BRxB/rRxB]+ xOff},
    
    	{
     		p,
    		 pminAddedAngleAlphaPrimeCases[x0, xD, th0, alpha, rRxB, BRxB, rD, {R, AlphaTrans, zt,AddedScale}],
    		 pmaxAddedAngleAlphaPrimeCases[x0, xD, th0, alpha, rRxB, BRxB, rD, {R, AlphaTrans, zt,AddedScale}]
    	 },
     
    PrecisionGoal -> IntPrec, AccuracyGoal->5,MaxPoints->150000,MaxRecursion->5(* Method -> {"AdaptiveMonteCarlo", "MaxPoints" -> 400000, 
    	Method -> {"MonteCarloRule", "AxisSelector" -> {"MinVariance", "SubsampleFraction" -> 1/2}}}*)
    ],
   {xD, XList}, Method -> "FinestGrained"];  