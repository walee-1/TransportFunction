(* Wolfram Language package *)

Get["TransportFunctions/Integrands.m"];


IntegrationTableArith[{b_, alpha_, BRxB_, rRxB_, rD_, rF_, 
    rA_}, {xA_, yA_, xOff_, yOff_, y0_},XList_,IntPrec_] := ParallelTable[
   NIntegrate[
   Re[
    IntegrandArith[b, {xi, p, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff,  yOff, x0, y0, rA}]
   ],
    
    	{th0, 0, thetamax[rF]},
    
    	{x0, -xA/2 - rG[pmax, theta2[th0,rA],rA* BRxB/rRxB], xA/2 + rG[pmax, theta2[th0,rA], rA*BRxB/rRxB]},
    
    	{
     		p,
    		 pminCases[x0, xi, theta2[th0, rD], rD*BRxB/rRxB, alpha, th0, rRxB, BRxB],
    		 pmaxCases[x0, xi, theta2[th0, rD], rD*BRxB/rRxB, alpha, th0, rRxB, BRxB]
    	 },
     
    PrecisionGoal -> IntPrec, AccuracyGoal->6,MaxRecursion->5,MaxPoints->120000(*, Method -> {"AdaptiveMonteCarlo", "MaxPoints" -> 100000, 
    	Method -> {"MonteCarloRule", "AxisSelector" -> {"MinVariance", "SubsampleFraction" -> 1/2}}}*)
    ],
   {xi, XList}, Method -> "FinestGrained"];



ProtonIntTableArith[{a_, alpha_, BRxB_, rRxB_, rD_, rF_, 
    rA_}, {xA_, yA_, xOff_, yOff_, y0_},XList_,IntPrec_] := ParallelTable[
   NIntegrate[
   Re[
    ProtonIntegrandArith[a, {xi, p, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff,  yOff, x0, y0, rA}]
   ],
    
    	{th0, 0, thetamax[rF]},
    
    	{x0, -xA/2 - rG[pmax, theta2[th0,rA],rA* BRxB/rRxB], xA/2 + rG[pmax, theta2[th0,rA], rA*BRxB/rRxB]},
    
    	{
     		p,
    		 pminCases[x0, xi, theta2[th0, rD], rD*BRxB/rRxB, alpha, th0, rRxB, BRxB],
    		 pmaxCases[x0, xi, theta2[th0, rD], rD*BRxB/rRxB, alpha, th0, rRxB, BRxB]
    	 },
     
    PrecisionGoal -> IntPrec, AccuracyGoal->6,MaxRecursion->5,MaxPoints->120000(*, Method -> {"AdaptiveMonteCarlo", "MaxPoints" -> 100000, 
    	Method -> {"MonteCarloRule", "AxisSelector" -> {"MinVariance", "SubsampleFraction" -> 1/2}}}*)
    ],
   {xi, XList}, Method -> "FinestGrained"];
   
   
