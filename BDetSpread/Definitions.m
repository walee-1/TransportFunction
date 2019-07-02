(* Wolfram Language package *)

Get["Spectra/ElectronSpectrum.m"];
Get["Bfield+Drift/Bfields.m"]
Get["Bfield+Drift/DetectorGyrationWeighting.m"];
Get["Bfield+Drift/Drifts.m"];
Get["Common/CommonFunctions.m"];
Get["Aperture/ArithApert.m"];

(*the introduction of this systematic makes most sense for an offset aperture, so that spread is symmetric on spectrum
this is important to consider, when integrating over x0, which should go from xOff - xA/2 -rmax to xOff + xA/2 + rmax*)

 IntegrandArith[
  b_?NumericQ, {xi_?NumericQ, p_?NumericQ, th0_?NumericQ, alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, rD_?NumericQ}, 
  	{xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, yOff_?NumericQ, x0_?NumericQ, y0_?NumericQ, rA_?NumericQ}] :=
    
    ArithApert[xA, yA, xOff, yOff, x0, y0, rG[ p,theta2[th0, rA], rA*BRxB/rRxB]]*
     wmomNormedWb[\[Lambda]0,\[Kappa]0,b, p]*
     Sin[th0]*
     WeightrG[
     	xi - Sqrt[rRxB/rD]*(x0 - D1stSimple[p, alpha, BRxB, th0, rRxB]),
     	p, theta2[th0, rD], rD*BRxB/rRxB
     	]/(xA + 2*rG[pmax, theta2[th0,rA], rA*BRxB/rRxB])
     	
     	(*of course, also pmin and pmax change*)
     	
     	
pminFuncDetSpread[x0_, xD_, th0_, alpha_, rRxB_, BRxB_, rD_] = 
 Solve[xD - Sqrt[rRxB/rD]*(x0 - D1stSimple[p, alpha, BRxB, th0, rRxB]) == -rG[p, theta2[th0, rD], rD*BRxB/rRxB], p][[1, 1, 2]]
      
pmaxFuncDetSpread[x0_, xD_, th0_, alpha_, rRxB_, BRxB_, rD_] = 
 Solve[xD - Sqrt[rRxB/rD]*(x0 - D1stSimple[p, alpha, BRxB, th0, rRxB]) == rG[p, theta2[th0, rD], rD*BRxB/rRxB], p][[1, 1, 2]]
     
          
(*still to change!*) 

    
pminDetSpreadCases[x0_, xD_, th0_, alpha_, rRxB_, BRxB_, rD_] := Which[
    pminFuncDetSpread[x0, xD, th0, alpha, rRxB, BRxB, rD] > pmax, pmax,
   pminFuncDetSpread[x0, xD, th0, alpha, rRxB, BRxB, rD] < 0, 0,
    True, pminFuncDetSpread[x0, xD, th0, alpha, rRxB, BRxB, rD]
    ]
    
pmaxDetSpreadCases[x0_, xD_, th0_, alpha_, rRxB_, BRxB_, rD_] := Which[
    pmaxFuncDetSpread[x0, xD, th0, alpha, rRxB, BRxB, rD] > pmax, pmax,
   pmaxFuncDetSpread[x0, xD, th0, alpha, rRxB, BRxB, rD] < 0, 0,
    True, pmaxFuncDetSpread[x0, xD, th0, alpha, rRxB, BRxB, rD]
    ]
    
    
    

IntTableDetSpreadArith[{b_, alpha_, BRxB_, rRxB_, rD_, rF_, 
    rA_}, {xA_, yA_, xOff_, yOff_, y0_},XList_,IntPrec_] := ParallelTable[
   NIntegrate[
   Re[
    IntegrandArith[b, {xi, p, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff,  yOff, x0, y0, rA}]
   ],
    
    	{th0, 0, thetamax[rF]},
    
    	{x0, -xA/2 - rG[pmax, theta2[th0,rA],rA* BRxB/rRxB]+xOff, xA/2 + rG[pmax, theta2[th0,rA], rA*BRxB/rRxB]+xOff},
    
    	{
     		p,
    		 pminDetSpreadCases[x0, xi, th0, alpha,  rRxB, BRxB, rD],
    		 pmaxDetSpreadCases[x0, xi, th0, alpha,  rRxB, BRxB, rD]
    	 },
     
    PrecisionGoal -> IntPrec, AccuracyGoal->6,MaxRecursion->5,MaxPoints->120000(*, Method -> {"AdaptiveMonteCarlo", "MaxPoints" -> 100000, 
    	Method -> {"MonteCarloRule", "AxisSelector" -> {"MinVariance", "SubsampleFraction" -> 1/2}}}*)
    ],
   {xi, XList}, Method -> "FinestGrained"];   
    
    
    