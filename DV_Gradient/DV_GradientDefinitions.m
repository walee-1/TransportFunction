(* Wolfram Language package *)

Get["Spectra/ElectronSpectrum.m"];
Get["Aperture/ApertureDefinition.m"];
Get["Bfield+Drift/Bfields.m"]
Get["Bfield+Drift/DetectorGyrationWeighting.m"];
Get["Bfield+Drift/Drifts.m"];
Get["Common/CommonFunctions.m"];
Get["Aperture/ArithApert.m"];
Get["Spectra/ProtonSpectrumNachtmann.m"];



BRatioDVGrad[ri0_,kDV_,z0_?NumericQ]:=ri0/(1+kDV*z0)


IntegrandArithDVGrad[
  b_?NumericQ, {xi_?NumericQ, p_?NumericQ, th0_?NumericQ, alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, rD_?NumericQ},
   {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, yOff_?NumericQ, x0_?NumericQ, y0_?NumericQ, rA_?NumericQ},
   {kDV_?NumericQ,z0_?NumericQ}] :=
    
    ArithApert[xA, yA, xOff, yOff, x0, y0, rG[ p,theta2[th0, BRatioDVGrad[rA,kDV,z0]],(*this doesnt change for DVGrad*) rA*BRxB/rRxB]]*
     wmomNormedWb[\[Lambda]0,\[Kappa]0,b, p]*
     Sin[th0]*
     WeightrG[xi - x0 + D1stSimple[p, alpha, BRxB, th0, BRatioDVGrad[rRxB,kDV,z0]], p, theta2[th0, BRatioDVGrad[rD,kDV,z0]], rD*BRxB/rRxB]/
     (xA + 2*rG[pmax, theta2[th0,BRatioDVGrad[rA,kDV,z0]], rA*BRxB/rRxB])
     
     
     (*for the integration, we don't actually have to recalc pmin pmax, because the factor just follows the ratios through the calculation
     and can be applied after the fact here*)
     
     
     (*but we do it because we suspect the evaluation of BRatioDVGrad not working, resulting in pminCases not working!*)
     
pminCases2[x0_?NumericQ, x_?NumericQ, alpha_?NumericQ, th0_?NumericQ, rRxB_?NumericQ,rD_?NumericQ, BRxB_?NumericQ,kDV_?NumericQ,z0_?NumericQ] := 
     pminCases[x0, x, theta2[th0,  BRatioDVGrad[rD,kDV,z0]], rD*BRxB/rRxB, alpha, th0, BRatioDVGrad[rRxB,kDV,z0], BRxB]
pmaxCases2[x0_?NumericQ, x_?NumericQ, alpha_?NumericQ, th0_?NumericQ, rRxB_?NumericQ,rD_?NumericQ, BRxB_?NumericQ,kDV_?NumericQ,z0_?NumericQ] := 
     pmaxCases[x0, x,  theta2[th0,  BRatioDVGrad[rD,kDV,z0]], rD*BRxB/rRxB, alpha, th0, BRatioDVGrad[rRxB,kDV,z0], BRxB]
     
     
rG2[th0_?NumericQ,rA_?NumericQ,BRxB_,rRxB_,kDV_?NumericQ,z0_?NumericQ]:=rG[pmax, theta2[th0, BRatioDVGrad[rA,kDV,z0]],rA* BRxB/rRxB] 

thetamax2[rF_?NumericQ,kDV_?NumericQ,z0_?NumericQ]:=Re[thetamax[BRatioDVGrad[rF,kDV,z0]]]
     
     
IntTableArithDVGrad[{b_, alpha_, BRxB_, rRxB_, rD_, rF_, rA_}, 
	{xA_, yA_, xOff_, yOff_, y0_},{kDV_,LDV_},XList_,IntPrec_] := ParallelTable[
   NIntegrate[
   Re[
    IntegrandArithDVGrad[b, {xi, p, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff,  yOff, x0, y0, rA},{kDV,z0}]
   ],
    
    	{z0,-LDV/2,LDV/2},
    	{th0, 0, thetamax2[rF,kDV,z0]},
    
    	{x0, -xA/2 - rG2[th0, rA, BRxB, rRxB, kDV, z0], xA/2 + rG2[th0, rA, BRxB, rRxB, kDV, z0]},
    
    	{
     		p,
    		 pminCases2[x0, xi, alpha, th0, rRxB,rD, BRxB,kDV,z0](*here we already implement the transformation of ratios inside*),
    		 pmaxCases2[x0, xi, alpha, th0, rRxB,rD, BRxB,kDV,z0]
    	 },
     
    PrecisionGoal -> IntPrec, AccuracyGoal->5,MaxRecursion->5,MaxPoints->180000(*, Method -> {"AdaptiveMonteCarlo", "MaxPoints" -> 100000, 
    	Method -> {"MonteCarloRule", "AxisSelector" -> {"MinVariance", "SubsampleFraction" -> 1/2}}}*)
    ],
   {xi, XList}, Method -> "FinestGrained"];
   
   
   
   