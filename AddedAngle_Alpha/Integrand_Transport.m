(* Wolfram Language package *)


Get["Spectra/ElectronSpectrum.m"];
Get["Aperture/ApertureDefinition.m"];
Get["Bfield+Drift/Bfields.m"]
Get["Bfield+Drift/DetectorGyrationWeighting.m"];
Get["Bfield+Drift/Drifts.m"];
Get["Common/CommonFunctions.m"];
Get["Aperture/ArithApert.m"];
Get["AddedAngle_Alpha/Definitions.m"];


(* ::Section:: *)
(*Integrand with only Transition region -> different alpha behaviour*)

 IntegrandArithAlphaPrime[
  b_?NumericQ, {xi_?NumericQ, p_?NumericQ, th0_?NumericQ, alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, rD_?NumericQ}, {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, 
   yOff_?NumericQ, x0_?NumericQ, y0_?NumericQ, rA_?NumericQ},{R_,AlphaTrans_,zt_}] :=
    
    ArithApert[xA, yA, xOff, yOff, x0, y0, rG[ p,theta2[th0, rA], rA*BRxB/rRxB]]*
     wmomNormedWb[\[Lambda]0,\[Kappa]0,b, p]*
     Sin[th0]*
     WeightrG[xi - x0 + DSumAlphaPrime[p, th0, alpha, BRxB, 
      rRxB, {R, AlphaTrans, zt}], p, theta2[th0, rD], rD*BRxB/rRxB]/(xA + 2*rG[pmax, th0, rA*BRxB/rRxB])
     
     (*now, pmin , pmax changes!*)
        
pminFuncAlphaPrime[x0_, xD_, th0_, alpha_, rRxB_, BRxB_, 
  rD_, {R_, AlphaTrans_, zt_}] = 
 Solve[xD - x0 + 
     DSumAlphaPrime[p, th0, alpha, BRxB, 
      rRxB, {R, AlphaTrans, zt}] == -rG[p, theta2[th0, rD], 
      rD*BRxB/rRxB], p][[1, 1, 2]]
      
pmaxFuncAlphaPrime[x0_, xD_, th0_, alpha_, rRxB_, BRxB_, 
  rD_, {R_, AlphaTrans_, zt_}] = 
 Solve[xD - x0 + 
     DSumAlphaPrime[p, th0, alpha, BRxB, rRxB, {R, AlphaTrans, zt}] ==
     rG[p, theta2[th0, rD], rD*BRxB/rRxB], p][[1, 1, 2]]
     
             
pminAlphaPrimeCases[x0_, xD_, th0_, alpha_, rRxB_, BRxB_, rD_, {R_, AlphaTrans_, zt_}] := Which[
    pminFuncAlphaPrime[x0, xD, th0, alpha, rRxB, BRxB, rD, {R, AlphaTrans, zt}] > pmax, pmax,
   pminFuncAlphaPrime[x0, xD, th0, alpha, rRxB, BRxB, rD, {R, AlphaTrans, zt}] < 0, 0,
    True, pminFuncAlphaPrime[x0, xD, th0, alpha, rRxB, BRxB, rD, {R, AlphaTrans, zt}]
    ]
    
pmaxAlphaPrimeCases[x0_, xD_, th0_, alpha_, rRxB_, BRxB_, rD_, {R_, AlphaTrans_, zt_}] := Which[
    pmaxFuncAlphaPrime[x0, xD, th0, alpha, rRxB, BRxB, rD, {R, AlphaTrans, zt}] > pmax, pmax,
   pmaxFuncAlphaPrime[x0, xD, th0, alpha, rRxB, BRxB, rD, {R, AlphaTrans, zt}] < 0, 0,
    True, pmaxFuncAlphaPrime[x0, xD, th0, alpha, rRxB, BRxB, rD, {R, AlphaTrans, zt}]
    ]
    
    
IntTableArithAlphaPrime[b_?NumericQ, {alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, rD_?NumericQ,rF_}, {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, 
   yOff_?NumericQ, y0_?NumericQ, rA_?NumericQ},{R_,AlphaTrans_,zt_},XList_,IntPrec_] := ParallelTable[
   NIntegrate[
   Re[
    IntegrandArithAlphaPrime[b, {xD, p, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff, yOff, x0, y0, rA},{R,AlphaTrans,zt}]
   ],
    
    	{th0, 0, thetamax[rF]},
    
    	{x0, -xA/2 - rG[pmax, th0, BRxB]+ xOff, xA/2 + rG[pmax, th0, BRxB]+ xOff},
    
    	{
     		p,
    		 pminAlphaPrimeCases[x0, xD, th0, alpha, rRxB, BRxB, rD, {R, AlphaTrans, zt}],
    		 pmaxAlphaPrimeCases[x0, xD, th0, alpha, rRxB, BRxB, rD, {R, AlphaTrans, zt}]
    	 },
     
    PrecisionGoal -> IntPrec, AccuracyGoal->5,MaxPoints->150000,MaxRecursion->5(* Method -> {"AdaptiveMonteCarlo", "MaxPoints" -> 400000, 
    	Method -> {"MonteCarloRule", "AxisSelector" -> {"MinVariance", "SubsampleFraction" -> 1/2}}}*)
    ],
   {xD, XList}, Method -> "FinestGrained"];   
    
    
    
    
    
    
    (* ::Section:: *)
(*Added angle in addition*)

IntegrandArithAddedAngleAlphaPrime[
  b_?NumericQ, {xi_?NumericQ, p_?NumericQ, th0_?NumericQ, alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, rD_?NumericQ},
   {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, yOff_?NumericQ, x0_?NumericQ, y0_?NumericQ, rA_?NumericQ},
   	{R_,AlphaTrans_,zt_,AddedScale_}] :=
    
    ArithApert[xA, yA, xOff, yOff, x0, y0, rG[ p,theta2[th0, rA], rA*BRxB/rRxB]]*
     wmomNormedWb[\[Lambda]0,\[Kappa]0,b, p]*
     Sin[th0]*
     WeightrG[xi - x0 + DSumAddedAngleAlpha[p,th0,alpha,BRxB,rRxB,{R,AlphaTrans,zt,AddedScale}], p, theta2[th0, rD], rD*BRxB/rRxB]/
     (xA + 2*rG[pmax, theta2[th0,rA], rA*BRxB/rRxB])
     
     (*now, pmin , pmax changes!*)
     
(* Solve can not analytically solve -> findroot *)

     
pminFuncAddedAngleAlphaPrime[x0_, xD_, th0_, alpha_, rRxB_, BRxB_, 
  rD_, {R_, AlphaTrans_, zt_,AddedScale_}] := 
 Re[FindRoot[xD - x0 + 
     DSumAddedAngleAlpha[p,th0,alpha,BRxB,rRxB,{R,AlphaTrans,zt,AddedScale}] == -rG[p, theta2[th0, rD], rD*BRxB/rRxB], {p,500000}][[1, 2]]]
      
pmaxFuncAddedAngleAlphaPrime[x0_, xD_, th0_, alpha_, rRxB_, BRxB_, 
  rD_, {R_, AlphaTrans_, zt_,AddedScale_}] := 
 Re[FindRoot[xD - x0 + 
     DSumAddedAngleAlpha[p,th0,alpha,BRxB,rRxB,{R,AlphaTrans,zt,AddedScale}] == rG[p, theta2[th0, rD], rD*BRxB/rRxB], {p,500000}][[1, 2]]]
     
          
      
pminAddedAngleAlphaPrimeCases[x0_?NumericQ, xD_?NumericQ, th0_?NumericQ, alpha_, rRxB_, BRxB_, rD_, {R_, AlphaTrans_, zt_,AddedScale_}] := Which[
    pminFuncAddedAngleAlphaPrime[x0, xD, th0, alpha, rRxB, BRxB, rD, {R, AlphaTrans, zt,AddedScale}] > pmax, pmax,
   pminFuncAddedAngleAlphaPrime[x0, xD, th0, alpha, rRxB, BRxB, rD, {R, AlphaTrans, zt,AddedScale}] < 0, 0,
    True, pminFuncAddedAngleAlphaPrime[x0, xD, th0, alpha, rRxB, BRxB, rD, {R, AlphaTrans, zt,AddedScale}]
    ]
    
pmaxAddedAngleAlphaPrimeCases[x0_?NumericQ, xD_?NumericQ, th0_?NumericQ, alpha_, rRxB_, BRxB_, rD_, {R_, AlphaTrans_, zt_,AddedScale_}] := Which[
    pmaxFuncAddedAngleAlphaPrime[x0, xD, th0, alpha, rRxB, BRxB, rD, {R, AlphaTrans, zt,AddedScale}] > pmax, pmax,
   pmaxFuncAddedAngleAlphaPrime[x0, xD, th0, alpha, rRxB, BRxB, rD, {R, AlphaTrans, zt,AddedScale}] < 0, 0,
    True, pmaxFuncAddedAngleAlphaPrime[x0, xD, th0, alpha, rRxB, BRxB, rD, {R, AlphaTrans, zt,AddedScale}]
    ]


IntTableArithAddedAngleAlphaPrime[b_?NumericQ, {alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, rD_?NumericQ,rF_}, {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, 
   yOff_?NumericQ, y0_?NumericQ, rA_?NumericQ},{R_,AlphaTrans_,zt_,AddedScale_},XList_,IntPrec_] := ParallelTable[
   NIntegrate[
   Re[
    IntegrandArithAddedAngleAlphaPrime[b, {xD, p, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff, yOff, x0, y0, rA},{R,AlphaTrans,zt,AddedScale}]
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



    