(* Wolfram Language package *)


Get["Spectra/ElectronSpectrum.m"];
Get["Aperture/ApertureDefinition.m"];
Get["Bfield+Drift/Bfields.m"]
Get["Bfield+Drift/DetectorGyrationWeighting.m"];
Get["Bfield+Drift/Drifts.m"];
Get["Common/CommonFunctions.m"];
Get["Aperture/ArithApert.m"];


(* ::Section:: *)
(* Simple Drift *)

Integrand2D[
  b_?NumericQ, {xD_,yD_,phi_?NumericQ, p_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rD_}, {xA_, yA_, xOff_, 
   yOff_, rA_,prec_}] :=
 Re[
 	
   wmomNormedWb[\[Lambda]0,\[Kappa]0,b, p]*
    Sin[th0]*
    ApertureFunc[
    	xA, yA, xOff, yOff, 
    	xD + rG[p,theta2[th0,rD],rD*BRxB/rRxB]*Sin[phi] + D1stSimple[p, alpha, BRxB, th0, rRxB],
    	yD + rG[p,theta2[th0,rD],rD*BRxB/rRxB]*Cos[phi],
    	theta2[th0, rA], p, BRxB,prec
    	 ]
    	 
   ]
   
   
pmin2D[{xD_, phi_, th0_, BRxB_, alpha_, rA_, rRxB_, rD_}, {xA_,xOff_}] = 
 Solve[
   
   -xA/2 - rG[pmin, th0, rA*BRxB/rRxB] + xOff == 
    xD + rG[pmin, th0, rD*BRxB/rRxB]*Sin[phi] + 
     D1stSimple[pmin, alpha, BRxB, th0, rRxB], pmin][[1, 1, 2]]
   
pmin2DCases[{xD_, phi_, th0_, BRxB_, alpha_, rA_, rRxB_, 
   rD_}, {xA_,xOff_}] := Piecewise[
  {
   {0., pmin2D[{xD, phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff}] < 0},
   {pmax, 
    pmin2D[{xD, phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff}] > pmax}
   },
  pmin2D[{xD, phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff}]
  ]
  
  
pmax2D[{xD_, phi_, th0_, BRxB_, alpha_, rA_, rRxB_, rD_}, {xA_,xOff_}] = 
 Solve[
   
   xA/2 + rG[pmin, th0, rA*BRxB/rRxB] +xOff == 
    xD + rG[pmin, th0, rD*BRxB/rRxB]*Sin[phi] + 
     D1stSimple[pmin, alpha, BRxB, th0, rRxB], pmin][[1, 1, 2]]
     
pmax2DCases[{xD_, phi_, th0_, BRxB_, alpha_, rA_, rRxB_, 
   rD_}, {xA_,xOff_}] := Piecewise[
  {
   {0., pmax2D[{xD, phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff}] < 0},
   {pmax, 
    pmax2D[{xD, phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff}] > pmax}
   },
  pmax2D[{xD, phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff}]
  ]   
  

  
  
  Integrand2DApertMC[
  b_?NumericQ, {xD_?NumericQ,yD_?NumericQ,phi_?NumericQ, p_?NumericQ, th0_?NumericQ, alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, rD_?NumericQ},
   {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, yOff_?NumericQ, rA_?NumericQ,steps_?NumericQ}] :=
 Re[
 	
   wmomNormedWb[\[Lambda]0,\[Kappa]0,b, p]*
    Sin[th0]*
    MonteCarloAperture[
    	xA, yA, xOff, yOff, 
    	xD + rG[p,theta2[th0,rD],rD*BRxB/rRxB]*Sin[phi] + D1stSimple[p, alpha, BRxB, th0, rRxB],
    	yD + rG[p,theta2[th0,rD],rD*BRxB/rRxB]*Cos[phi],
    	theta2[th0, rA], p, rA*BRxB/rRxB,steps
    ]
    	 
   ] 

   
Integrand2DArith[
  b_?NumericQ, {xD_,yD_,phi_?NumericQ, p_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rD_}, {xA_, yA_, xOff_, yOff_, rA_}] :=
 	
   wmomNormedWb[\[Lambda]0,\[Kappa]0,b, p]*
    Sin[th0]*
    ArithApert[
    	xA, yA, xOff, yOff, 
    	xD + rG[p,theta2[th0,rD],rD*BRxB/rRxB]*Sin[phi] + D1stSimple[p, alpha, BRxB, th0, rRxB],
    	yD + rG[p,theta2[th0,rD],rD*BRxB/rRxB]*Cos[phi],
    	rG[p, theta2[th0, rA], BRxB]
]
   
   
     
  (* ::Section:: *)
(* Integrand with B gradient *)
   
   
   Integrand2DBGrad[
  b_, {xD_,yD_,phi_, p_, th0_, alpha_, BRxB_, rRxB_, rD_,R0_,G1_,G2_}, 
  	{xA_, yA_, xOff_, yOff_, rA_,prec_}] :=
 
 Re[
   wmomNormedWb[\[Lambda]0,\[Kappa]0,b, p]*
    Sin[th0]*
    ApertureFunc[
    	xA, yA, xOff, yOff, 
    	xD + rG[p,theta2[th0,rD],rD*BRxB/rRxB]*Sin[phi] 
    			+ D1stBPolyGrad[p, alpha, th0, BRxB, rRxB, yD + rG[p,theta2[th0,rD],rD*BRxB/rRxB]*Cos[phi] , R0, G1, G2],
    	yD + rG[p,theta2[th0,rD],rD*BRxB/rRxB]*Cos[phi],
    	theta2[th0, rA], p, BRxB,prec
    	 ]   	 
   ]

pmin2DBPoly[{xD_, phi_, th0_, BRxB_, alpha_, rA_, rRxB_, rD_}, {xA_}] = 
 Solve[
   
   -xA/2 - rG[pmin, th0, rA*BRxB/rRxB] == 
    xD + rG[pmin, th0, rD*BRxB/rRxB]*Sin[phi] + 
     D1stBPolyGrad[pmin, alpha, th0, BRxB, rRxB, yD + rG[pmin,theta2[th0,rD],rD*BRxB/rRxB]*Cos[phi] , R0, G1, G2]
     , pmin][[1, 1, 2]]
   
  
  