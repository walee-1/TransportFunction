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
    	rG[p, theta2[th0, rA], rA*BRxB/rRxB]
]
   
   
   
   
     
  (* ::Section:: *)
(* Integrand with B gradient *)
  

   
Integrand2DArithBGrad[
  b_?NumericQ, {xD_?NumericQ,yD_?NumericQ,phi_?NumericQ, p_?NumericQ, th0_?NumericQ, alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, rD_?NumericQ},
   {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, yOff_?NumericQ, rA_?NumericQ},{R0_?NumericQ,G1_?NumericQ,G2_?NumericQ}] :=
 	
   wmomNormedWb[\[Lambda]0,\[Kappa]0,b, p]*
    Sin[th0]*
    ArithApert[
    	xA, yA, xOff, yOff, 
    	xD + rG[p,theta2[th0,rD],rD*BRxB/rRxB]*Sin[phi] 
					+ D1stBPolyGrad[p, alpha, th0, BRxB, rRxB, yD + rG[p,theta2[th0,rD],rD*BRxB/rRxB]*Cos[phi] , R0, G1, G2],
    	yD + rG[p,theta2[th0,rD],rD*BRxB/rRxB]*Cos[phi],
    	rG[p, theta2[th0, rA], rA*BRxB/rRxB]
]

   

pmin2DBPoly[{xD_,yD_, phi_, th0_, BRxB_, alpha_, rA_, rRxB_, rD_}, {xA_,xOff_},{R0_,G1_,G2_}] = 
 Solve[
   
   -xA/2 - rG[pmin, th0, rA*BRxB/rRxB] + xOff == 
    xD + rG[pmin, theta2[th0,rD], rD*BRxB/rRxB]*Sin[phi] + 
     D1stBPolyGrad[pmin, alpha, th0, BRxB, rRxB, yD + rG[pmin,theta2[th0,rD],rD*BRxB/rRxB]*Cos[phi] , R0, G1, G2]
     , pmin][[1, 1, 2]];
     
pmax2DBPoly[{xD_,yD_, phi_, th0_, BRxB_, alpha_, rA_, rRxB_, rD_}, {xA_,xOff_},{R0_,G1_,G2_}] = 
 Solve[
   
   xA/2 + rG[pmin, th0, rA*BRxB/rRxB] + xOff == 
    xD + rG[pmin, theta2[th0,rD], rD*BRxB/rRxB]*Sin[phi] + 
     D1stBPolyGrad[pmin, alpha, th0, BRxB, rRxB, yD + rG[pmin,theta2[th0,rD],rD*BRxB/rRxB]*Cos[phi] , R0, G1, G2]
     , pmin][[3, 1, 2]];

pmin2DBGradCases[{xD_?NumericQ,yD_?NumericQ, phi_?NumericQ, th0_?NumericQ, BRxB_?NumericQ, alpha_?NumericQ, rA_?NumericQ, rRxB_?NumericQ, rD_?NumericQ},
	 {xA_?NumericQ,xOff_?NumericQ},{R0_?NumericQ,G1_?NumericQ,G2_?NumericQ}] := Piecewise[
  {
   {0., Re[pmin2DBPoly[{xD,yD, phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff},{R0,G1,G2}]] < 0},
   {pmax,Re[pmin2DBPoly[{xD,yD, phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff},{R0,G1,G2}]] > pmax}
   },
  Re[pmin2DBPoly[{xD,yD, phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff},{R0,G1,G2}]]
  ]
  
  
pmax2DBGradCases[{xD_?NumericQ,yD_?NumericQ, phi_?NumericQ, th0_?NumericQ, BRxB_?NumericQ, alpha_?NumericQ, rA_?NumericQ, rRxB_?NumericQ, rD_?NumericQ},
	 {xA_?NumericQ,xOff_?NumericQ},{R0_?NumericQ,G1_?NumericQ,G2_?NumericQ}] := Piecewise[
  {
   {0., Re[pmax2DBPoly[{xD,yD, phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff},{R0,G1,G2}]] < 0},
   {pmax,Re[pmax2DBPoly[{xD,yD, phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff},{R0,G1,G2}]] > pmax}
   },
  Re[pmax2DBPoly[{xD,yD, phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff},{R0,G1,G2}]]
  ]
  
  
  
  
  
(* ::Section:: *)
(* 2D integrand with BGrad and det spread*)



Integrand2DArithBGradDetSpread[
  b_?NumericQ, {xD_?NumericQ,yD_?NumericQ,phi_?NumericQ, p_?NumericQ, th0_?NumericQ, alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, rD_?NumericQ},
   {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, yOff_?NumericQ, rA_?NumericQ},{R0_?NumericQ,G1_?NumericQ,G2_?NumericQ}] :=
 	
   wmomNormedWb[\[Lambda]0,\[Kappa]0,b, p]*
    Sin[th0]*
    ArithApert[
    	xA, yA, xOff, yOff, 
    	(xD + rG[p,theta2[th0,rD],rD*BRxB/rRxB]*Sin[phi])*Sqrt[rD/rRxB]
					+ D1stBPolyGrad[p, alpha, th0, BRxB, rRxB, (yD + rG[p,theta2[th0,rD],rD*BRxB/rRxB]*Cos[phi])*Sqrt[rD/rRxB] , R0, G1, G2],
    	(yD + rG[p,theta2[th0,rD],rD*BRxB/rRxB]*Cos[phi])*Sqrt[rD/rRxB],
    	rG[p, theta2[th0, rA], rA*BRxB/rRxB]
]

(*we use NSolve here*)

pmin2DBGradDetSpreadNSolve[{xD_, yD_, phi_, th0_, BRxB_, alpha_, rA_, rRxB_, rD_}, {xA_, xOff_}, {R_, G1_, G2_}] := 
 NSolve[
 	-xA/2 - rG[pmin, theta2[th0, rA], rA*BRxB/rRxB] + xOff == (xD + rG[pmin, theta2[th0, rD], rD*BRxB/rRxB]*Sin[phi])*Sqrt[rD/rRxB] + 
     D1stBPolyGrad[pmin, alpha, th0, BRxB, rRxB, (yD + rG[pmin, theta2[th0, rD], rD*BRxB/rRxB]*Cos[phi])*Sqrt[rD/rRxB], R, G1, G2], 
     pmin, Reals][[1, 1, 2]]
     
pmax2DBGradDetSpreadNSolve[{xD_, yD_, phi_, th0_, BRxB_, alpha_, rA_, rRxB_, rD_}, {xA_, xOff_}, {R_, G1_, G2_}] := 
 NSolve[xA/2 + rG[pmin, theta2[th0, rA], rA*BRxB/rRxB] + xOff == (xD + rG[pmin, theta2[th0, rD], rD*BRxB/rRxB]*Sin[phi])*Sqrt[rD/rRxB] + 
     D1stBPolyGrad[pmin, alpha, th0, BRxB, rRxB, (yD + rG[pmin, theta2[th0, rD], rD*BRxB/rRxB]*Cos[phi])*Sqrt[rD/rRxB], R, G1, G2],
      pmin, Reals][[-1, 1, 2]]


pmin2DBGradDetSpreadCases[{xD_?NumericQ,yD_?NumericQ, phi_?NumericQ, th0_?NumericQ, BRxB_?NumericQ, alpha_?NumericQ, rA_?NumericQ, rRxB_?NumericQ, rD_?NumericQ},
	 {xA_?NumericQ,xOff_?NumericQ},{R_?NumericQ,G1_?NumericQ,G2_?NumericQ}] := Piecewise[
  {
   {0., pmin2DBGradDetSpreadNSolve[{xD,yD, phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff},{R,G1,G2}] < 0},
   {pmax,pmin2DBGradDetSpreadNSolve[{xD,yD, phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff},{R,G1,G2}] > pmax}
   },
  pmin2DBGradDetSpreadNSolve[{xD,yD, phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff},{R,G1,G2}]
  ]

pmax2DBGradDetSpreadCases[{xD_?NumericQ,yD_?NumericQ, phi_?NumericQ, th0_?NumericQ, BRxB_?NumericQ, alpha_?NumericQ, rA_?NumericQ, rRxB_?NumericQ, rD_?NumericQ},
	 {xA_?NumericQ,xOff_?NumericQ},{R_?NumericQ,G1_?NumericQ,G2_?NumericQ}] := Piecewise[
  {
   {0., pmax2DBGradDetSpreadNSolve[{xD,yD, phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff},{R,G1,G2}] < 0},
   {pmax,pmax2DBGradDetSpreadNSolve[{xD,yD, phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff},{R,G1,G2}] > pmax}
   },
  pmax2DBGradDetSpreadNSolve[{xD,yD, phi, th0, BRxB, alpha, rA, rRxB, rD}, {xA,xOff},{R,G1,G2}]
  ]



(*you can additionally have an estimation for phi limits, either to indicate a gap in the integration range or a min,max restriction*)

phiGapyDneg[yD_, th0_, BRxB_, rA_, rRxB_, yA_, yOff_, rD_] = 
  Solve[-yA/2 - rG[pmax, th0, rA*BRxB/rRxB] + 
      yOff == (yD + rG[pmax, th0, rD*BRxB/rRxB]*Cos[phi])*
      Sqrt[rD/rRxB], phi][[2, 1, 2]];
      
phiBoundyDpos[yD_, th0_, BRxB_, rA_, rRxB_, yA_, yOff_, rD_] = 
 Solve[yA/2 + rG[pmax, th0, rA*BRxB/rRxB] + 
     yOff == (yD + rG[pmax, th0, rD*BRxB/rRxB]*Cos[phi])*
     Sqrt[rD/rRxB], phi][[2, 1, 2]];

(* whether or not the integrand is separated in phi, is determined by yA and yOff! *)

PhiLimits[yD_?NumericQ, th0_?NumericQ, BRxB_, rA_, rRxB_, yA_, yOff_, rD_] :=
(Piecewise[
   {
    (*separated case*)
    (*{
    	{0, phiGapyDneg[yD, th0, BRxB, rA, rRxB, yA, yOff, rD], 2*Pi - phiGapyDneg[yD, th0, BRxB, rA, rRxB, yA, yOff, rD],2*Pi},
    	yD < -yA/2 + yOff
   	},*)
    (*bunched case*)
    {
    	{phiBoundyDpos[yD, th0, BRxB, rA, rRxB, yA, yOff, rD], 2*Pi - phiBoundyDpos[yD, th0, BRxB, rA, rRxB, yA, yOff, rD]}, 
     	yD > yA/2 + yOff
    }
   },
   (*when entire phase space*)
   {0, 2*Pi}
   ])


  