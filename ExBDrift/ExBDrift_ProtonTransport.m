(* Wolfram Language package *)

Get["Spectra/ElectronSpectrum.m"];
Get["Spectra/ProtonSpectrumNachtmann.m"];
Get["Aperture/ApertureDefinition.m"];
Get["Bfield+Drift/Bfields.m"]
Get["Bfield+Drift/DetectorGyrationWeighting.m"];
Get["Bfield+Drift/Drifts.m"];
Get["Common/CommonFunctions.m"];
Get["Aperture/ArithApert.m"];


k[p_, scale_] := scale*(20 - 6.5*p/pmax)/50


equ1 = x0 + De - kk*(y0 - r*Cos[phi]) - r*Sin[phi] == xD

equ2 = y0 + kk*(x0 + De - r*Sin[phi]) - r*Cos[phi] == yD

x0fromExB[kk_, xD_, yD_, r_, phi_, De_] = 
 Solve[equ1 && equ2, {x0, y0}][[1, 1, 2]];
 
 y0fromExB[kk_, xD_, yD_, r_, phi_] = 
 Solve[equ1 && equ2, {x0, y0}][[1, 2, 2]];


Integrand2DExB[
  a_?NumericQ, {xD_?NumericQ, yD_?NumericQ, phi_?NumericQ, 
   p_?NumericQ, th0_?NumericQ, alpha_?NumericQ, BRxB_?NumericQ, 
   rRxB_?NumericQ, rD_?NumericQ}, {xA_?NumericQ, yA_?NumericQ, 
   xOff_?NumericQ, yOff_?NumericQ, rA_?NumericQ, prec_?NumericQ}, 
  kscale_?NumericQ] :=
 Re[
  
  pmomNormed[p, a]*
   Sin[th0]*
   
   ApertureFunc[xA, yA, xOff, yOff,
    
    (*D1stSimple[p,alpha,BRxB,th0,rRxB]*)
    
    x0fromExB[k[p, kscale], xD, yD, rG[p, theta2[th0, rD], rD*BRxB/rRxB], phi, -D1stSimple[p, alpha, BRxB, th0, rRxB]],
    
    y0fromExB[k[p, kscale], xD, yD, rG[p, theta2[th0, rD], rD*BRxB/rRxB], phi],
    
    theta2[th0, rA], p, rA*BRxB/rRxB, prec]
  ]
  
  
  Integrand2DExBMC[
  a_?NumericQ, {xD_?NumericQ, yD_?NumericQ, phi_?NumericQ, 
   p_?NumericQ, th0_?NumericQ, alpha_?NumericQ, BRxB_?NumericQ, 
   rRxB_?NumericQ, rD_?NumericQ}, {xA_?NumericQ, yA_?NumericQ, 
   xOff_?NumericQ, yOff_?NumericQ, rA_?NumericQ, MCSteps_?IntegerQ}, 
  kscale_?NumericQ] :=
 Re[
  
  pmomNormed[p, a]*
   Sin[th0]*
   
   MonteCarloAperture[xA, yA, xOff, yOff,
    
    (*D1stSimple[p,alpha,BRxB,th0,rRxB]*)
    
    x0fromExB[k[p, kscale], xD, yD, rG[p, theta2[th0, rD], rD*BRxB/rRxB], phi, -D1stSimple[p, alpha, BRxB, th0, rRxB]],
    
    y0fromExB[k[p, kscale], xD, yD, rG[p, theta2[th0, rD], rD*BRxB/rRxB], phi],
    
    theta2[th0, rA], p, rA*BRxB/rRxB, MCSteps]
  ]
  
  Integrand2DExBArith[
  a_?NumericQ, {xD_?NumericQ, yD_?NumericQ, phi_?NumericQ, p_?NumericQ, th0_?NumericQ, alpha_?NumericQ, BRxB_?NumericQ, 
   rRxB_?NumericQ, rD_?NumericQ}, {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, yOff_?NumericQ, rA_?NumericQ}, 
  kscale_?NumericQ] :=
 Re[
  
  pmomNormed[p, a]*
   Sin[th0]*
   
   ArithApert[xA, yA, xOff, yOff,
    
    x0fromExB[k[p, kscale], xD, yD, rG[p, theta2[th0, rD], rD*BRxB/rRxB], phi, -D1stSimple[p, alpha, BRxB, th0, rRxB]],
    
    y0fromExB[k[p, kscale], xD, yD, rG[p, theta2[th0, rD], rD*BRxB/rRxB], phi],
    
    rG[ p, theta2[th0, rA], rA*BRxB/rRxB]
    ]
  ]
  
  
  
pminExB[{xD_?NumericQ, yD_, phi_?NumericQ, th0_, alpha_, BRxB_, rRxB_,rD_}, {xA_, yA_, xOff_, yOff_, rA_}, kscale_] = 
  Solve[
  	x0fromExB[k[p, kscale], xD, yD, rG[p, theta2[th0, rD], rD*BRxB/rRxB], phi, -D1stSimple[p, alpha, BRxB, th0, rRxB]] ==
     -xA/2 - rG[p, theta2[th0, rA], rA*BRxB/rRxB] + xOff, p][[3, 1, 2]];
    
    
pmaxExB[{xD_?NumericQ, yD_, phi_?NumericQ, th0_, alpha_, BRxB_, rRxB_,rD_}, {xA_, yA_, xOff_, yOff_, rA_}, kscale_] = 
  Solve[
  	x0fromExB[k[p, kscale], xD, yD, rG[p, theta2[th0, rD], rD*BRxB/rRxB], phi, -D1stSimple[p, alpha, BRxB, th0, rRxB]] ==
     xA/2 + rG[p, theta2[th0, rA], rA*BRxB/rRxB] + xOff, p][[3, 1, 2]];
    
pminExBCases[{xD_?NumericQ, yD_?NumericQ, phi_?NumericQ, 
   th0_?NumericQ, alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, 
   rD_?NumericQ}, {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, 
   yOff_?NumericQ, rA_?NumericQ}, kscale_?NumericQ] := Piecewise[
  {
   {0., pminExB[{xD, yD, phi, th0, alpha, BRxB, rRxB, rD}, {xA, yA, 
       xOff, yOff, rA}, kscale] < 0.},
   {pmax, 
    pminExB[{xD, yD, phi, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff, 
       yOff, rA}, kscale] > pmax}
   },
  pminExB[{xD, yD, phi, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff, 
    yOff, rA}, kscale]
]
  
pmaxExBCases[{xD_?NumericQ, yD_?NumericQ, phi_?NumericQ, 
   th0_?NumericQ, alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, 
   rD_?NumericQ}, {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, 
   yOff_?NumericQ, rA_?NumericQ}, kscale_?NumericQ] := Piecewise[
  {
   {0., pmaxExB[{xD, yD, phi, th0, alpha, BRxB, rRxB, rD}, {xA, yA, 
       xOff, yOff, rA}, kscale] < 0.},
   {pmax, 
    pmaxExB[{xD, yD, phi, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff, 
       yOff, rA}, kscale] > pmax}
   },
  pmaxExB[{xD, yD, phi, th0, alpha, BRxB, rRxB, rD}, {xA, yA, xOff, 
    yOff, rA}, kscale]
]
  
  
  IntegrationExB[
  a_?NumericQ, {alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, 
   rD_?NumericQ, rF_?NumericQ}, {xA_?NumericQ, yA_?NumericQ, 
   xOff_?NumericQ, yOff_?NumericQ, rA_?NumericQ, prec_?NumericQ}, 
  kscale_?NumericQ, BinN_?NumericQ, BinList_List,IntPrec_] :=
 ParallelTable[
  NIntegrate[
   Integrand2DExB[
    a,{BinList[[bin, 1]], BinList[[bin, 2]], phi, p, th0, alpha, 
     BRxB, rRxB, rD}, {xA, yA, xOff, yOff, rA, prec}, kscale],
   {th0, 0, thetamax[rF]},
   {phi, 0, 2*Pi},
   {
    p,
    pminExBCases[{BinList[[bin, 1]], BinList[[bin, 2]], phi, th0, 
      alpha, BRxB, rRxB, rD}, {xA, yA, xOff, yOff, rA}, kscale], 
    pmaxExBCases[{BinList[[bin, 1]], BinList[[bin, 2]], phi, th0, 
      alpha, BRxB, rRxB, rD}, {xA, yA, xOff, yOff, rA}, kscale]
    },
   PrecisionGoal -> IntPrec
   ],
  {bin, 1, BinN}, Method -> "FinestGrained"]
  
  
  IntegrationExBMC[
  a_?NumericQ, {alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, 
   rD_?NumericQ, rF_?NumericQ}, {xA_?NumericQ, yA_?NumericQ, 
   xOff_?NumericQ, yOff_?NumericQ, rA_?NumericQ, MCSteps_?NumericQ}, 
  kscale_?NumericQ, BinN_?NumericQ, BinList_List,IntPrec_] :=
 ParallelTable[
  NIntegrate[
   Integrand2DExBMC[
    a, {BinList[[bin, 1]], BinList[[bin, 2]], phi, p, th0, alpha, 
     BRxB, rRxB, rD}, {xA, yA, xOff, yOff, rA, MCSteps}, kscale],
   {th0, 0, thetamax[rF]},
   {phi, 0, 2*Pi},
   {
    p,
    pminExBCases[{BinList[[bin, 1]], BinList[[bin, 2]], phi, th0, 
      alpha, BRxB, rRxB, rD}, {xA, yA, xOff, yOff, rA}, kscale], 
    pmaxExBCases[{BinList[[bin, 1]], BinList[[bin, 2]], phi, th0, 
      alpha, BRxB, rRxB, rD}, {xA, yA, xOff, yOff, rA}, kscale]
    },
   PrecisionGoal -> IntPrec
   ],
  {bin, 1, BinN}, Method -> "FinestGrained"]
  
  
   IntegrationExBArith[
  a_?NumericQ, {alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, rD_?NumericQ, rF_?NumericQ}, {xA_?NumericQ, yA_?NumericQ, 
   xOff_?NumericQ, yOff_?NumericQ, rA_?NumericQ}, kscale_?NumericQ, BinN_?NumericQ, BinList_List,IntPrec_] :=
 ParallelTable[
  NIntegrate[
   Integrand2DExBArith[
    a,{BinList[[bin, 1]], BinList[[bin, 2]], phi, p, th0, alpha, 
     BRxB, rRxB, rD}, {xA, yA, xOff, yOff, rA}, kscale],
   {th0, 0, thetamax[rF]},
   {phi, 0, 2*Pi},
   {
    p,
    pminExBCases[{BinList[[bin, 1]], BinList[[bin, 2]], phi, th0, 
      alpha, BRxB, rRxB, rD}, {xA, yA, xOff, yOff, rA}, kscale], 
    pmaxExBCases[{BinList[[bin, 1]], BinList[[bin, 2]], phi, th0, 
      alpha, BRxB, rRxB, rD}, {xA, yA, xOff, yOff, rA}, kscale]
    },
   PrecisionGoal -> IntPrec
   ],
  {bin, 1, BinN}, Method -> "FinestGrained"]
  
  
  