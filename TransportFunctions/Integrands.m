(* Wolfram Language package *)

Get["Spectra/ElectronSpectrum.m"];
Get["Aperture/ApertureDefinition.m"];
Get["Bfield+Drift/Bfields.m"]
Get["Bfield+Drift/DetectorGyrationWeighting.m"];
Get["Bfield+Drift/Drifts.m"];
Get["Common/CommonFunctions.m"];
Get["Aperture/ArithApert.m"];
Get["Spectra/ProtonSpectrumNachtmann.m"];


  (*::Section::*)
 (**Simplest Integrands**)


Integrand[
  b_, {xi_, p_, th0_, alpha_, BRxB_, rRxB_, rD_, rF_}, {xA_, 
   x0_}] :=
 Re[
   wmomNormedWb[\[Lambda]0,\[Kappa]0,b, p]*
    Sin[th0]*
    WeightrG[xi - x0 + D1stSimple[p, alpha, BRxB, th0, rRxB], p, 
     theta2[th0, rD], rD*BRxB/rRxB]
   ]/(xA + 2*rG[pmax, thetamax[rF], BRxB])


IntegrandPiece[
  b_, {xi_, p_, th0_, alpha_, BRxB_, rRxB_, rD_, rF_}, {xA_, x0_}] :=

  Piecewise[
  {
   {0, x0 < -xA/2 - rG[p, th0, BRxB]},
   {0, x0 > xA/2 + rG[p, th0, BRxB]}
   },
  Re[
    wmomNormedWb[\[Lambda]0,\[Kappa]0,b, p]*
     Sin[th0]*
     WeightrG[xi - x0 + D1stSimple[p, alpha, BRxB, th0, rRxB], p, 
      theta2[th0, rD], rD*BRxB/rRxB]
    ]/(xA + 2*rG[pmax, th0, BRxB])
  ]
  

IntegrandBoole[
  b_, {xi_, p_, th0_, alpha_, BRxB_, rRxB_, rD_, rF_}, {xA_, x0_}] :=

 Boole[Abs[x0] < xA/2 + rG[p, th0, BRxB]]*
  Re[
    wmomNormedWb[\[Lambda]0,\[Kappa]0,b, p]*
     Sin[th0]*
     WeightrG[xi - x0 + D1stSimple[p, alpha, BRxB, th0, rRxB], p, 
      theta2[th0, rD], rD*BRxB/rRxB]
    ]/(xA + 2*rG[pmax, th0, BRxB])
    
 
    
  
  
    (*::Section::*)
 (**Integrands with Aperture**)   
    
    
IntegrandApert[
   b_?NumericQ, {xi_?NumericQ, p_?NumericQ, th0_?NumericQ, alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, rD_?NumericQ}, {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, 
   yOff_?NumericQ, x0_?NumericQ, y0_?NumericQ, rA_?NumericQ,prec_?NumericQ}] :=
  Re[
    ApertureFunc[xA, yA, xOff, yOff, x0, y0, theta2[th0, rA], p, rA*BRxB/rRxB,prec]*
     wmomNormedWb[\[Lambda]0,\[Kappa]0,b, p]*
     Sin[th0]*
     WeightrG[xi - x0 + D1stSimple[p, alpha, BRxB, th0, rRxB], p, 
      theta2[th0, rD], rD*BRxB/rRxB]
    ]/(xA + 2*rG[pmax, th0, rA*BRxB/rRxB])
    
    
    
IntegrandApertPiecewise[
  b_, {xi_, p_, th0_, alpha_, BRxB_, rRxB_, rD_}, {xA_, yA_, xOff_, 
   yOff_, x0_, y0_, rA_,prec_}] :=
 Piecewise[
  {
   {0, x0 <= -xA/2 - rG[p, th0, BRxB]},
   {0, x0 >= xA/2 + rG[p, th0, BRxB]}
   },
  Re[
    ApertureFunc[xA, yA, xOff, yOff, x0, y0, theta2[th0, rA], p, BRxB,prec]*
     wmomNormedWb[\[Lambda]0,\[Kappa]0,b, p]*
     Sin[th0]*
     WeightrG[xi - x0 + D1stSimple[p, alpha, BRxB, th0, rRxB], p, 
      theta2[th0, rD], rD*BRxB/rRxB]
    ]/(xA + 2*rG[pmax, th0, BRxB])
  ]  
  
  
IntegrandApertBoole[
  b_, {xi_, p_, th0_, alpha_, BRxB_, rRxB_, rD_}, {xA_, yA_, xOff_, 
   yOff_, x0_, y0_, rA_,prec_}] :=
 Boole[ Abs[x0] <= xA/2 + rG[p, th0, BRxB]]*
  Re[
    ApertureFunc[xA, yA, xOff, yOff, x0, y0, theta2[th0, rA], p, BRxB,prec]*
     wmomNormedWb[\[Lambda]0,\[Kappa]0,b, p]*
     Sin[th0]*
     WeightrG[xi - x0 + D1stSimple[p, alpha, BRxB, th0, rRxB], p, 
      theta2[th0, rD], rD*BRxB/rRxB]
    ]/(xA + 2*rG[pmax, th0, rA*BRxB/rRxB])
  
  
  
  IntegrandArith[
  b_?NumericQ, {xi_?NumericQ, p_?NumericQ, th0_?NumericQ, alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, rD_?NumericQ}, {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, 
   yOff_?NumericQ, x0_?NumericQ, y0_?NumericQ, rA_?NumericQ}] :=
    
    ArithApert[xA, yA, xOff, yOff, x0, y0, rG[ p,theta2[th0, rA], rA*BRxB/rRxB]]*
     wmomNormedWb[\[Lambda]0,\[Kappa]0,b, p]*
     Sin[th0]*
     WeightrG[xi - x0 + D1stSimple[p, alpha, BRxB, th0, rRxB], p, theta2[th0, rD], rD*BRxB/rRxB]/(xA + 2*rG[pmax, theta2[th0,rA], rA*BRxB/rRxB])
    
    
 
    
ProtonIntegrandArith[
  a_?NumericQ, {xi_?NumericQ, p_?NumericQ, th0_?NumericQ, alpha_?NumericQ, BRxB_?NumericQ, rRxB_?NumericQ, rD_?NumericQ}, {xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, 
   yOff_?NumericQ, x0_?NumericQ, y0_?NumericQ, rA_?NumericQ}] :=
    
    ArithApert[xA, yA, xOff, yOff, x0, y0, rG[ p,theta2[th0, rA], rA*BRxB/rRxB]]*
     pmomNormed[p, a]*
     Sin[th0]*
     WeightrG[xi - x0 + D1stSimple[p, alpha, BRxB, th0, rRxB], p, theta2[th0, rD], rD*BRxB/rRxB]/(xA + 2*rG[pmax, theta2[th0,rA], rA*BRxB/rRxB])
     
     

    