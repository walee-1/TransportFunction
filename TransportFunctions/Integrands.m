(* Wolfram Language package *)

Get["Spectra/ElectronSpectrum.m"];
Get["Aperture/ApertureDefinition.m"];
Get["Bfield+Drift/Bfields.m"]
Get["Bfield+Drift/DetectorGyrationWeighting.m"];
Get["Bfield+Drift/Drifts.m"];
Get["CommonFunctions.m"];

Integrand[
  b_, {xi_, p_, th0_, alpha_, BRxB_, rRxB_, rD_, rF_}, {xA_, 
   x0_}] :=
 Re[
   wmomInterNormedWb[b, p]*
    Sin[th0]*
    WeightrG[xi - x0 + D1stSimple[p, alpha, BRxB, th0, rRxB], p, 
     theta2[th0, rD], rD*BRxB/rRxB]
   ]/(xA + 2*rG[pmax, thetamax[rF], BRxB])
(*]*)


IntegrandPiece[
  b_, {xi_, p_, th0_, alpha_, BRxB_, rRxB_, rD_, rF_}, {xA_, x0_}] :=

  Piecewise[
  {
   {0, x0 < -xA/2 - rG[p, th0, BRxB]},
   {0, x0 > xA/2 + rG[p, th0, BRxB]}
   },
  Re[
    wmomInterNormedWb[b, p]*
     Sin[th0]*
     WeightrG[xi - x0 + D1stSimple[p, alpha, BRxB, th0, rRxB], p, 
      theta2[th0, rD], rD*BRxB/rRxB]
    ]/(xA + 2*rG[pmax, th0, BRxB])
  ]
  
  
  
IntegrandBoole[
  b_, {xi_, p_, th0_, alpha_, BRxB_, rRxB_, rD_, rF_}, {xA_, x0_}] :=

 Boole[Abs[x0] < xA/2 + rG[p, th0, BRxB]]*
  Re[
    wmomInterNormedWb[b, p]*
     Sin[th0]*
     WeightrG[xi - x0 + D1stSimple[p, alpha, BRxB, th0, rRxB], p, 
      theta2[th0, rD], rD*BRxB/rRxB]
    ]/(xA + 2*rG[pmax, th0, BRxB])
    
    
IntegrandApert[
  b_, {xi_, p_, th0_, alpha_, BRxB_, rRxB_, rD_}, {xA_, yA_, xOff_, 
   yOff_, x0_, y0_, rA_,prec_}] :=
  Re[
    ApertureFunc[xA, yA, xOff, yOff, x0, y0, theta2[th0, rA], p, BRxB,prec]*
     wmomInterNormedWb[b, p]*
     Sin[th0]*
     WeightrG[xi - x0 + D1stSimple[p, alpha, BRxB, th0, rRxB], p, 
      theta2[th0, rD], rD*BRxB/rRxB]
    ]/(xA + 2*rG[pmax, th0, BRxB])
    
    
IntegrandApertPiecewise[
  b_, {xi_, p_, th0_, alpha_, BRxB_, rRxB_, rD_}, {xA_, yA_, xOff_, 
   yOff_, x0_, y0_, rA_,prec_}] :=
 Piecewise[
  {
   {0, x0 < -xA/2 - rG[p, th0, BRxB]},
   {0, x0 > xA/2 + rG[p, th0, BRxB]}
   },
  Re[
    ApertureFunc[xA, yA, xOff, yOff, x0, y0, theta2[th0, rA], p, BRxB,prec]*
     wmomInterNormedWb[b, p]*
     Sin[th0]*
     WeightrG[xi - x0 + D1stSimple[p, alpha, BRxB, th0, rRxB], p, 
      theta2[th0, rD], rD*BRxB/rRxB]
    ]/(xA + 2*rG[pmax, th0, BRxB])
  ]  
  
  
IntegrandApertBoole[
  b_, {xi_, p_, th0_, alpha_, BRxB_, rRxB_, rD_}, {xA_, yA_, xOff_, 
   yOff_, x0_, y0_, rA_,prec_}] :=
 Boole[ Abs[x0] < xA/2 + rG[p, th0, BRxB]]*
  Re[
    ApertureFunc[xA, yA, xOff, yOff, x0, y0, theta2[th0, rA], p, BRxB,prec]*
     wmomInterNormedWb[b, p]*
     Sin[th0]*
     WeightrG[xi - x0 + D1stSimple[p, alpha, BRxB, th0, rRxB], p, 
      theta2[th0, rD], rD*BRxB/rRxB]
    ]/(xA + 2*rG[pmax, th0, BRxB])
  
    
    
    