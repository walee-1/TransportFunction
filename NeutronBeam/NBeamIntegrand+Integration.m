(* Wolfram Language package *)
Get["Spectra/ElectronSpectrum.m"];
Get["Bfield+Drift/Bfields.m"]
Get["Bfield+Drift/Drifts.m"];
Get["Common/CommonFunctions.m"];
Get["NeutronBeam/NeutronBeamDefinition.m"];


ApertBoole[xtraj_, ytraj_, {xA_, yA_, xOff_, yOff_}] := 
 Boole[xtraj < xOff + xA/2 && xtraj > xOff - xA/2 && 
   ytraj < yOff + yA/2 && ytraj > yOff - yA/2]
   
xGCA[xDV_, phiDV_, p_, th0_, BRxB_, rRxB_, 
  rA_] := (xDV - rG[p, th0, BRxB/rRxB]*Cos[phiDV])*Sqrt[1/rA]
yGCA[yDV_, phiDV_, p_, th0_, BRxB_, rRxB_, 
  rA_] := (yDV - rG[p, th0, BRxB/rRxB]*Sin[phiDV])*Sqrt[1/rA]
  
  
DeltayDV[phiDV_, yD_, p_, th0_, BRxB_, rRxB_, rA_, rD_, 
  phiDet_] := (yD - rG[p, th0, rD*BRxB/rRxB]*Sin[phiDet])*
   Sqrt[rD/rRxB]*Sqrt[rA] + rG[p, th0, BRxB/rRxB]*Sin[phiDV]
DeltaxDV[phiDV_, yDV_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, 
  phiDet_, R_, G1_, 
  G2_] := ((xD - rG[p, th0, rD*BRxB/rRxB]*Cos[phiDet])*Sqrt[rD/rRxB] -
      D1stBPolyGrad[p, alpha, th0, BRxB, rRxB, 
      yGCA[yDV, phiDV, p, th0, BRxB, rRxB, rA], R, G1, G2])*Sqrt[rA] +
   rG[p, th0, BRxB/rRxB]*Cos[phiDV]
   
   
DeltaxDV2[phiDV_, yD_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, 
  phiDet_, R_, G1_, G2_] := 
 DeltaxDV[phiDV, 
  DeltayDV[phiDV, yD, p, th0, BRxB, rRxB, rA, rD, phiDet], xD, p, th0,
   alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2]
   
   
xA[phiA_, yD_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_,
   R_, G1_, G2_] := 
 xGCA[DeltaxDV2[phiDV, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, 
    phiDet, R, G1, G2], phiDV, p, th0, BRxB, rRxB, rA] + 
  rG[p, th0, rA*BRxB/rRxB]*Cos[phiA]
  
yA[phiA_, yD_, p_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_] := 
 yGCA[DeltayDV[phiDV, yD, p, th0, BRxB, rRxB, rA, rD, phiDet], phiDV, 
   p, th0, BRxB, rRxB, rA] + rG[p, th0, rA*BRxB/rRxB]*Sin[phiA]
   
   
   
   
Integrand2DwNBeamv31[b_, yD_, 
  xD_, {phiDV_?NumericQ, phiA_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, {alpha_, BRxB_, rRxB_, rA_,
    rD_, R_, G1_, G2_}, {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, 
   k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] :=
 
 wmomNormedWb[\[Lambda]0, \[Kappa]0, b, p]*Sin[th0]*
  
  TrapezNBeamNewNormed[
   DeltaxDV2[phiDV, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,
     R, G1, G2]-xOff, {twx, plx, k1x, k2x, k3x}]*
  TrapezNBeamNewNormed[
   DeltayDV[phiDV, yD, p, th0, BRxB, rRxB, rA, rD, phiDet]-yOff, {twy, ply,
     k1y, k2y, k3y}]*
  
  ApertBoole[
   xA[phiA, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, 
    G2], yA[phiA, yD, p, th0, BRxB, rRxB, rA, rD, phiDet], {xAA, yAA, 
    xOff, yOff}]   
   
   
Integration2DwNBeamc31[b_, XYData_List, {alpha_, BRxB_,rF_, rA_,rRxB_,rD_, R_, G1_, G2_},
	 {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},IntPrec_]:=

ParallelTable[
  NIntegrate[
   Integrand2DwNBeamv31[b, XYData[[bin,2]], XYData[[bin,1]],
   	 {phiDV, phiA, phiDet, p, th0},
   	 {alpha,BRxB, rRxB, rA, rD, R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA,yAA, xOff, yOff}
   	 ],
     
     {th0, 0, thetamax[rF]}, {phiDet, 0, 2*Pi}, {phiA, 0, 2 Pi}, {phiDV, 0, 2*Pi},{p, 0, pmax},
   PrecisionGoal -> IntPrec, Method -> "AdaptiveMonteCarlo", AccuracyGoal -> 5, MaxPoints -> 2500000]
	, {bin,1,Length[XYData]}, Method -> "FinestGrained"]
  
  
  
  
  
     
   