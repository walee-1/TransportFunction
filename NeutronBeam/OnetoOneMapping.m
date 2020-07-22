(* Wolfram Language package *)

(*Get["NeutronBeam/ManualShift_NBeamInt.m"];*)


yDofDV[phiDV_, yDV_, p_, th0_, BRxB_, rRxB_, rD_, phiDet_, 
  YDetShift_] = 
 Solve[yDV == 
    DeltayDVShift[phiDV, yD, p, th0, BRxB, rRxB, rD, phiDet, 
     YDetShift], yD][[1, 1, 2]];
     
xDofDV[phiDV_, yDV_, xDV_, p_, th0_, alpha_, BRxB_, rRxB_, rD_, 
  phiDet_, R_, G1_, G2_, {YRxBShift_, XDetShift_, YDetShift_}] = 
 Solve[xDV == 
    DeltaxDV2Shift[phiDV, 
     yDofDV[phiDV, yDV, p, th0, BRxB, rRxB, rD, phiDet, YDetShift], 
     xD, p, th0, alpha, BRxB, rRxB, rD, phiDet, R, G1, 
     G2, {YRxBShift, XDetShift, YDetShift}], xD][[1, 1, 2]];     
     
     