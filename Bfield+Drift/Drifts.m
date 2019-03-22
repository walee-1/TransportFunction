(* Wolfram Language package *)

Get["Bfields.m"];
Get["Constants.m"];


fReduced[th0_, 
  rRxB_] := ((2 - Sin[th0]^2*rRxB)/Sqrt[1 - Sin[th0]^2*rRxB])/2

D1stSimple[p_, alpha_, B_, th0_, rRxB_] := 
 p*alpha/c/B*fReduced[th0, rRxB]
 
 
 D1stBPoly[p_, alpha_, BRxB_, th0_, rRxB_, y0_, R0_] := 
 p*alpha/c*fReduced[th0, rRxB]/Bpolynom[y0, BRxB, R0]
D1stBPolyGrad[p_, alpha_, th0_, BRxB_, rRxB_, y0_, R0_, G1_, G2_] := 
 p*alpha/c*fReduced[th0, rRxB]/BpolynomGrad[y0, BRxB, R0, G1, G2]
 
 
 DriftPolywrGWg[p_?NumericQ, alpha_?NumericQ, th0_?NumericQ, 
  rRxB_?NumericQ, BRxB_?NumericQ, R_?NumericQ, y0_?NumericQ, 
  G1_?NumericQ, G2_?NumericQ] := 
 p*alpha*fReduced[th0, rRxB]/c/
   BBarPolyWg[p, theta2[th0, rRxB], BRxB, R, G1, G2, y0]