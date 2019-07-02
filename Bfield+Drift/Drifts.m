(* Wolfram Language package *)

Get["Bfield+Drift/Bfields.m"];
Get["Common/Constants.m"];
Get["Common/CommonFunctions.m"];

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
   
 
pminFunc[x0_, x_, thD_, BD_, alpha_, th0_, rRxB_, 
  BRxB_] := (x0 - x)/(Sin[thD]/c/BD + alpha*fReduced[th0, rRxB]/c/BRxB)
pmaxFunc[x0_, x_, thD_, BD_, alpha_, th0_, rRxB_, 
  BRxB_] := (x - x0)/(Sin[thD]/c/BD - 
    alpha*fReduced[th0, rRxB]/c/BRxB)
    
pminCases[x0_, x_, thD_, BD_, alpha_, th0_, rRxB_, BRxB_] := Which[
    pminFunc[x0, x, thD, BD, alpha, th0, rRxB, BRxB] > pmax, pmax,
    pminFunc[x0, x, thD, BD, alpha, th0, rRxB, BRxB] < 0, 0,
    True, pminFunc[x0, x, thD, BD, alpha, th0, rRxB, BRxB]
    ]
pmaxCases[x0_, x_, thD_, BD_, alpha_, th0_, rRxB_, BRxB_] := Which[
    pmaxFunc[x0, x, thD, BD, alpha, th0, rRxB, BRxB] > pmax, pmax,
    pmaxFunc[x0, x, thD, BD, alpha, th0, rRxB, BRxB] < 0, 0,
    True, pmaxFunc[x0, x, thD, BD, alpha, th0, rRxB, BRxB]
    ]
    