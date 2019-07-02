(* Wolfram Language package *)

Get["Common/CommonFunctions.m"];

Bpolynom[y0_, B0_, R0_] := 3*B0 - 3*B0*(R0+y0)/R0 + B0*(R0+y0)^2/R0^2
BpolynomGrad[y0_, BRxB_, R0_, G1_, G2_] := 
 BRxB - G1*BRxB/R0*(y0) + G2*BRxB/R0^2*(y0)^2

B1overr[y0_, B0_, R_] := B0*R/y0
B1overrWgW\[CapitalDelta]G[B0_, R0_, y0_, \[CapitalDelta]G1_] := 
 B0*R0/y0 + \[CapitalDelta]G1*B0/R0*(y0 - R0)
 
 
 (* ::Section:: *)
 (* Different Gyration effects in RxB *)
 
 (* ::Subsection:: *)
 (* polynom B without r(g) *)

CirclexBPoly[r_, g_] := Sqrt[r^2 - g^2]

BPolygLimits[p_, th2_, B0_, R0_, G1_, G2_, y0_] = 
 Solve[CirclexBPoly[r, g] == 0, g][[All, 1, 2]];
 
 dxdg[r_, g_] = D[CirclexBPoly[r, g], g];

 dsdxNormed[r_, g_] = Sqrt[1 + dxdg[r, g]^2]/Pi/r;
 
 BBarPoly[y0_, BRxB_, R_, G1_, G2_, r_] = 
 Integrate[
   dsdxNormed[r, g]*BpolynomGrad[y0 + g, BRxB, R, G1, G2], {g, -r, 
    r}][[1]];
    
    (* ::Subsection:: *)
    (* polynom with r(g) *)
    
    CirclexBPolyWg[p_, th2_, B0_, R0_, G1_, G2_, y0_, g_] := 
 Sqrt[rG[p, th2, BpolynomGrad[y0 + g, B0, R0, G1, G2]]^2 - g^2]

BPolygLimitsWg[p_, th2_, B0_, R_, G1_, G2_, y0_] = 
  Solve[CirclexBPolyWg[p, th2, B0, R, G1, G2, y0, g] == 0, g][[{1, 4},
    1, 2]];
    
 dxdgBPolyWg[p_, th2_, B_, R_, G1_, G2_, y0_, g_] = 
 D[CirclexBPolyWg[p, th2, B, R, G1, G2, y0, g], g];   
    
  dsdxNormBPolyWg[p_, th2_, B_, R_, G1_, G2_, y0_] :=
 
 dsdxNormBPolyWg[p, th2, B, R, G1, G2, y0] = 
  Re[NIntegrate[
    Sqrt[1 + dxdgBPolyWg[p, th2, B, R, G1, G2, y0, g]^2], {g, 
     BPolygLimitsWg[p, th2, B, R, G1, G2, y0][[1]], 
     BPolygLimitsWg[p, th2, B, R, G1, G2, y0][[2]]}, 
    PrecisionGoal -> 4]]
    
 dsdxNormedBPolyWg[p_?NumericQ, th2_?NumericQ, B_?NumericQ, 
  R_?NumericQ, G1_?NumericQ, G2_?NumericQ, y0_?NumericQ, 
  g_?NumericQ] := 
 Sqrt[1 + dxdgBPolyWg[p, th2, B, R, G1, G2, y0, g]^2]/
  dsdxNormBPolyWg[p, th2, B, R, G1, G2, y0]
  
 BBarPolyWg[p_?NumericQ, th2_?NumericQ, B_?NumericQ, R_?NumericQ, 
  G1_?NumericQ, G2_?NumericQ, y0_?NumericQ] := Re[Quiet[NIntegrate[
    BpolynomGrad[y0 + g, B, R, G1, G2]*
     dsdxNormedBPolyWg[p, th2, B, R, G1, G2, y0, g],
    {g, BPolygLimitsWg[p, th2, B, R, G1, G2, y0][[1]], 
     BPolygLimitsWg[p, th2, B, R, G1, G2, y0][[2]]}, 
    PrecisionGoal -> 4]]]
    
    BBarPolyWg::usage="Gyration Mean B-field of polynomial form  with arguments: p_?NumericQ, th2_?NumericQ, B_?NumericQ, R_?NumericQ, 
  G1_?NumericQ, G2_?NumericQ, y0_?NumericQ"; 
  
  
     
    
    
    
    
      
    
    
    