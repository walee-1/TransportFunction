(* Wolfram Language package *)

Get["Common/CommonFunctions.m"];

ds[x_, r_] := Sqrt[1/(r^2 - x^2)]

dsNorm = Integrate[ds[xi, r], {xi, -r, r}, Assumptions -> r > 0]

(*we introduce g here as local coordinate - 
  it's not x, its local!*)
WeightrG[g_?NumericQ, p_?NumericQ, thetaDet_?NumericQ, 
  BDet_?NumericQ] :=
 ds[g, rG[p, thetaDet, BDet]]/dsNorm
