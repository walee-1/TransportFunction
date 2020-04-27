(* Wolfram Language package *)

Get["Common/Constants.m"];
Get["Common/CommonFunctions.m"];

Gw[lambda_] := 1.2755/lambda;
\[Epsilon][lambda_] := 1 + 3 *lambda^2;

e2Min[T_] := 1/2 (del - pofTClassic[T,mn] + me^2/(del - pofTClassic[T,mn]))
e2Max[T_] := 1/2 (del + pofTClassic[T,mn] + me^2/(del + pofTClassic[T,mn]))

WpGlueck[e2_, T_, a_, b_] := 1/2*(1 + a)*e2^2*(del - 2/3*e2) + a*mn*e2*(T - tpMax) + b*me*e2 (del - 1/2*e2);

wpGlueck[T_, lambda_, a_, b_] := mn*Gw[lambda]^2*\[Epsilon][lambda]*(WpGlueck[e2Max[T], T, a, b] - WpGlueck[e2Min[T], T, a, b]);

wpNormGlueck[lambda_, a_, b_]:= wpNormGlueck[lambda, a, b] = NIntegrate[wpGlueck[T, lambda, a, b],{T, 0., tpMax},PrecisionGoal->6]

wpNormedGlueck[T_, lambda_, a_, b_] := wpGlueck[T, lambda, a, b]/wpNormGlueck[lambda, a, b]


(* ::Section:: *)
(* Momentum part *)

ppmax = Reduce[TofPClassical[p] == tpMax, p][[2, 2]]

pmomGlueck[p_, lambda_, a_, b_] := wpNormedGlueck[TofPClassical[p], lambda, a, b]*p/mp

pmomNormGlueck[lambda_, a_, b_] := 
 pmomNormGlueck[lambda, a, b] = 
  NIntegrate[pmomGlueck[p, lambda, a, b], {p, 0., ppmax}, PrecisionGoal -> 6]

pmomNormedGlueck[p_, lambda_, a_, b_] := pmomGlueck[p, lambda, a, b]/pmomNormGlueck[lambda, a, b]