(* Wolfram Language package *)

(* ::Section:: *)
(* Energy *)


sigma[T_] := 1 - 2*T*mn/del^2;

g1[T_] := (1 - xr^2/sigma[T])^2*
   Sqrt[1 - 
     sigma[T]]*(4*(1 + 
        xr^2/sigma[T]) - (4/3*(sigma[T] - xr^2)/sigma[T])*(1 - 
        sigma[T]));

g2[T_] := (1 - xr^2/sigma[T])^2*
   Sqrt[1 - sigma[T]]*(4*(1 + xr^2/sigma[T] - 2*sigma[T]) - 
     4/3*(sigma[T] - xr^2)/sigma[T]*(1 - sigma[T]));

g1Inter = 
  Interpolation[Table[{T, g1[T]}, {T, 0, tpMax, tpMax/4000}]];
g2Inter = Interpolation[Table[{T, g2[T]}, {T, 0, tpMax, tpMax/4000}]];

dwdt[T_, a_] := g1[T] + a*g2[T];

dwdtInter[T_, a_] := g1Inter[T] + a*g2Inter[T];

dwdtNorm[a_] := 
 dwdtNorm[a] = 
  NIntegrate[dwdt[T, a], {T, 0, tpMax}, PrecisionGoal -> 3]
  
dwdtNormed[T_, a_] := dwdt[T, a]/dwdtNorm[a]

(* ::Section:: *)
(* momentum *)


ppmax = Reduce[TofPClassical[p] == tpMax, p][[2, 2]]

pmom[p_, a_] := dwdtNormed[TofPClassical[p], a]*p/mp

pmomNorm[a_] := 
 pmomNorm[a] = 
  NIntegrate[pmom[p, a], {p, 0., ppmax}, PrecisionGoal -> 3]

pmomNormed[p_, a_] := pmom[p, a]/pmomNorm[a]




