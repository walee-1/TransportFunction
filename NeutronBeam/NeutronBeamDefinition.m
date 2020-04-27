
(*to show how I come to the description*)
FullSimplify[
 Solve[{0 == k1*x1 + d1 && y2 == k1*x2 + d1 && pl == x3 - x2 && 
    tw == pl + x4 - x3 + x2 - x1 && y2 == k2*x2 + d2 && 
    y3 == k2*x3 + d2 && y3 == k3*x3 + d3 && 0 == k3*x4 + d3 && 
    x1 == -tw/2 && x4 == tw/2 && k1 > 0 && 
    k2 < 0 && {k1, x1, d1, y2, k1, x2, pl, x3, tw, x4, d2, y3, 
      k3} \[Element] Reals && pl < tw}, {x1, x2, x3, x4, y2, y3, d1, 
   d2, d3}]];
   
   
   TrapezNBeamNew[xn_, {tw_, pl_, k1_, k2_, k3_}] := Piecewise[
  {
   {0, xn < -tw/2},
   {k1*xn + k1*tw/2, -tw/2 <= 
     xn < -((2 k2 pl - 2 k3 pl + (k1 + k3) tw)/(2 (k1 - k3)))},
   {k2*xn + (-2 (k1 - k2) (k2 - k3) pl + (k1 (k2 - 2 k3) + 
         k2 k3) tw)/(
     2 (k1 - k3)), -((2 k2 pl - 2 k3 pl + (k1 + k3) tw)/(
      2 (k1 - k3))) <= 
     xn < -((2 (-k1 + k2) pl + (k1 + k3) tw)/(2 (k1 - k3)))},
   {k3*xn - (k3 tw)/
     2, -((2 (-k1 + k2) pl + (k1 + k3) tw)/(2 (k1 - k3))) <= xn <= 
     tw/2},
   {0, xn > tw/2}
   }
  ]
  
  TrapezNBeamNewNorm[{tw_, pl_, k1_, k2_, k3_}] := 
 TrapezNBeamNewNorm[{tw, pl, k1, k2, k3}] = 
  NIntegrate[
   TrapezNBeamNew[xn, {tw, pl, k1, k2, k3}], {xn, -tw/2, tw/2},PrecisionGoal->7,AccuracyGoal->8]
   
   
   TrapezNBeamNewNormed[xn_, {tw_, pl_, k1_, k2_, k3_}] := 
 TrapezNBeamNew[xn, {tw, pl, k1, k2, k3}]/
  TrapezNBeamNewNorm[{tw, pl, k1, k2, k3}]



(* ::Section:: *)
(* Compiled version *)

TrapezNBeamCompiled = Compile[
  {{xn, _Real}, {tw, _Real}, {pl, _Real}, {k1, _Real}, {k2, _Real}, {k3, _Real}},
  Piecewise[
   {
    {0, xn < -tw/2},
    {k1*xn + k1*tw/2, -tw/2 <= 
      xn < -((2 k2 pl - 2 k3 pl + (k1 + k3) tw)/(2 (k1 - k3)))}, {k2*
       xn + (-2 (k1 - k2) (k2 - k3) pl + (k1 (k2 - 2 k3) + 
            k2 k3) tw)/(2 (k1 - k3)), -((2 k2 pl - 
           2 k3 pl + (k1 + k3) tw)/(2 (k1 - k3))) <= 
      xn < -((2 (-k1 + k2) pl + (k1 + k3) tw)/(2 (k1 - k3)))},
    {k3*xn - (k3 tw)/
       2, -((2 (-k1 + k2) pl + (k1 + k3) tw)/(2 (k1 - k3))) <= xn <= 
      tw/2},
    {0, xn > tw/2}
    }], CompilationTarget -> "C", RuntimeAttributes -> {Listable}(*, RuntimeOptions -> "Speed"*)
  ];

TrapezNBeamCompiledNorm[{tw_, pl_, k1_, k2_, k3_}] := 
 TrapezNBeamCompiledNorm[{tw, pl, k1, k2, k3}] = 
  NIntegrate[TrapezNBeamNew[xn, {tw, pl, k1, k2, k3}], {xn, -tw/2, tw/2},PrecisionGoal->7,AccuracyGoal->8]
  
   
TrapezNBeamCompiledNormed[xn_, {tw_, pl_, k1_, k2_, k3_}] := 
 TrapezNBeamCompiled[xn, tw, pl, k1, k2, k3]/
  TrapezNBeamCompiledNorm[{tw, pl, k1, k2, k3}]

