(* Wolfram Language package *)

Get["Constants.m"];
Get["CommonFunctions.m"];




(* ::Section:: *)
(* electron energy *)

\[Beta][Ee_] := Sqrt[1 - me^2/Ee^2];

F[Ee_] := 
 2*Pi*\[Alpha]/\[Beta][Ee]/(1 - Exp[-2*Pi*\[Alpha]/\[Beta][Ee]]);

Fapprox[Ee_] := 
 1 + Pi*\[Alpha]/\[Beta][Ee] + \[Alpha]^2*(11/4 - 0.5772 - 
     Log[2*\[Beta][Ee]*Ee*0.01/4/me] + Pi^2/3/\[Beta][Ee]^2);

L[z_]:=Integrate[Log[Abs[1-t]]/t,{t,0,z}];

dR[Ee_] := \[Alpha]/(2*Pi)*(
   3*Log[mp/me] - 3/4
    + 4*(ArcTanh[\[Beta][Ee]]/\[Beta][Ee] - 1)*((E0 - Ee)/(3*Ee) - 
       3/2 + Log[2*(E0 - Ee)/me])
    + 4*L[2*\[Beta][Ee]/(1 + \[Beta][Ee])]/\[Beta][Ee]
    + ArcTanh[\[Beta][Ee]]/\[Beta][
       Ee]*(2*(1 + \[Beta][Ee]^2) + (E0 - Ee)^2/(6*Ee^2) - 
       4*ArcTanh[\[Beta][Ee]])
   );

dRInter = 
 Interpolation[
  ParallelTable[{Ee, dR[Ee]}, {Ee, me + 2, 
    E0 - 1, (E0 - 1 - me - 2)/1000}]];

R0[\[Lambda]_, Ee_] :=  1/(1 + 3*\[Lambda]^2)*(2*Ee/mn + \[Lambda]^2*(10*Ee/mn - 2*me^2/(mn*Ee) - 2*E0/mn) + \[Lambda]*(1 + 2*1.85)*(-4*Ee/mn + 2*me^2/mn/Ee + 2*E0/mn));

\[Lambda]0 = -1.2723;

R0Inter = 
 Interpolation[
  ParallelTable[{Ee, R0[\[Lambda]0, Ee]}, {Ee, me, 
    E0, (E0 - me)/1000}]];

we[Ee_] :=(*(4*Pi)^2/(2*Pi*hq)^6**)
 F[Ee]*Sqrt[Ee^2 - me^2]*
  Ee*(E0 - Ee)^2*(1 + dRInter[Ee])*(1 + R0Inter[Ee])/10^23;


weInter = 
 Interpolation[
  ParallelTable[{Ee, we[Ee]}, {Ee, me, E0, (E0 - me)/2000}]];

weNorm[b_] := weNorm[b] = NIntegrate[(1 + b*me/Ee)*weInter[Ee], {Ee, me, E0}];

weNormed[b_,Ee_]:=(1+b*me/Ee)*weInter[Ee]/weNorm[b];

(* ::Section:: *)
(* Electron Momentum *)

pmax = Solve[Eofp[p] == E0, p][[2, 1, 2]];

wmom[p_] := we[Eofp[p]]*p*c/Eofp[p]

wmomPiecewise[p_] := Piecewise[{
   {0, p <= 0},
   {0, p >= pmax}
   },
  wmom[p]
  ]

wmomInter =  Interpolation[  ParallelTable[{p, wmomPiecewise[p]}, {p, 0, pmax, pmax/4000}]]
  
wmomInterNorm[b_] := 
 wmomInterNorm[b] = 
  NIntegrate[(1 + b*me/Eofp[p])*wmomInter[p], {p, 0, pmax}, 
   PrecisionGoal -> 3]

wmomInterNormedWb[b_, p_] := (1 + b*me/Eofp[p])*
  wmomInter[p]/wmomInterNorm[b]
  
  
