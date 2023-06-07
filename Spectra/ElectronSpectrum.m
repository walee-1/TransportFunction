(* ::Package:: *)

(* Wolfram Language package *)

Get["Common/Constants.m"];
Get["Common/CommonFunctions.m"];




(* ::Section:: *)
(* electron energy *)


\[Beta][Ee_] := Sqrt[1 - me^2/Ee^2];

F[Ee_] := 
 2*Pi*\[Alpha]/\[Beta][Ee]/(1 - Exp[-2*Pi*\[Alpha]/\[Beta][Ee]]);

FPiecewise[Ee_]:=Piecewise[{{F[me+0.01],Ee<me+0.01}},F[Ee]]

Fapprox[Ee_] := 
 1 + Pi*\[Alpha]/\[Beta][Ee] + \[Alpha]^2*(11/4 - 0.5772 - 
     Log[2*\[Beta][Ee]*Ee*0.01/4/me] + Pi^2/3/\[Beta][Ee]^2);

L[z_]=Integrate[Log[Abs[1-t]]/t,{t,0,z}];

(*Get["Spectra/LInter.mx"];*)
(*LInter=Interpolation[ParallelTable[{2*\[Beta][Ee]/(1 + \[Beta][Ee]),L[2*\[Beta][Ee]/(1 + \[Beta][Ee])]},{Ee,me,E0,(E0-me)/2000}]];*)

dR[Ee_] := \[Alpha]/(2*Pi)*(
   3*Log[mp/me] - 3/4
    + 4*(ArcTanh[\[Beta][Ee]]/\[Beta][Ee] - 1)*((E0 - Ee)/(3*Ee) - 
       3/2 + Log[2*(E0 - Ee)/me])
    + 4*L[2*\[Beta][Ee]/(1 + \[Beta][Ee])]/\[Beta][Ee]
    + ArcTanh[\[Beta][Ee]]/\[Beta][
       Ee]*(2*(1 + \[Beta][Ee]^2) + (E0 - Ee)^2/(6*Ee^2) - 
       4*ArcTanh[\[Beta][Ee]])
   );

Get["Spectra/dRInter.mx"];
(*dRInter = 
 Interpolation[
  ParallelTable[{Ee, dR[Ee]}, {Ee, me + 2, 
    E0 - 1, (E0 - 1 - me - 2)/2000}]];*)

R0[\[Lambda]_,\[Kappa]_,Ee_] :=  1/(1 + 3*\[Lambda]^2)*
(2*Ee/mn + \[Lambda]^2*(10*Ee/mn - 2*me^2/(mn*Ee) - 2*E0/mn) +
 \[Lambda]*(1 + 2*\[Kappa])*(-4*Ee/mn + 2*me^2/mn/Ee + 2*E0/mn));


(*R0Inter = 
 Interpolation[
  ParallelTable[{Ee, R0[\[Lambda]0,\[Kappa]0, Ee]}, {Ee, me, 
    E0, (E0 - me)/2000}]];*)

we[\[Lambda]_,\[Kappa]_,Ee_] :=(*(4*Pi)^2/(2*Pi*hq)^6**)
 FPiecewise[Ee]*Sqrt[Ee^2 - me^2]*
  Ee*(E0 - Ee)^2*(1 + dRInter[Ee])*(1 + R0[\[Lambda],\[Kappa],Ee])/10^23;


(*weInter = 
 Interpolation[
  ParallelTable[{Ee, we[\[Lambda],\[Kappa],Ee]}, {Ee, me, E0, (E0 - me)/2000}]];*)

weNorm[\[Lambda]_,\[Kappa]_,b_] := weNorm[\[Lambda],\[Kappa],b] = NIntegrate[(1 + b*me/Ee)*we[\[Lambda],\[Kappa],Ee], {Ee, me, E0},PrecisionGoal->6];

weNormed[\[Lambda]_,\[Kappa]_,b_,Ee_]:=(1+b*me/Ee)*we[\[Lambda],\[Kappa],Ee]/weNorm[\[Lambda],\[Kappa],b];



(* ::Section:: *)
(* Electron Momentum *)


pmax = Solve[Eofp[p] == E0, p][[2, 1, 2]];

SetAttributes[pmax,{Constant,Protected}];

wmom[\[Lambda]_,\[Kappa]_,p_] := we[\[Lambda],\[Kappa],Eofp[p]]*p/Eofp[p]

wmomPiecewise[\[Lambda]_,\[Kappa]_,p_] := Piecewise[{
   {0, p <= 0},
   {0, p >= pmax}
   },
  wmom[\[Lambda],\[Kappa],p]
  ]
  
wmomNorm[\[Lambda]_,\[Kappa]_,b_] := 
 wmomNorm[\[Lambda],\[Kappa],b] = 
  NIntegrate[(1 + b*me/Eofp[p])*wmom[\[Lambda],\[Kappa],p], {p, 0, pmax}, PrecisionGoal -> 8]
  
wmomNormedWb[\[Lambda]_,\[Kappa]_,b_, p_] := (1 + b*me/Eofp[p])*
  wmomPiecewise[\[Lambda],\[Kappa],p]/wmomNorm[\[Lambda],\[Kappa],b]


(*wmomInter =  Interpolation[  ParallelTable[{p, wmomPiecewise[p]}, {p, 0, pmax, pmax/4000}]]
  
wmomInterNorm[b_] := 
 wmomInterNorm[b] = 
  NIntegrate[(1 + b*me/Eofp[p])*wmomInter[p], {p, 0, pmax}, 
   PrecisionGoal -> 3]

wmomInterNormedWb[b_, p_] := (1 + b*me/Eofp[p])*
  wmomInter[p]/wmomInterNorm[b]*)
  
  
