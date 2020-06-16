(* Wolfram Language package *)


Get["Common/Constants.m"];
Get["Common/CommonFunctions.m"];


beta[Ee_]:=Sqrt[Ee^2-me^2]/Ee


Fapprox[Ee_] := 
 1 + Pi*\[Alpha]/\[Beta][Ee] + \[Alpha]^2*(11/4 - 0.5772 - 
     Log[2*\[Beta][Ee]*Ee*0.01/4/me] + Pi^2/3/\[Beta][Ee]^2);

Fapproxapprox[Ee_]:=1+ \[Alpha]*Pi/beta[Ee]



Enumax = del - (del^2 + me^2)/(2*mn);

Enu[Ee_,Ep_]:=mn - Ee - Ep

DV[Ee_,Ep_]:=Ee*(E0-Ee) + Enu[Ee,Ep]*(Enumax-Enu[Ee,Ep]) - mp*(epMax-Ep)
DA[Ee_,Ep_]:=Ee*(E0-Ee) + Enu[Ee,Ep]*(Enumax-Enu[Ee,Ep]) + mp*(epMax-Ep)
DI[Ee_,Ep_]:=2*(Ee*(E0-Ee)-Enu[Ee,Ep]*(Enumax-Enu[Ee,Ep]))


W0[Ee_,Ep_,lambda_,kappa_]:=mn* 	(
			DV[Ee,Ep]
			+ lambda^2 *DA[Ee,Ep]
			+ lambda * (1+2*kappa) * DI[Ee,Ep]
			)
			
(*min = -1, max = 1*)
Eeminmax[Ep_,minmax_]:= 1/2* ( mn - Ep + minmax * pofTClassic[Ep-mp,mp] + me^2/(mn - Ep + minmax * pofTClassic[Ep-mp,mp]) )
			
omega0C[Ep_?NumericQ,lambda_,kappa_]:=NIntegrate[
				W0[Ee,Ep,lambda,kappa] * Fapproxapprox[Ee],
				{Ee,Eeminmax[Ep,-1],Eeminmax[Ep,1]}
]

(*rC correction*)
rCTable={{0.1,0.07},{0.2,0.08},{0.3,0.09},{0.4,0.11},{0.42,0.12},{0.43,0.14},{0.45,0.11},{0.5,0.09},{0.6,0.08},{0.8,0.06},{0.9,0.06}};
rCInter[y_]=Interpolation[rCTable,InterpolationOrder->1,"ExtrapolationHandler"->"WarningMessage"->False][y];
y[Ep_]:=(Ep-mp)/(epMax-mp)
rCInterE[Ep_]:=rCInter[y[Ep]]


(*rC correction*)
rpTable={
	{0.1,0.12},{0.2,0.11},{0.3,0.1},{0.4,0.08},{0.5,0.05},{0.55,0.04},{0.6,0.01},{0.65,-0.02},{0.7,-0.06},{0.75,-0.12},
	{0.78,-0.16},{0.8,-0.2},{0.83,-0.26},{0.85,-0.32},{0.88,-0.43},{0.9,-0.52},{0.92,-0.63},{0.94,-0.79},{0.96,-1.},{0.98,-1.34}};
rpInter[y_]=Interpolation[rpTable,InterpolationOrder->2,"ExtrapolationHandler"->"WarningMessage"->False][y];
rpInterE[Ep_]:=rpInter[y[Ep]]

rroh=1.505;



omega0Calpha[Ep_,lambda_,kappa_]:=omega0C[Ep,lambda,kappa]*(1+0.01*rCInterE[Ep])*(1+0.01*rroh + 0.01*rpInterE[Ep])


(*below formula is the approximation of the electron energy integration in the paper, equ. 3.12*)
(*Omega[Ee_, Ep_, lambda_]:= (1+lambda^2) * 
	( 
		E0 * Ee^2 * (1+Pi*\[Alpha]*beta[Ee]) 
		- 2/3*Ee^3 
		+ Pi*\[Alpha]*E0*me^2*Log[(1+beta[Ee])/(1-beta[Ee])]/2 
		- 2*Pi*\[Alpha]*beta[Ee]*Ee *(me^2+beta[Ee]^2*Ee^2/3)		
	) - (1-lambda^2) * mp * (epMax - Ep) * Ee * (1 + Pi*\[Alpha]*beta[Ee])
	*)
	

wpGlueck93Norm[lambda_,kappa_]:=wpGlueck93Norm[lambda,kappa]= NIntegrate[omega0Calpha[Ep,lambda,kappa],{Ep,mp,epMax},PrecisionGoal->6]


wpNormedGlueck93[Tp_, lambda_, kappa_] := omega0Calpha[Tp + mp, lambda, kappa]/wpGlueck93Norm[lambda, kappa]




ppmax = Reduce[TofPClassical[p] == tpMax, p][[2, 2]]

pmomGlueck93[p_, lambda_,kappa_] := omega0Calpha[mp+TofPClassical[p], lambda, kappa]*p/mp

pmomNormGlueck93[lambda_, kappa_] := 
 pmomNormGlueck93[lambda, kappa] = 
  NIntegrate[pmomGlueck93[p, lambda, kappa], {p, 0., ppmax}, PrecisionGoal -> 6]

pmomNormedGlueck93[p_, lambda_, kappa_] := pmomGlueck93[p, lambda, kappa]/pmomNormGlueck93[lambda, kappa]

 