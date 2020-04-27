(* Wolfram Language package *)


Get["Common/Constants.m"];
Get["Common/CommonFunctions.m"];

beta[Ee_]:=Sqrt[Ee^2-me^2]/Ee

Omega[Ee_, Ep_, lambda_]:= (1+lambda^2) * 
	( 
		E0 * Ee^2 * (1+Pi*\[Alpha]*beta[Ee]) 
		- 2/3*Ee^3 
		+ Pi*\[Alpha]*E0*me^2*Log[(1+beta[Ee])/(1-beta[Ee])]/2 
		- 2*Pi*\[Alpha]*beta[Ee]*Ee *(me^2+beta[Ee]^2*Ee^2/3)		
	) - (1-lambda^2) * mp * (epMax - Ep) * Ee * (1 + Pi*\[Alpha]*beta[Ee])
	

(*min = -1, max = 1*)
Eeminmax[Ep_,minmax_]:= 1/2* ( mn - Ep + minmax * pofTClassic[Ep-mp,mp] + me^2/(mn - Ep + minmax * pofTClassic[Ep-mp,mp]) )

wpGlueck93[]:= 