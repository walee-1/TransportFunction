(* Wolfram Language package *)


Get["Common/Constants.m"];
Get["Common/CommonFunctions.m"];
Get["Bfield+Drift/Drifts.m"];


(*Definitions for transition circle. we need those for calculating AlphaTransition*)

(* ::Section:: *)
(* just transition region with different circle*)

CartCircle[zm_, ym_, rc_, z_] := -Sqrt[rc^2 - (z - zm)^2] - ym

rc[r_, zt_, phit_] = Quiet[Solve[
 CartCircle[zt, r - rc, rc, r*Sin[phit]] == -r*Cos[phit] &&
  {rc, zt,
     r, phit} \[Element] Reals && 0 < phit < Pi/4 && -0.2 < zt < 0. &&
   0 < r,
 rc, Reals
 ][[1, 1, 2, 1]]] // FullSimplify;

AlphaPrime[R_,AlphaTrans_,zt_]:=ArcSin[R*Sin[AlphaTrans]/rc[R,zt,AlphaTrans]]

DSumAlphaPrime[p_,th0_,alpha_,BRxB_,rRxB_,{R_,AlphaTrans_,zt_}]:= 
	2*D1stSimple[p,AlphaPrime[R,AlphaTrans,zt],BRxB,th0,rRxB] + D1stSimple[p, alpha-2*AlphaTrans, BRxB, th0, rRxB] 



(* ::Section:: *)
(* added angle in addition*)

AddedAngleScale0 = 0.00132/180*Pi;

AddedAnglePeak[scale_, p_] := scale*p/1000

DriftAddedAngleAlpha[p_, alphaPrime_, BRxB_, thRxB_, AddedScale_] := 
 p*alphaPrime/c/BRxB/Cos[thRxB + AddedAnglePeak[AddedScale, p]]*
 (3/4 + ( Sin[2*(thRxB + AddedAnglePeak[AddedScale, p])] - Sin[2*thRxB] ) / (8*AddedAnglePeak[AddedScale, p]) )
       
       (*still keep the general parameter alpha for the total curvature of the system*)
       
DSumAddedAngleAlpha[p_,th0_,alpha_,BRxB_,rRxB_,{R_,AlphaTrans_,zt_,AddedScale_}]:= 
	2*DriftAddedAngleAlpha[p,AlphaPrime[R,AlphaTrans,zt],BRxB,theta2[th0,rRxB],AddedScale] + 
							D1stSimple[p, alpha-2*AlphaTrans, BRxB, th0, rRxB]
							
							
 



