(* Wolfram Language package *)

Get["Spectra/GlueckProtonSpectrum.m"];
Get["Bfield+Drift/Drifts.m"];
Get["Common/CommonFunctions.m"];
Get["NeutronBeam/NeutronBeamDefinition.m"];
Get["Aperture/BooleApert.m"];


(*Acceleration formulas*)
pacc[p_,U_]:=Sqrt[2*mp*(p^2/(2*mp)-U)]

thDetAcc[p_,th0_,rD_,U_]:=ArcSin[Sqrt[ rD*p^2*Sin[th0]^2/(2*mp*(p^2/(2*mp)-U))]]

(*ExB drift formulas*)
kExB[p_, scale_] := scale*(20 - 6.5*p/ppmax)/50

(* ::Section:: *)
(*forward transfer definitions*)
xGCAF[xDV_, phiDV_, p_, th0_, BRxB_, rRxB_, rA_] := 
	(xDV - rG[p, th0, BRxB/rRxB]*Cos[phiDV])*Sqrt[1/rA]
yGCAF[yDV_, phiDV_, p_, th0_, BRxB_, rRxB_, rA_] := 
	(yDV - rG[p, th0, BRxB/rRxB]*Sin[phiDV])*Sqrt[1/rA]
  

xRxBGCF[phiDV_,xDV_, yDV_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, R_, G1_, G2_]:=
  xGCAF[xDV,phiDV,p,th0,BRxB,rRxB,rA] + D1stBPolyGrad[p, alpha, th0, BRxB, rRxB, yGCAF[yDV, phiDV, p, th0, BRxB, rRxB, rA], R, G1, G2]
  
(*ExB drift applied as well as Det Spread*)
xDetGCF[phiDV_,xDV_, yDV_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, R_, G1_, G2_, rD_, ExBscale_]:=
	(
	xRxBGCF[phiDV,xDV, yDV, p, th0, alpha, BRxB, rRxB, rA, R, G1, G2] 
	- kExB[p, ExBscale]*yGCAF[yDV, phiDV, p, th0, BRxB, rRxB, rA]
	)*Sqrt[rRxB/rD]
	
yDetGCF[phiDV_,xDV_, yDV_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, R_, G1_, G2_, rD_, ExBscale_]:=
	(
	yGCAF[yDV, phiDV, p, th0, BRxB, rRxB, rA] 
	+ kExB[p, ExBscale]*xRxBGCF[phiDV,xDV, yDV, p, th0, alpha, BRxB, rRxB, rA, R, G1, G2]
	)*Sqrt[rRxB/rD]
	
xDetF[phiDV_,xDV_, yDV_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, R_, G1_, G2_, rD_, ExBscale_, U_, phiDet_]:=
	xDetGCF[phiDV,xDV, yDV, p, th0, alpha, BRxB, rRxB, rA, R, G1, G2, rD, ExBscale] + rG[ pacc[p,U], thDetAcc[p,th0,rD,U], rD*BRxB/rRxB] * Cos[phiDet]
	
yDetF[phiDV_,xDV_, yDV_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, R_, G1_, G2_, rD_, ExBscale_, U_, phiDet_]:=
	yDetGCF[phiDV,xDV, yDV, p, th0, alpha, BRxB, rRxB, rA, R, G1, G2, rD, ExBscale] + rG[ pacc[p,U], thDetAcc[p,th0,rD,U], rD*BRxB/rRxB] * Sin[phiDet]
	
(* ::Section:: *)	
(* backwards transfer definitions *)
xDetGCB[xDet_, p_, th0_, phiDet_, BRxB_, rRxB_, rD_, U_]:= xDet - rG[ pacc[p,U], thDetAcc[p,th0,rD,U], rD*BRxB/rRxB ] * Cos[phiDet]
yDetGCB[yDet_, p_, th0_, phiDet_, BRxB_, rRxB_, rD_, U_]:= yDet - rG[ pacc[p,U], thDetAcc[p,th0,rD,U], rD*BRxB/rRxB ] * Sin[phiDet]

xExBGCB[xDet_, p_, th0_, phiDet_, BRxB_, rRxB_, rD_, U_]:= xDetGCB[xDet, p, th0, phiDet, BRxB, rRxB, rD, U]/Sqrt[rRxB/rD]
yExBGCB[yDet_, p_, th0_, phiDet_, BRxB_, rRxB_, rD_, U_]:= yDetGCB[yDet, p, th0, phiDet, BRxB, rRxB, rD, U]/Sqrt[rRxB/rD]

(*xRxBGCB[]:= xExBGCB[] + kExB[p, ExBscale] * yAGCB*)

(* ::Section:: *)
(*hand calculated xDV yDV*)

yDVHandSolve[phiDV_, xDet_, yDet_, p_, th0_, BRxB_, rRxB_, rA_, rD_, ExBscale_, U_, phiDet_] := 
(yExBGCB[yDet, p, th0, phiDet, BRxB, rRxB, rD, U] - 
      kExB[p, ExBscale]*
       xExBGCB[xDet, p, th0, phiDet, BRxB, rRxB, rD, U])/(1 + 
      kExB[p, ExBscale]^2)/Sqrt[1/rA] + 
  rG[p, th0, BRxB/rRxB]*Sin[phiDV]
  
xDVHandSolve[phiDV_, xDet_, yDet_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, R_, G1_, G2_, ExBscale_, U_, phiDet_] :=
 (
    xExBGCB[xDet, p, th0, phiDet, BRxB, rRxB, rD, U]
     + kExB[p, ExBscale]/Sqrt[1/rA]*
      (
       yDVHandSolve[phiDV, xDet, yDet, p, th0, BRxB, rRxB, rA, rD, 
         ExBscale, U, phiDet] - rG[p, th0, BRxB/rRxB]*Sin[phiDV]
       )
     - D1stBPolyGrad[p, alpha, th0, BRxB, rRxB, 
      yGCAF[yDVHandSolve[phiDV, xDet, yDet, p, th0, BRxB, rRxB, rA, 
        rD, ExBscale, U, phiDet], phiDV, p, th0, BRxB, rRxB, rA], R, 
      G1, G2]
    )/Sqrt[1/rA] +
  rG[p, th0, BRxB/rRxB]*Cos[phiDV]

 (*xA for Aperture function*)
 
 
xAAccel[phiA_, xDet_, yDet_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_,R_, G1_, G2_, ExBscale_, U_] := 
	xGCAF[xDVHandSolve[phiDV, xDet, yDet, p, th0, alpha, BRxB, rRxB, rA, rD, R, G1, G2, ExBscale, U, phiDet], phiDV, p, th0, BRxB, rRxB, rA] + 
	rG[p, th0, rA*BRxB/rRxB]*Cos[phiA]
	
yAAccel[phiA_, xDet_, yDet_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_,R_, G1_, G2_, ExBscale_, U_] := 
	yGCAF[yDVHandSolve[phiDV, xDet, yDet, p, th0, BRxB, rRxB, rA, rD, ExBscale, U, phiDet], phiDV, p, th0, BRxB, rRxB, rA] + 
	rG[p, th0, rA*BRxB/rRxB]*Sin[phiA]
  
  
  
  
(*Compiled Integrand*)

IntegrandProton2DwNBeamAccelCompiled[
	lambda_, a_, b_, xDet_, yDet_, {phiDV_?NumericQ, phiA_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, 
	{alpha_, BRxB_, rRxB_, rA_,rD_, R_, G1_, G2_, ExBscale_,U_}, 
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] :=
 
 pmomNormedGlueck[p, lambda, a, b]*Sin[th0]*
  
  TrapezNBeamCompiledNormed[
   xDVHandSolve[phiDV, xDet, yDet, p, th0, alpha, BRxB, rRxB, rA, rD, R, G1, G2, ExBscale, U, phiDet]-xOff, {twx, plx, k1x, k2x, k3x}]*
  TrapezNBeamCompiledNormed[
   yDVHandSolve[phiDV, xDet, yDet, p, th0, BRxB, rRxB, rA, rD, ExBscale, U, phiDet]-yOff, {twy, ply,k1y, k2y, k3y}]*
  
  ApertBooleCompiled[
   xAAccel[phiA, xDet, yDet, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, ExBscale, U], 
   yAAccel[phiA, xDet, yDet, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, ExBscale, U], 
   	xAA, yAA, xOff, yOff] 
  
  
(* ::Section:: *)
(* get limits of integration *)

phiAminApertAccel[xDet_, yDet_, plocal_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_, ExBscale_,U_] = 
  Solve[xOff + xAA/2 == xAAccel[phiAlocal, xDet, yDet, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, ExBscale, U], phiAlocal][[All, 1, 2, 1, 1]];
phiAmaxApertAccel[xDet_, yDet_, plocal_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_, ExBscale_,U_] = 
   Solve[xOff - xAA/2 == xAAccel[phiAlocal, xDet, yDet, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, ExBscale, U], phiAlocal][[All, 1, 2, 1, 1]];


phiACasesAccel[yD_?NumericQ, xD_?NumericQ, plocal_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_, ExBscale_,U_] :=
 {phiA,
 	Sequence@@Re[{
  phiAmaxApertAccel[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA, ExBscale, U][[1]],
  phiAminApertAccel[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA, ExBscale, U][[1]],
  phiAminApertAccel[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA, ExBscale, U][[2]],
  phiAmaxApertAccel[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA, ExBscale, U][[2]]
  }]
 }

(*phiA limits from Y Aperture*)

phiAminApertfromYAccel[xDet_, yDet_, plocal_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_, R_, G1_, G2_, ExBscale_, U_] = 
 Solve[yOff + yAA/2 == yAAccel[phiAlocal, xDet, yDet, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, ExBscale, U], phiAlocal][[All, 1, 2, 1, 1 ;; -2]]
 
phiAmaxApertfromYAccel[xDet_, yDet_, plocal_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_, R_, G1_, G2_, ExBscale_, U_] = 
 Solve[yOff - yAA/2 == yAAccel[phiAlocal, xDet, yDet, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, ExBscale, U], phiAlocal][[All, 1, 2, 1, 1 ;; -2]]

  
(*new manual phiA integration with all phi limits*)
ManualphiAIntegrandAllLimitsAccel[
	lambda_, a_, b_, xDet_, yDet_, {phiDV_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, 
	{alpha_, BRxB_, rRxB_, rA_, rD_, R_, G1_, G2_, ExBscale_, U_}, {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] := 
	Module[
  		{
  			phimin = Re[phiAminApertAccel[xDet, yDet, If[p==0.,0.00000001,p], th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, ExBscale, U]],
  			phimax = Re[phiAmaxApertAccel[xDet, yDet, If[p==0.,0.00000001,p], th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, ExBscale, U]], 
  			phiminY = Re[phiAminApertfromYAccel[xDet, yDet, If[p==0.,0.00000001,p], th0, alpha, BRxB, rRxB, rA, rD, phiDet, yOff, yAA, R, G1, G2, ExBscale, U]], 
  			phimaxY = Re[phiAmaxApertfromYAccel[xDet, yDet, If[p==0.,0.00000001,p], th0, alpha, BRxB, rRxB, rA, rD, phiDet, yOff, yAA, R, G1, G2, ExBscale, U]],
			IntValues, philimitlist, philimitlistSorted, Deltalimits, phiCenters, Integrated
		},
  		
  		philimitlist = Flatten[{-Pi//N, phimin, phimax, phiminY, phimaxY, Pi//N}];
		philimitlist = If[# > Pi, # - 2 Pi, #] & /@ philimitlist;
		philimitlistSorted = Union[philimitlist];
		Deltalimits = Table[philimitlistSorted[[i + 1]] - philimitlistSorted[[i]], {i, 1, Length[philimitlistSorted] - 1}];
		phiCenters = Table[philimitlistSorted[[i]] + Deltalimits[[i]]/2, {i, 1, Length[philimitlistSorted] - 1}];
		IntValues = IntegrandProton2DwNBeamAccelCompiled[
			lambda, a, b, xDet, yDet, {phiDV, #, phiDet, p, th0}, 
			{alpha, BRxB, rRxB, rA,rD, R, G1, G2, ExBscale,U}, 
			{twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, 
			{xAA, yAA, xOff, yOff}] & /@ phiCenters;
		Integrated = Total[Deltalimits*IntValues];
  		Integrated
  		
  ]

  
  