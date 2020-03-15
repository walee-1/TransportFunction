(* Wolfram Language package *)
Get["Spectra/ElectronSpectrum.m"];
Get["Bfield+Drift/Bfields.m"]
Get["Bfield+Drift/Drifts.m"];
Get["Common/CommonFunctions.m"];
Get["NeutronBeam/NeutronBeamDefinition.m"];
Get["Aperture/BooleApert.m"];
   
   
   
   
xGCA[xDV_, phiDV_, p_, th0_, BRxB_, rRxB_, 
  rA_] := (xDV - rG[p, th0, BRxB/rRxB]*Cos[phiDV])*Sqrt[1/rA]
yGCA[yDV_, phiDV_, p_, th0_, BRxB_, rRxB_, 
  rA_] := (yDV - rG[p, th0, BRxB/rRxB]*Sin[phiDV])*Sqrt[1/rA]
  
  
DeltayDV[phiDV_, yD_, p_, th0_, BRxB_, rRxB_, rA_, rD_, 
  phiDet_] := (yD - rG[p, th0, rD*BRxB/rRxB]*Sin[phiDet])*
   Sqrt[rD/rRxB]*Sqrt[rA] + rG[p, th0, BRxB/rRxB]*Sin[phiDV]
DeltaxDV[phiDV_, yDV_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, 
  phiDet_, R_, G1_, 
  G2_] := ((xD - rG[p, th0, rD*BRxB/rRxB]*Cos[phiDet])*Sqrt[rD/rRxB] -
      D1stBPolyGrad[p, alpha, th0, BRxB, rRxB, 
      yGCA[yDV, phiDV, p, th0, BRxB, rRxB, rA], R, G1, G2])*Sqrt[rA] +
   rG[p, th0, BRxB/rRxB]*Cos[phiDV]
   
   
DeltaxDV2[phiDV_, yD_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, 
  phiDet_, R_, G1_, G2_] := 
 DeltaxDV[phiDV, 
  DeltayDV[phiDV, yD, p, th0, BRxB, rRxB, rA, rD, phiDet], xD, p, th0,
   alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2]
   
   
xA[phiA_, yD_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_,
   R_, G1_, G2_] := 
 xGCA[DeltaxDV2[phiDV, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, 
    phiDet, R, G1, G2], phiDV, p, th0, BRxB, rRxB, rA] + 
  rG[p, th0, rA*BRxB/rRxB]*Cos[phiA]
  
yA[phiA_, yD_, p_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_] := 
 yGCA[DeltayDV[phiDV, yD, p, th0, BRxB, rRxB, rA, rD, phiDet], phiDV, 
   p, th0, BRxB, rRxB, rA] + rG[p, th0, rA*BRxB/rRxB]*Sin[phiA]
   
   
   
   
Integrand2DwNBeamv31[b_, yD_, 
  xD_, {phiDV_?NumericQ, phiA_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, {alpha_, BRxB_, rRxB_, rA_,
    rD_, R_, G1_, G2_}, {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, 
   k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] :=
 
 wmomNormedWb[\[Lambda]0, \[Kappa]0, b, p]*Sin[th0]*
  
  TrapezNBeamNewNormed[
   DeltaxDV2[phiDV, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,
     R, G1, G2]-xOff, {twx, plx, k1x, k2x, k3x}]*
  TrapezNBeamNewNormed[
   DeltayDV[phiDV, yD, p, th0, BRxB, rRxB, rA, rD, phiDet]-yOff, {twy, ply,
     k1y, k2y, k3y}]*
  
  ApertBoole[
   xA[phiA, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, 
    G2], yA[phiA, yD, p, th0, BRxB, rRxB, rA, rD, phiDet], {xAA, yAA, 
    xOff, yOff}]   
   
   

(*Compiled Integrand*)

Integrand2DwNBeamCompiled[
	b_, yD_, xD_, {phiDV_?NumericQ, phiA_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, {alpha_, BRxB_, rRxB_, rA_,rD_, R_, G1_, G2_}, 
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] :=
 
 wmomNormedWb[\[Lambda]0, \[Kappa]0, b, p]*Sin[th0]*
  
  TrapezNBeamCompiledNormed[
   DeltaxDV2[phiDV, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2]-xOff, {twx, plx, k1x, k2x, k3x}]*
  TrapezNBeamCompiledNormed[
   DeltayDV[phiDV, yD, p, th0, BRxB, rRxB, rA, rD, phiDet]-yOff, {twy, ply,k1y, k2y, k3y}]*
  
  ApertBooleCompiled[
   xA[phiA, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2], yA[phiA, yD, p, th0, BRxB, rRxB, rA, rD, phiDet], 
   	xAA, yAA, xOff, yOff]  



(* ::Section:: *)
(* Integrand with phiA limits and compiled integrand *)




phiAminApert[yD_, xD_, plocal_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_] = 
  Solve[xOff + xAA/2 == xA[phiAlocal, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2], phiAlocal][[All, 1, 2, 1, 1]];
phiAmaxApert[yD_, xD_, plocal_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_] = 
   Solve[xOff - xAA/2 == xA[phiAlocal, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2], phiAlocal][[All, 1, 2, 1, 1]];


phiACases[yD_?NumericQ, xD_?NumericQ, plocal_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_] :=
 {phiA,
 	Sequence@@Re[{
  phiAmaxApert[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[1]],
  phiAminApert[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[1]],
  phiAminApert[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[2]],
  phiAmaxApert[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[2]]
  }]
 }



(*Integrand with manual summation of phiA integration*)
ManualphiAIntegrand[b_, yD_, xD_, {phiDV_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, 
	{alpha_, BRxB_, rRxB_, rA_, rD_, R_, G1_, G2_}, {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] := 
	Module[
  		{phimin2, phimax2, phimean, int5Dvalue},
  		phimin2 = Re[phiAminApert[yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA]];
 		phimax2 = Re[phiAmaxApert[yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA]];
  		phimean = (phimax2[[1]] + phimin2[[1]])/2;
  		int5Dvalue = Integrand2DwNBeamCompiled[b, yD, xD, {phiDV, phimean, phiDet, p, th0}, {alpha, BRxB, rRxB, rA, rD, 
     R, G1, G2}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, 
     k3y}, {xAA, yAA, xOff, yOff}];
  		(phimin2[[1]] - phimax2[[1]])* int5Dvalue + (phimax2[[2]] - phimin2[[2]])*int5Dvalue
  ]
  

(*phiA limits from Y Aperture*)

phiAminApertfromY[yD_, plocal_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_] = 
 Solve[yOff + yAA/2 == yA[phiAlocal, yD, plocal, th0, BRxB, rRxB, rA, rD, phiDet], phiAlocal][[All, 1, 2, 1, 1 ;; -2]]
 
phiAmaxApertfromY[yD_, plocal_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_] = 
 Solve[yOff - yAA/2 == yA[phiAlocal, yD, plocal, th0, BRxB, rRxB, rA, rD, phiDet], phiAlocal][[All, 1, 2, 1, 1 ;; -2]]

  
(*new manual phiA integration with all phi limits*)
ManualphiAIntegrandAllLimits[
	b_, yD_, xD_, {phiDV_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, 
	{alpha_, BRxB_, rRxB_, rA_, rD_, R_, G1_, G2_}, {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] := 
	Module[
  		{
  			phimin = Re[phiAminApert[yD, xD, If[p==0.,0.00000001,p], th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA]],
  			phimax = Re[phiAmaxApert[yD, xD, If[p==0.,0.00000001,p], th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA]], 
  			phiminY = Re[phiAminApertfromY[yD, If[p==0.,0.00000001,p], th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA]], 
  			phimaxY = Re[phiAmaxApertfromY[yD, If[p==0.,0.00000001,p], th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA]],
			IntValues, philimitlist, philimitlistSorted, Deltalimits, phiCenters, Integrated
		},
  		
  		philimitlist = Flatten[{-Pi//N, phimin, phimax, phiminY, phimaxY, Pi//N}];
		philimitlist = If[# > Pi, # - 2 Pi, #] & /@ philimitlist;
		philimitlistSorted = DeleteDuplicates[Sort[philimitlist]];
		Deltalimits = Table[philimitlistSorted[[i + 1]] - philimitlistSorted[[i]], {i, 1, Length[philimitlistSorted] - 1}];
		phiCenters = Table[philimitlistSorted[[i]] + Deltalimits[[i]]/2, {i, 1, Length[philimitlistSorted] - 1}];
		IntValues = Integrand2DwNBeamCompiled[b, yD, xD, {phiDV, #, phiDet, p, th0}, {alpha, BRxB, rRxB, rA, rD, R, G1,G2},
			{twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA, yAA, xOff, yOff}] & /@ phiCenters;
		Integrated = Total[Deltalimits*IntValues];
  		Integrated
  		
  ]  




(* ::Section:: *)
(* Integrand with p limits from aperture *)


pminApert[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_] :=
Module[
	{plocal}, 
 	NSolve[xOff + xAA/2 == xA[phiA, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1,G2], plocal][[All, 1, 2]]
]
      
      
pminApertCases[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_] := 
Module[
  {pmin = Cases[pminApert[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA], _?(# \[Element] Reals && 0. <= # <= 2*pmax &)]},
  Piecewise[
   {
    {0, Length[pmin] == 0},
    {pmax, Length[pmin] == 1 && pmin[[1]] > pmax},
    {pmin[[1]], Length[pmin] == 1},
    {
    	Print["pmin: Length =",Length[pmin],pmin,phiA," ",th0," ",phiDet," ",pminApert[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA]], 
    	Length[pmin]>1(*Length[pmin]!=1 && Length[pmin]!=0*) 
    }
    },
    Print["pmin: other Case"]
   ]
 
]
  
pmaxApert[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_] :=
	Module[
		{plocal}, 
 		NSolve[xOff - xAA/2 == xA[phiA, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1,G2], plocal][[All, 1, 2]]
	]


pmaxApertCases[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_] := 
Module[
  {pmaxx =Cases[pmaxApert[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA], _?(# \[Element] Reals && 0. <= # <= 2*pmax &)]},
  
  Piecewise[
   {
    {pmax, Length[pmaxx] == 0},
    {pmax, Length[pmaxx] == 1 && pmaxx[[1]] > pmax},
    {pmaxx[[1]], Length[pmaxx] == 1},
    {	
    	Print["pmax: Length =",Length[pmaxx],pmaxx,phiA," ",th0," ",phiDet," ",pmaxApert[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA]],
    	Length[pmaxx]>1(*Length[pmaxx]!=1 && Length[pmaxx]!=0*) 
    }
    },
    Print["pmax: other case"]
   ]

  ]
  
  
pLimitsApertXList[yD_?NumericQ, xD_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_] := 
(*pLimitsApertXList[yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA] =*)
Module[
  {
   pminXmPi = pminApertCases[-Pi, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA],
   pminX0 = pminApertCases[0., yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA],
   pmaxXmPi = pmaxApertCases[-Pi, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA],
   pmaxX0 = pmaxApertCases[0., yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA]
   },
 (*Sequence@@*)Union[{pminXmPi, pminX0, pmaxXmPi, pmaxX0}]
  ]


(*p limits from apertY*)

pminApertY[phiA_, yD_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_] =
 Module[{plocal}, Solve[yOff + yAA/2 == yA[phiA, yD, plocal, th0, BRxB, rRxB, rA, rD, phiDet], plocal, Reals][[1, 1, 2, 1]]]
 
pmaxApertY[phiA_, yD_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_] =
Module[{plocal}, Solve[yOff - yAA/2 == yA[phiA, yD, plocal, th0, BRxB, rRxB, rA, rD, phiDet], plocal][[1, 1, 2]]]

pminApertYCases[yD_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_] := Module[
  {phiAlocal = If[yD > yAA/2 + yOff, -Pi/2, Pi/2], pminlocal},
  pminlocal = pminApertY[phiAlocal, yD, th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA];
  Piecewise[{{pminlocal, 0. < pminlocal < pmax}}, 0.]
  ]

pmaxApertYCases[yD_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_] := Module[
  {phiAlocal = If[yD < -yAA/2 + yOff, Pi/2, -Pi/2], pmaxlocal},
  pmaxlocal = pmaxApertY[phiAlocal, yD, th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA];
  Piecewise[{{pmaxlocal, 0. < pmaxlocal < pmax}}, 0.]
  ]

pLimitsApertAllList[yD_?NumericQ, xD_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, {xAA_, xOff_,yAA_, yOff_}] := 
Module[
  {
   pminXmPi = pminApertCases[-Pi, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA],
   pminX0 = pminApertCases[0., yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA],
   pmaxXmPi = pmaxApertCases[-Pi, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA],
   pmaxX0 = pmaxApertCases[0., yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA],
   pminY = pminApertYCases[yD, th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA],
   pmaxY = pmaxApertYCases[yD, th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA]
   },
 	Sort[{pminXmPi, pminX0, pmaxXmPi, pmaxX0,pminY,pmaxY}]
  ]



     
   