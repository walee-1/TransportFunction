(* Wolfram Language package *)



(* Wolfram Language package *)
Get["Spectra/ElectronSpectrum.m"];
Get["Bfield+Drift/Bfields.m"]
Get["Bfield+Drift/Drifts.m"];
Get["Common/CommonFunctions.m"];
Get["NeutronBeam/NeutronBeamDefinition.m"];
Get["Aperture/BooleApert.m"];


DY[p_,alpha_,BRxB_,th0_,rRxB_,R_]:=2*D1stSimple[p, alpha, th0, BRxB, rRxB]^2/Pi^2/R
   
charge= -1;   
   
xGCAShift[xDV_, phiDV_, p_, th0_, BRxB_, rRxB_, rA_, XAShift_] := (xDV - rG[p, th0, BRxB/rRxB]*Cos[phiDV])*Sqrt[1/rA] + XAShift
yGCAShift[yDV_, phiDV_, p_, th0_, BRxB_, rRxB_, rA_, YAShift_] := (yDV - rG[p, th0, BRxB/rRxB]*Sin[phiDV])*Sqrt[1/rA] + YAShift
yRxBGCShift[yDV_, phiDV_, p_, th0_, BRxB_, rRxB_, YRxBShift_]:= (yDV - rG[p, th0, BRxB/rRxB]*Sin[phiDV]) *Sqrt[1/rRxB] + YRxBShift  
  
DeltayDVShift[phiDV_, yD_, p_,alpha_, th0_, BRxB_, rRxB_, rD_, phiDet_, YDetShift_,R_] := 
(	
	(yD - rG[p, th0, rD*BRxB/rRxB]*Sin[phiDet] - YDetShift)*Sqrt[rD/rRxB] -
	DY[p,alpha,BRxB,th0,rRxB,R]
)*Sqrt[rRxB] + rG[p, th0, BRxB/rRxB]*Sin[phiDV]

DeltaxDVShift[phiDV_, yDV_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rD_, phiDet_, R_, G1_, G2_, YRxBShift_, XDetShift_] := 

(
	(xD - rG[p, th0, rD*BRxB/rRxB]*Cos[phiDet] - XDetShift)*Sqrt[rD/rRxB] -
   	charge*D1stBPolyGrad[p, alpha, th0, BRxB, rRxB, yRxBGCShift[yDV, phiDV, p, th0, BRxB, rRxB, YRxBShift], R, G1, G2]
)*Sqrt[rRxB] + rG[p, th0, BRxB/rRxB]*Cos[phiDV]
   
   
DeltaxDV2Shift[phiDV_, yD_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rD_, phiDet_, R_, G1_, G2_, {YRxBShift_, XDetShift_, YDetShift_}] := 
 	DeltaxDVShift[phiDV, 
  		DeltayDVShift[phiDV, yD, p, alpha, th0, BRxB, rRxB, rD, phiDet, YDetShift, R],
  		xD, p, th0, alpha, BRxB, rRxB, rD, phiDet, R, G1, G2, YRxBShift, XDetShift
  	]
   
   
xAShift[phiA_, yD_, xD_, p_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, {XAShift_, YRxBShift_, XDetShift_, YDetShift_}] := 
 	xGCAShift[
 		DeltaxDV2Shift[phiDV, yD, xD, p, th0, alpha, BRxB, rRxB, rD, phiDet, R, G1, G2, {YRxBShift, XDetShift, YDetShift}], 
 		phiDV, p, th0, BRxB, rRxB, rA, XAShift
 	] + rG[p, th0, rA*BRxB/rRxB]*Cos[phiA]
  
yAShift[phiA_, yD_, p_, alpha_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, YAShift_, YDetShift_, R_] := 
 	yGCAShift[
 		DeltayDVShift[phiDV, yD, p, alpha, th0, BRxB, rRxB, rD, phiDet, YDetShift, R], 
 		phiDV, p, th0, BRxB, rRxB, rA, YAShift
 	] + rG[p, th0, rA*BRxB/rRxB]*Sin[phiA]
   
   
   
   
(*Integrand2DwNBeamv31[b_, yD_, 
  xD_, {phiDV_?NumericQ, phiA_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, {alpha_, BRxB_, rRxB_, rA_,
    rD_, R_, G1_, G2_}, {twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, 
   k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] :=
 
 wmomNormedWb[\[Lambda]0, \[Kappa]0, b, p]*Sin[th0]*
  
  TrapezNBeamNewNormed[
   	DeltaxDV2[phiDV, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2]-xOff,
    {twx, plx, k1x, k2x, k3x}]*
  TrapezNBeamNewNormed[
   	DeltayDV[phiDV, yD, p, th0, BRxB, rRxB, rA, rD, phiDet]-yOff, 
   	{twy, ply, k1y, k2y, k3y}]*
  
  ApertBoole[
   xA[phiA, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2], 
   yA[phiA, yD, p, th0, BRxB, rRxB, rA, rD, phiDet], 
   {xAA, yAA, xOff, yOff}
  ]   *)
   
   

(*Compiled Integrand*)

Integrand2DwNBeamCompiledShift[
	b_, yD_, xD_, {phiDV_?NumericQ, phiA_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, {alpha_, BRxB_, rRxB_, rA_,rD_, R_, G1_, G2_},
	{XAShift_, YAShift_, YRxBShift_, XDetShift_, YDetShift_}, 
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] :=
 
 wmomNormedWb[\[Lambda]0, \[Kappa]0, b, p]*Sin[th0]*
  
  TrapezNBeamCompiledNormed[
   DeltaxDV2Shift[phiDV, yD, xD, p, th0, alpha, BRxB, rRxB, rD, phiDet, R, G1, G2, {YRxBShift, XDetShift, YDetShift}]-xOff, {twx, plx, k1x, k2x, k3x}]*
  TrapezNBeamCompiledNormed[
   DeltayDVShift[phiDV, yD, p, alpha, th0, BRxB, rRxB, rD, phiDet, YDetShift, R]-yOff, {twy, ply,k1y, k2y, k3y}]*
  
  ApertBooleCompiled[
   xAShift[phiA, yD, xD, p, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, {XAShift, YRxBShift, XDetShift, YDetShift}], 
   yAShift[phiA, yD, p, alpha, th0, BRxB, rRxB, rA, rD, phiDet, YAShift, YDetShift, R], 
   	xAA, yAA, xOff, yOff]  



(* ::Section:: *)
(* Integrand with phiA limits and compiled integrand *)




phiAminApertShift[yD_, xD_, plocal_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_,{XAShift_, YRxBShift_, XDetShift_, YDetShift_}] = 
  Solve[xOff + xAA/2 == xAShift[phiAlocal, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, {XAShift, YRxBShift, XDetShift, YDetShift}], phiAlocal][[All, 1, 2, 1, 1]];
phiAmaxApertShift[yD_, xD_, plocal_, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_,{XAShift_, YRxBShift_, XDetShift_, YDetShift_}] = 
   Solve[xOff - xAA/2 == xAShift[phiAlocal, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, {XAShift, YRxBShift, XDetShift, YDetShift}], phiAlocal][[All, 1, 2, 1, 1]];


(*phiACases[yD_?NumericQ, xD_?NumericQ, plocal_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_] :=
 {phiA,
 	Sequence@@Re[{
  phiAmaxApert[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[1]],
  phiAminApert[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[1]],
  phiAminApert[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[2]],
  phiAmaxApert[yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet,R, G1, G2, xOff, xAA][[2]]
  }]
 }*)



(*Integrand with manual summation of phiA integration*)
(*ManualphiAIntegrand[b_, yD_, xD_, {phiDV_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, 
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
  ]*)
  

(*phiA limits from Y Aperture*)

phiAminApertfromYShift[yD_, plocal_, alpha_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_, YAShift_, YDetShift_, R_] = 
 Solve[yOff + yAA/2 == yAShift[phiAlocal, yD, plocal, alpha, th0, BRxB, rRxB, rA, rD, phiDet, YAShift, YDetShift, R], phiAlocal][[All, 1, 2, 1, 1 ;; -2]]
 
phiAmaxApertfromYShift[yD_, plocal_, alpha_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_, YAShift_, YDetShift_,R_] = 
 Solve[yOff - yAA/2 == yAShift[phiAlocal, yD, plocal, alpha, th0, BRxB, rRxB, rA, rD, phiDet, YAShift, YDetShift, R], phiAlocal][[All, 1, 2, 1, 1 ;; -2]]

  
(*new manual phiA integration with all phi limits*)
ManualphiAIntegrandAllLimitsShift[
	b_, yD_, xD_, {phiDV_?NumericQ, phiDet_?NumericQ, p_?NumericQ, th0_?NumericQ}, 
	{alpha_, BRxB_, rRxB_, rA_, rD_, R_, G1_, G2_}, {XAShift_, YAShift_, YRxBShift_, XDetShift_, YDetShift_}, 
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_}] := 
	Module[
  		{
  			phimin = Re[phiAminApertShift[yD, xD, If[p==0.,0.00000001,p], th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}]],
  			phimax = Re[phiAmaxApertShift[yD, xD, If[p==0.,0.00000001,p], th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}]], 
  			phiminY = Re[phiAminApertfromYShift[yD, If[p==0.,0.00000001,p], alpha, th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA, YAShift, YDetShift, R]], 
  			phimaxY = Re[phiAmaxApertfromYShift[yD, If[p==0.,0.00000001,p], alpha, th0, BRxB, rRxB, rA, rD, phiDet, yOff, yAA, YAShift, YDetShift, R]],
			IntValues, philimitlist, philimitlistSorted, Deltalimits, phiCenters, Integrated
		},
  		
  		philimitlist = Flatten[{-Pi//N, phimin, phimax, phiminY, phimaxY, Pi//N}];
		philimitlist = If[# > Pi, # - 2 Pi, #] & /@ philimitlist;
		philimitlistSorted = Union[philimitlist];
		Deltalimits = Table[philimitlistSorted[[i + 1]] - philimitlistSorted[[i]], {i, 1, Length[philimitlistSorted] - 1}];
		phiCenters = Table[philimitlistSorted[[i]] + Deltalimits[[i]]/2, {i, 1, Length[philimitlistSorted] - 1}];
		IntValues = Integrand2DwNBeamCompiledShift[b, yD, xD, {phiDV, #, phiDet, p, th0}, {alpha, BRxB, rRxB, rA, rD, R, G1,G2}, 
			{XAShift, YAShift, YRxBShift, XDetShift, YDetShift},
			{twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA, yAA, xOff, yOff}] & /@ phiCenters;
		Integrated = Total[Deltalimits*IntValues];
  		Integrated
  		
  ]  




(* ::Section:: *)
(* Integrand with p limits from aperture *)


pminApertShift[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_,{XAShift_, YRxBShift_, XDetShift_, YDetShift_}] :=
Module[
	{plocal}, 
 	NSolve[xOff + xAA/2 == xAShift[phiA, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1,G2, {XAShift, YRxBShift, XDetShift, YDetShift}], plocal][[All, 1, 2]]
]
      
      
pminApertCasesShift[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_,{XAShift_, YRxBShift_, XDetShift_, YDetShift_}] := 
Module[
  {pmin = Cases[pminApertShift[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}], _?(# \[Element] Reals && 0. <= # <= 2*pmax &)]},
  Piecewise[
   {
    {0, Length[pmin] == 0},
    {pmax, Length[pmin] == 1 && pmin[[1]] > pmax},
    {pmax, Length[pmin] == 2 && pmin[[1]] > pmax && pmin[[2]] > pmax},
    {pmin[[1]], Length[pmin] == 1},
    {
    	Print[
    		"pmin: Length =",
    		Length[pmin],
    		pmin,phiA," ",th0," ",phiDet," ",
    		pminApertShift[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}]
    	], 
    	Length[pmin]>1(*Length[pmin]!=1 && Length[pmin]!=0*) 
    }
    },
    Print["pmin: other Case"]
   ]
 
]
  
pmaxApertShift[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_, R_, G1_, G2_, xOff_, xAA_,{XAShift_,YRxBShift_, XDetShift_, YDetShift_}] :=
	Module[
		{plocal}, 
 		NSolve[xOff - xAA/2 == xAShift[phiA, yD, xD, plocal, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1,G2, {XAShift, YRxBShift, XDetShift, YDetShift}], plocal][[All, 1, 2]]
	]


pmaxApertCasesShift[phiA_?NumericQ, yD_?NumericQ, xD_?NumericQ, th0_, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_,{XAShift_, YRxBShift_, XDetShift_, YDetShift_}] := 
Module[
  {pmaxx =Cases[pmaxApertShift[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}], _?(# \[Element] Reals && 0. <= # <= 2*pmax &)]},
  
  Piecewise[
   {
    {pmax, Length[pmaxx] == 0},
    {pmax, Length[pmaxx] == 1 && pmaxx[[1]] > pmax},
    {pmax, Length[pmaxx] == 2 && pmaxx[[1]] > pmax && pmaxx[[2]] > pmax},
    {pmaxx[[1]], Length[pmaxx] == 1},
    {	
    	Print["pmax: Length =",Length[pmaxx],pmaxx,phiA," ",th0," ",phiDet," ",
    		pmaxApertShift[phiA, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}]
    	],
    	Length[pmaxx]>1(*Length[pmaxx]!=1 && Length[pmaxx]!=0*) 
    }
    },
    Print["pmax: other case"]
   ]

  ]
  
  
pLimitsApertXListShift[yD_?NumericQ, xD_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_,
	{XAShift_, YRxBShift_, XDetShift_, YDetShift_}] := 
(*pLimitsApertXList[yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA] =*)
Module[
  {
   pminXmPi = pminApertCasesShift[-Pi, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}],
   pminX0 = pminApertCasesShift[0., yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}],
   pmaxXmPi = pmaxApertCasesShift[-Pi, yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}],
   pmaxX0 = pmaxApertCasesShift[0., yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}],
   plimitlist
   },
	plimitlist=Union[{pminXmPi, pminX0, pmaxXmPi, pmaxX0}];
	Which[
		Length[plimitlist]==1,Join[plimitlist,plimitlist],
		True,plimitlist
	]
]
  


(*p limits from apertY*)

(*pminApertY[phiA_, yD_, th0_, BRxB_, rRxB_, rA_, rD_, phiDet_, yOff_, yAA_] =
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
  ]*)



(* ::Section:: *)
(* 2 step-Integration *)

pXDomainShift[yD_?NumericQ, xD_?NumericQ, th0_?NumericQ, alpha_, BRxB_, rRxB_, rA_, rD_, phiDet_?NumericQ, R_, G1_, G2_, xOff_, xAA_,
	{XAShift_, YRxBShift_, XDetShift_, YDetShift_}]:=
{p,Sequence@@pLimitsApertXListShift[yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}]}


Step1IntShift[b_, yD_, xD_, {phiDet_?NumericQ, th0_?NumericQ}, 
	{alpha_, BRxB_, rRxB_, rA_, rD_, R_, G1_, G2_},{XAShift_, YAShift_, YRxBShift_, XDetShift_, YDetShift_}, 
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_}, {xAA_, yAA_, xOff_, yOff_},
	method1_,PrecGoal1_,AccGoal1_,MinRec1_,MaxRec1_]:=
	NIntegrate[
		ManualphiAIntegrandAllLimitsShift[b, yD, xD, {phiDV, phiDet, p, th0}, 
			{alpha, BRxB, rRxB, rA, rD, R, G1, G2}, {XAShift, YAShift, YRxBShift, XDetShift, YDetShift}, {twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA, yAA, xOff, yOff}],
		Evaluate[pXDomainShift[yD, xD, th0, alpha, BRxB, rRxB, rA, rD, phiDet, R, G1, G2, xOff, xAA, {XAShift, YRxBShift, XDetShift, YDetShift}]],
		{phiDV,-Pi,Pi},
		Method->method1,PrecisionGoal->PrecGoal1,AccuracyGoal->AccGoal1, MinRecursion->MinRec1, MaxRecursion->MaxRec1	
]


Step2IntShift[b_, xD_?NumericQ,yD_?NumericQ, {alpha_, BRxB_, rF_, rRxB_, rA_, rD_, R_, G1_, G2_},{XAShift_, YAShift_, YRxBShift_, XDetShift_, YDetShift_}, 
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_},
	{xAA_, yAA_, xOff_, yOff_}, method1_,PrecGoal1_,AccGoal1_,MinRec1_,MaxRec1_, method2_,PrecGoal2_,AccGoal2_,MinRec2_,MaxRec2_]:=

		NIntegrate[
		Step1IntShift[b, yD, xD, {phiDet, th0}, 
			{alpha, BRxB, rRxB, rA, rD, R, G1, G2}, {XAShift, YAShift, YRxBShift, XDetShift, YDetShift}, 
			{twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA, yAA, xOff, yOff},
			method1,PrecGoal1,AccGoal1,MinRec1,MaxRec1],
			{th0,0.,thetamax[rF]},{phiDet,-Pi,Pi},Method->method2,PrecisionGoal->PrecGoal2,AccuracyGoal->AccGoal2,MinRecursion->MinRec2, MaxRecursion->MaxRec2		
		]
		
		
(* ::Section:: *)
(* Add xD and yD integration *)

BinIntShift[b_, OneBinList_, {alpha_, BRxB_, rF_, rRxB_, rA_, rD_, R_, G1_, G2_},{XAShift_, YAShift_, YRxBShift_, XDetShift_, YDetShift_}, 
	{twx_, plx_, k1x_, k2x_, k3x_, twy_, ply_, k1y_, k2y_, k3y_},
	{xAA_, yAA_, xOff_, yOff_}, 
	{method1_,PrecGoal1_,AccGoal1_,MinRec1_,MaxRec1_}, 
	{method2_,PrecGoal2_,AccGoal2_,MinRec2_,MaxRec2_}, 
	{method3_,PrecGoal3_,AccGoal3_,MinRec3_,MaxRec3_}]:=
	NIntegrate[
		
		Step2IntShift[b, xD, yD, {alpha, BRxB, rF, rRxB, rA, rD, R, G1, G2}, {XAShift, YAShift, YRxBShift, XDetShift, YDetShift}, 
			{twx, plx, k1x, k2x, k3x, twy, ply, k1y, k2y, k3y}, {xAA, yAA, xOff, yOff},
			method1,PrecGoal1,AccGoal1,MinRec1,MaxRec1,method2,PrecGoal2,AccGoal2,MinRec2,MaxRec2],
			
		{xD, OneBinList[[1,1]], OneBinList[[1,2]]}, {yD, OneBinList[[2,1]], OneBinList[[2,2]]},
		PrecisionGoal->PrecGoal3,Method->method3,MinRecursion->MinRec3,MaxRecursion->MaxRec3,AccuracyGoal->AccGoal3
	]


     
