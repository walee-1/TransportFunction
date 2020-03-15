(* Wolfram Language package *)

Get["Common/Constants.m"];

Eofp[p_] := Sqrt[me^2 + p^2]

TofPClassical[p_] := p^2/(2*mp);

pofT[T_, m_] := Sqrt[T^2 + 2*T*m]
pofTClassic[T_, m_] := Sqrt[2*T*m]

thetamax[rB1_] := ArcSin[Sqrt[1/rB1]]

rG[p_, th2_, B_] := p*Sin[th2]/c/B

theta2[th0_, rB2_] := ArcSin[Sin[th0]*Sqrt[rB2]]

Chi2FitandPlot::usage = "Takes a drift distance data1, a variation of data2 (for 1 varying parameter b),
	 a list of values of varied parameter b and optionally plot options.
	 It creates chi2 values comparing all varied datasets with dataset1, fits a parable onto it and gives the results."
Chi2FitandPlot[Data1_, VariedData2_, Parameters_, plotopts : OptionsPattern[ListPlot]] := 
Module[
  {chi2list = Total[(Data1 - #[[All, 2]])^2] & /@ VariedData2,chi2data,fitfunc,fitresult},
  chi2data=Transpose[{Parameters, chi2list}];
  fitfunc[x_,x0_,a_,b_]:=a+b*(x-x0)^2;
  fitresult=NonlinearModelFit[
  	chi2data,
  	{fitfunc[b, b0, yoffset, scale], -0.002 < b0 < 0.002 && 0.00001 < scale < 0.001 && 0. < yoffset < 10^-8},

  	{{b0, -0.001}, {scale, 0.0005}, {yoffset, Min[chi2list]}},
  	b, Method -> "NMinimize"];
  {
  	{
  		{"b0 = ",fitresult["BestFitParameters"][[1,2]]},
  		fitresult,
  		fitresult["ParameterTable"]
  	}//TableForm,
  	Show[
  		ListPlot[chi2data, Sequence@plotopts],
  		Plot[fitresult[b],{b,-0.002,0.002},PlotStyle->Red]
  	]
  }
]