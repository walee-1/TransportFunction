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

Chi2FitandPlot[Data1_, VariedData2_, Parameters_,precision_, plotopts : OptionsPattern[ListPlot]] := 
Module[
  {chi2list = Total[(Data1 - #[[All, 2]])^2] & /@ VariedData2,deltachi2Index,deltachi2list,chi2data,fitfunc,fitresult,errorplotdata},
  
  chi2data=Transpose[{Parameters, chi2list}];
  
  deltachi2Index=Table[
  	2*Abs[Data1[[i]]-VariedData2[[b,i,2]]]*10^-precision*Sqrt[Data1[[i]]^2+VariedData2[[b,i,2]]^2],
  		{b,1,Length[Parameters]},{i,1,Length[Data1]}
  		];
  deltachi2list=Table[
  	Sqrt[
  		Sum[deltachi2Index[[b,i]]^2,{i,1,Length[Data1]}]
  		],
  		{b,1,Length[Parameters]}];
  errorplotdata=
  	Transpose[
  		{
  			Parameters,
  			Table[Around[chi2list[[b]],deltachi2list[[b]]],{b,1,Length[Parameters]}]
  		}
  	];
  
  fitfunc[x_,x0_,a_,b_]:=a+b*(x-x0)^2;
  fitresult=NonlinearModelFit[
  	chi2data,
  	{fitfunc[b, b0, yoffset, scale], -0.005 < b0 < 0.005 && 0.00001 < scale < 0.01 && 0. < yoffset < Max[chi2list]},

  	{{b0, -0.001}, {scale, 0.0005}, {yoffset, Min[chi2list]}},
  	b, Method -> "NMinimize",Weights->1/deltachi2list^2,VarianceEstimatorFunction->(1&)
  	];
  {
  	{
  		{"b0 = ",fitresult["BestFitParameters"][[1,2]]},
  		fitresult,
  		fitresult["ParameterTable"],
  		fitresult["ANOVATable"]
  	}//TableForm,
  	Show[
  		ListPlot[errorplotdata, Sequence@plotopts],
  		Plot[fitresult[b],{b,-0.003,0.003},PlotStyle->Red]
  	]
  }
]