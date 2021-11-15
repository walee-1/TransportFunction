(* Wolfram Language package *)

Get["Common/Constants.m"];

Eofp[p_] := Sqrt[me^2 + p^2]

Tofpe[p_] := Eofp[p] - me

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


TransferToMCParableFit[MCData_, TransferDataList_, ParaList_] := 
 Module[
  {N = Total[MCData], NonZeroPos, MCDataNormed, SigmaNormed, TransferDataListNormed, Chi2List, dChi2List, errorplotdata, fitresult, fitfunc},
  NonZeroPos=Flatten[Position[MCData,_?(#!=0&)],1];
  MCDataNormed = MCData[[NonZeroPos]]/N;
  SigmaNormed = Sqrt[MCData[[NonZeroPos]]]/N;
  TransferDataListNormed = ((#/Total[#])[[NonZeroPos]]) &/@ TransferDataList;
  Chi2List = Table[Total[(MCDataNormed - TransferDataListNormed[[b]])^2/SigmaNormed^2]/(Length[NonZeroPos-1]), {b, 1, Length[ParaList]}];
  dChi2List = Table[2*Sqrt[Total[(MCDataNormed - TransferDataListNormed[[b]])^2/SigmaNormed^2]]/(Length[NonZeroPos-1]), {b, 1, Length[ParaList]}];
  errorplotdata = Transpose[{ParaList, Table[Around[Chi2List[[b]], dChi2List[[b]]], {b, 1, Length[ParaList]}]}];
  fitfunc[x_, x0_, f_, g_] := f + g*(x - x0)^2; 
  fitresult = 
   NonlinearModelFit[
    Transpose[{ParaList, Chi2List}],
    {fitfunc[b, b0, yoffset, scale], -0.01 < b0 < 0.01 && 0.001 < scale < 1 && 0. < yoffset < Max[Chi2List]}, 
    {{b0, -0.001}, {scale, 0.0005}, {yoffset, Min[Chi2List]}}, b,
    Method -> "NMinimize"(*, Weights -> 1/dChi2List^2, VarianceEstimatorFunction -> (1 &)*)];
  {
  	{Chi2List,dChi2List},
   {{"b0 = ", fitresult["BestFitParameters"][[1, 2]]}, fitresult, 
     fitresult["ParameterTable"], fitresult["ANOVATable"]} // 
    TableForm,
   Show[
     ListPlot[errorplotdata,PlotRange->{{-0.005,0.005},Automatic}],
     Plot[fitresult[b], {b, -0.003, 0.003}, PlotStyle -> Red]
    ]
   }
  ]
  
aOfLambda[lambda_] := (1 - Abs[lambda]^2)/(1 + 3 Abs[lambda]^2)

lambdaOfa[a_] := -(Sqrt[1 - a])/(Sqrt[1 + 3 a]) 

(* ::Subsubsection:: *)
(* Region Title *)

(* ::Subsubsection:: *)
(* Region Title *)
TransferToRefParableFit[RefTransfer_, TransferDataList_, ParaList_, Prec_] := 
 Module[
  {RefNorm = Total[RefTransfer], (*NonZeroPos,*) RefNormed, Sigma, TransferDataListNormed, Chi2List, dChi2List, errorplotdata, fitresult, fitfunc},
 
  RefNormed = RefTransfer/RefNorm;
  Sigma = If[# ==0, 10000, 10^-Prec*#]&/@RefNormed;
  TransferDataListNormed = (#/Total[#]) &/@ TransferDataList;
  Chi2List = Table[
  	Total[(RefNormed - TransferDataListNormed[[b]])^2], 
  		{b, 1, Length[ParaList]}];
  dChi2List = Table[
  	10^-Prec*Sqrt[
  		Total[(RefNormed - TransferDataListNormed[[b]])^2*(RefNormed^2 + TransferDataListNormed[[b]]^2)]
  		], 
  		{b, 1, Length[ParaList]}];
  errorplotdata = Transpose[{
  	ParaList, 
  	Table[10^9*Around[Chi2List[[b]], dChi2List[[b]]], {b, 1, Length[ParaList]}]
  }];
  fitfunc[x_, x0_, f_, g_] := f + g*(x - x0)^2; 
  fitresult = 
   NonlinearModelFit[
    Transpose[{ParaList, Chi2List}],
    {fitfunc[b, b0, yoffset, scale], -0.001 < b0 < 0.005 && 0. < scale < 1*10^-1 && 0. < yoffset <= Min[Chi2List]}, 
    {{b0, 0.0005}, {scale, 5*10^-4}, {yoffset, Min[Chi2List]}}, b,
    Method -> "NMinimize", Weights -> 1/dChi2List^2, VarianceEstimatorFunction -> (1 &),PrecisionGoal->10];
  {
  	{Chi2List,dChi2List},
   {{"b0 = ", fitresult["BestFitParameters"][[1, 2]]}, fitresult, 
     fitresult["ParameterTable"], fitresult["ANOVATable"]} // 
    TableForm,
   Show[
     ListPlot[errorplotdata,PlotRange->{{-0.002,0.002},Automatic}],
     Plot[10^9*fitresult[b], {b, -0.003, 0.003}, PlotStyle -> Red]
    ]
   }
  ]
  
 TransferToRefParableFitForProtons[RefTransfer_, TransferDataList_, ParaList_, Prec_] := 
 Module[
  {RefNorm = Total[RefTransfer], (*NonZeroPos,*) RefNormed, Sigma, TransferDataListNormed, Chi2List, dChi2List, errorplotdata, fitresult, fitfunc,bMin,bMax},
  bMin=Min[ParaList]-0.0005;
  bMax=Max[ParaList]+0.0005;
  RefNormed = RefTransfer/RefNorm;
  Sigma = If[# ==0, 10000, 10^-Prec*#]&/@RefNormed;
  TransferDataListNormed = (#/Total[#]) &/@ TransferDataList;
  Chi2List = Table[
  	Total[(RefNormed - TransferDataListNormed[[b]])^2], 
  		{b, 1, Length[ParaList]}];
  dChi2List = Table[
  	10^-Prec*Sqrt[
  		Total[(RefNormed - TransferDataListNormed[[b]])^2*(RefNormed^2 + TransferDataListNormed[[b]]^2)]
  		], 
  		{b, 1, Length[ParaList]}];
  errorplotdata = Transpose[{
  	ParaList, 
  	Table[10^9*Around[Chi2List[[b]], dChi2List[[b]]], {b, 1, Length[ParaList]}]
  }];
  fitfunc[x_, x0_, f_, g_] := f + g*(x - x0)^2; 
  fitresult = 
   NonlinearModelFit[
    Transpose[{ParaList, Chi2List}],
    {fitfunc[b, b0, yoffset, scale], bMin < b0 < bMax && 0. < scale < 1*10^-1 && 0. < yoffset <= Min[Chi2List]}, 
    {{b0, -0.1051}, {scale, 5*10^-4}, {yoffset, Min[Chi2List]}}, b,
    Method -> "NMinimize", Weights -> 1/dChi2List^2, VarianceEstimatorFunction -> (1 &),PrecisionGoal->10];
  {
  	{Chi2List,dChi2List},
   {{"a0 = ", fitresult["BestFitParameters"][[1, 2]]}, fitresult, 
     fitresult["ParameterTable"], fitresult["ANOVATable"]} // 
    TableForm,
   Show[
     ListPlot[errorplotdata,PlotRange->{{bMin,bMax},Automatic}],
     Plot[10^9*fitresult[b], {b, bMin, bMax}, PlotStyle -> Red]
    ]
   }
  ]
  
 

  
 TransferToRefParableFitForProtonsUnlimited[RefTransfer_, TransferDataList_, ParaList_, Prec_] := 
 Module[
  {RefNorm = Total[RefTransfer], (*NonZeroPos,*) RefNormed, Sigma, TransferDataListNormed, Chi2List, dChi2List, errorplotdata, fitresult, fitfunc,bMin,bMax},
  bMin=-0.2;
  bMax=0;
  RefNormed = RefTransfer/RefNorm;
  Sigma = If[# ==0, 10000, 10^-Prec*#]&/@RefNormed;
  TransferDataListNormed = (#/Total[#]) &/@ TransferDataList;
  Chi2List = Table[
  	Total[(RefNormed - TransferDataListNormed[[b]])^2], 
  		{b, 1, Length[ParaList]}];
  dChi2List = Table[
  	10^-Prec*Sqrt[
  		Total[(RefNormed - TransferDataListNormed[[b]])^2*(RefNormed^2 + TransferDataListNormed[[b]]^2)]
  		], 
  		{b, 1, Length[ParaList]}];
  errorplotdata = Transpose[{
  	ParaList, 
  	Table[10^9*Around[Chi2List[[b]], dChi2List[[b]]], {b, 1, Length[ParaList]}]
  }];
  fitfunc[x_, x0_, f_, g_] := f + g*(x - x0)^2; 
  fitresult = 
   NonlinearModelFit[
    Transpose[{ParaList, Chi2List}],
    {fitfunc[b, b0, yoffset, scale], bMin < b0 < bMax && 0. < scale < 1*10^-1 && 0. < yoffset <= Min[Chi2List]}, 
    {{b0, -0.1051}, {scale, 5*10^-4}, {yoffset, Min[Chi2List]}}, b,
    Method -> "NMinimize", Weights -> 1/dChi2List^2, VarianceEstimatorFunction -> (1 &),PrecisionGoal->10];
  {
  	{Chi2List,dChi2List},
   {{"a0 = ", fitresult["BestFitParameters"][[1, 2]]}, fitresult, 
     fitresult["ParameterTable"], fitresult["ANOVATable"]} // 
    TableForm,
   Show[
     ListPlot[errorplotdata,PlotRange->{{Min[ParaList]-0.0005,bMax=Max[ParaList]+0.0005},Automatic}],
     Plot[10^9*fitresult[b], {b, bMin, bMax}, PlotStyle -> Red]
    ]
   }
  ]
  