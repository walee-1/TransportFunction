(* Wolfram Language package *)


(*first, a function that selects the correct integrand dependent on the optionvalue given*)
Get["Integrands.m"];


IntegrandSelect[string_]=Which[
		string=="DSimpleApertMC",
		IntegrandDSimple[xi,x0,{b,...}]
]



(* we need to give integrand (spectrum, drift (bfield), aperture, nbeam etc.) either directly or as option, integration coords + limits, *)
(* x fix or integrated should be a flag*)

Options[TransportIntegration]={Integrand->"DSimpleApertMC",xBinInt->True,PrecisionGoal->3,IntMethod->Automatic};

TransportIntegration[XBinList_List,IntVarList_List,OptionsPattern[TransportIntegration]]:=
Module[{
		integrand=	IntegrandSelect[OptionValue[Integrand]]},
	{
		ParallelTable[
			NIntegrate[
				integrand[y,xi]/.(*If[OptionValue[xBinInt],xi->xDet,*)xi->XBinList[[xbin]](*]*),
				(*the following joins the other integrations with the xbin int and makes a sequence out of it, if flag is true, otherwise without xInt*)
				Sequence@@(*If[
					OptionValue[xBinInt],
					Join[IntVarList,{xDet,XBinList[[xbin]],XBinList[[xbin+1]]}],*)
					IntVarList
					(*]*),
				PrecisionGoal->OptionValue[PrecisionGoal],Method->OptionValue[IntMethod]
			],
			{xbin,1,(*If[OptionValue[xBinInt],Length[XBinList]-1,*)Length[XBinList](*]*)}
		],
		Print[OptionValue[PrecisionGoal],OptionValue[xBinInt]]
	}
]