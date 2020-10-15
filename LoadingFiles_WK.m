(* Wolfram Language package *)


Get["Efficiency/Efficiency.m"]
Get["Common/physMathForms.m"]
Get["Common/CommonFunctions.m"]
Get["General/easeFunctions.m"]


additionalDependencies::usage="Pure option function to choose if to load extra dependency files.
	Usage:
		additionalDependencies[option Value->Bool]
		By default all values are set to false
	Options:
		LoadAll = Loads all the files in the function
		LoadSpectras = Loads the Nachtmann and the Electron Spectra
		Nachtmann = Loads the Nachtmann Spectrum
		Glueck93 = Loads the Glueck 93 Proton Spectrum
		Glueck = Loads the Glueck 95 Proton Spectrum
		Electron = Loads only the Electron Spectrum
		MyColors = Loads my color shceme
"
Options[additionalDependencies]={Nachtmann->False,Glueck93->False,Glueck->False,Electron->False,LoadAll->False,MyColors->False,LoadSpectras->False}

additionalDependencies[opts:OptionsPattern[]]:=Module[{},
	If[OptionValue[LoadAll],OptionValue[Nachtmann]=True;OptionValue[Glueck93]=True;OptionValue[Glueck]=True;OptionValue[Electron]=True;OptionValue[MyColors]=True];
	If[OptionValue[LoadSpectras],OptionValue[Nachtmann]=True;OptionValue[Electron]=True];
	If[OptionValue[Nachtmann],Get["Spectra/ProtonSpectrumNachtmann.m"]];
	If[OptionValue[Glueck93],Get["Spectra/Glueck93ProtonSpectrum.m"]];
	If[OptionValue[Glueck],Get["Spectra/GlueckProtonSpectrum.m"]];
	If[OptionValue[Electron],Get["Spectra/ElectronSpectrum.m"]];
	If[OptionValue[MyColors],Get["General/easeFunctions.m"]];
];

