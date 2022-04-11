(* Wolfram Language package *)
t0=AbsoluteTime[];
Print["Started on "<>DateString[]];
SetDirectory["/users/waleed.khalid/Mma/"];

Get["Efficiency/Efficiency.m"]

en = 4;

enTable={15,15.1,15.2,15.3,15.4,15.5,15.6,15.7,15.8};

resultDir="/users/waleed.khalid/plgad_results/";
Print["File will be: "<>resultDir<>"rangeHisto"<>ToString[enTable[[en]]]<>"AllAng.mx"];

Print["Energy Value "<>ToString[enTable[[en]]]];

$Pre=Function[Null,MemoryConstrained[#,32000000000-MemoryInUse[]],HoldAll];


If[FileExistsQ[resultDir<>"ImsilPlgadDataEn"<>ToString[enTable[[en]]]<>"AllAngIonList.mx"],
ImsilPlgadDataEnAllAngIonList=Import[resultDir<>"ImsilPlgadDataEn"<>ToString[enTable[[en]]]<>"AllAngIonList.mx"];,
Print["ERROR No File there, Please put the all ang list there"];Quit[]]




Export[resultDir<>"maxDepthHisto"<>ToString[enTable[[en]]]<>"AllAng.mx",
  Table[maxDepthModImsil[ImsilPlgadDataEnAllAngIonList[[ang]]],{ang, 
    Length[ImsilPlgadDataEnAllAngIonList]}]
    ];
    
    Export[resultDir<>"avgEHpairHisto"<>ToString[enTable[[en]]]<>"AllAngDefault.mx",
  Table[avgEhpairDepthModImsil[IImsilPlgadDataEnAllAngIonList[[ang]]],{ang, 
    Length[ImsilPlgadDataEnAllAngIonList]}]
    ];
    
Print["Time consumption: ",(AbsoluteTime[]-t0)/3600.];
Print["MaxMemUse :",MaxMemoryUsed[]/1024/1024/1024.," GB"];
CloseKernels[];
Quit[];