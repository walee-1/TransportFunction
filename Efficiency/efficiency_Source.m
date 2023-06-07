(* Wolfram Language package *)(* Wolfram Language package *)
t0=AbsoluteTime[];
Print["Started on "<>DateString[]];
SetDirectory["/users/waleed.khalid/Mma/"];

Get["Efficiency/Efficiency.m"]

entoFind = 4;

enTable = {5,10,15,15.1,15.2,15.3,15.4,15.5,15.6,15.7,15.8,17,19,20,25,30};

binsToJoin = 12;

pos=Position[enTable,entoFind];
en=If[pos=={},Print["Energy not Found, expand energy table in CCEMaker_Source.m"];Quit[],pos[[1,1]]];
angTab = {5};

folderName = "Al2O3_4nm_Efficiency";
resultDir="/users/waleed.khalid/plgad_results/";
Print["File will be: "<>resultDir<>folderName<>"_efficiency"<>ToString[enTable[[en]]]<>"Ang"<>ToString[angTab[[1]]]<>"_"<>ToString[angTab[[-1]]]<>".mx"];



$Pre=Function[Null,MemoryConstrained[#,32000000000-MemoryInUse[]],HoldAll];



If[!FileExistsQ[resultDir<>folderName<>"_cceAdjustIonListEn"<>ToString[enTable[[en]]]<>"Ang"<>ToString[angTab[[1]]]<>"_"<>ToString[angTab[[-1]]]<>".mx"],
Print["ERROR No File there, Please put the cce list there"];Quit[]]

If[!FileExistsQ[resultDir<>folderName<>"_cceAdjustIonListEn"<>ToString[enTable[[en]]]<>"Ang"<>ToString[angTab[[1]]]<>"_"<>ToString[angTab[[-1]]]<>"Worse.mx"],
Print["ERROR No Worse File there, Please put the cce list there"];Quit[]]

ImsilPlgadDataEnAllAngIonList=Import[resultDir<>folderName<>"_cceAdjustIonListEn"<>ToString[enTable[[en]]]<>"Ang"<>ToString[angTab[[1]]]<>"_"<>ToString[angTab[[-1]]]<>".mx"];


Export[resultDir <> folderName <> "_avgEHpairHisto" <> 
    ToString[enTable[[en]]]<> "Ang" <> ToString[angTab[[1]]] <> "_" <> 
    ToString[angTab[[-1]]] <> ".mx", 
  Table[avgEhpairDepthModImsil[
    ImsilPlgadDataEnAllAngIonList[[ang]]], {ang, 
    Length[ImsilPlgadDataEnAllAngIonList]}]];
    
Export[resultDir <> folderName <> "_maxDepthHisto" <> 
   ToString[enTable[[en]]] <> "Ang" <> ToString[angTab[[1]]] <> "_" <>
    ToString[angTab[[-1]]] <> ".mx", 
  Table[maxDepthModImsil[ImsilPlgadDataEnAllAngIonList[[ang]]], {ang, 
    Length[ImsilPlgadDataEnAllAngIonList]}]];



   Export[resultDir <> folderName <> "_HistoData"<>ToString[enTable[[en]]]<>"Ang"<> ToString[angTab[[1]]]<> "_" <> 
    ToString[angTab[[-1]]] <> ".mx", 
  histoDataAll = 
   Table[histogramModImsilErr[ImsilPlgadDataEnAllAngIonList[[ang]], 
     binsToJoin], {ang, Length[ImsilPlgadDataEnAllAngIonList]}]];
     
Export[resultDir <> folderName <> "_EfficiencyData" <> 
    ToString[enTable[[en]]]<> "Ang" <> ToString[angTab[[1]]] <> "_" <> 
    ToString[angTab[[-1]]] <> ".mx", 
  Table[effiencyModImsil[histoDataAll[[ang]][[2]], 
    Length[ImsilPlgadDataEnAllAngIonList]], {ang, 
    Length[ImsilPlgadDataEnAllAngIonList]}]];
    
ClearAll[ImsilPlgadDataEnAllAngIonList, 
    histoDataAll]
    
    
     ImsilPlgadDataEnAllAngIonList = 
  Import[resultDir <> folderName <> "_cceAdjustIonListEn" <> 
    ToString[enTable[[en]]] <> "Ang" <> ToString[angTab[[1]]] <> "_" <>
     ToString[angTab[[-1]]] <> "Worse.mx"];
     
Export[resultDir <> folderName <> "_avgEHpairHisto" <> 
    ToString[enTable[[en]]]<> "Ang" <> ToString[angTab[[1]]] <> "_" <> 
    ToString[angTab[[-1]]] <> "Worse.mx", 
  Table[avgEhpairDepthModImsil[
    ImsilPlgadDataEnAllAngIonList[[ang]]], {ang, 
    Length[ImsilPlgadDataEnAllAngIonList]}]];
    
Export[resultDir <> folderName <> "_maxDepthHisto" <> 
   ToString[enTable[[en]]] <> "Ang" <> ToString[angTab[[1]]] <> "_" <>
    ToString[angTab[[-1]]] <> "Worse.mx", 
  Table[maxDepthModImsil[ImsilPlgadDataEnAllAngIonList[[ang]]], {ang, 
    Length[ImsilPlgadDataEnAllAngIonList]}]];
    
Export[resultDir <> folderName <> "_HistoData" <> 
    ToString[enTable[[en]]]<> "Ang" <> ToString[angTab[[1]]] <> "_" <> 
    ToString[angTab[[-1]]] <> "Worse.mx", 
  histoDataAll = 
   Table[histogramModImsilErr[ImsilPlgadDataEnAllAngIonList[[ang]], 
     binsToJoin], {ang, Length[ImsilPlgadDataEnAllAngIonList]}]];
     
Export[resultDir <> folderName <> "_EfficiencyData" <> 
    ToString[enTable[[en]]]<> "Ang" <> ToString[angTab[[1]]] <> "_" <> 
    ToString[angTab[[-1]]] <> "Worse.mx", 
  Table[effiencyModImsil[histoDataAll[[ang]][[2]], 
    Length[ImsilPlgadDataEnAllAngIonList]], {ang, 
    Length[ImsilPlgadDataEnAllAngIonList]}]];



    
Print["Time consumption: ",(AbsoluteTime[]-t0)/3600.];
Print["MaxMemUse :",MaxMemoryUsed[]/1024/1024/1024.," GB"];
CloseKernels[];
Quit[];