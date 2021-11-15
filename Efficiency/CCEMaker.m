(* Wolfram Language package *)
t0=AbsoluteTime[];
SetDirectory["/users/waleed.khalid/Mma/"];

Get["Efficiency/Efficiency.m"]

kernels=1;

LaunchKernels[kernels];

Print["Kernels Launched"]

$Pre=Function[Null,MemoryConstrained[#,32000000000-MemoryInUse[]],HoldAll];


ImsilPlgadDataEnFineAllEnAllAngIonList=Table[If[en==0,Import["/users/waleed.khalid/plgad_results/ImsilPlgadDataEnFine15EnAllAngIonList.mx"],Import["/users/waleed.khalid/plgad_results/ImsilPlgadDataEnFine15p"<>ToString[en-1]<>"EnAllAngIonList.mx"]],{en,9}];

Export["/users/waleed.khalid/plgad_results/ImsilPlgadDataEnFineAllEnAllAngIonList.mx",ImsilPlgadDataEnFineAllEnAllAngIonList];

Export["/users/waleed.khalid/plgad_results/cceAdjustIonListEn15FineAllEnAllAng.mx",
  Table[cceModImsil[
    ImsilPlgadDataEnFineAllEnAllAngIonList[[en, ang]][[i]], 
    gammaVals[[1]], tauVals[[1]]], {en, 1}, {ang, 
    Dimensions[ImsilPlgadDataEnFineAllEnAllAngIonList][[2]]}, {i, 
    Length[ImsilPlgadDataEnFineAllEnAllAngIonList[[en, ang]]]}]];

Print["Time consumption: ",(AbsoluteTime[]-t0)/3600.];
Print["MaxMemUse :",MaxMemoryUsed[]/1024/1024/1024.," GB"];
CloseKernels[];
Quit[];