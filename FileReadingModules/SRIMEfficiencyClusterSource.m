t0=AbsoluteTime[];
Print["Started on "<>DateString[]];
SetDirectory["/users/waleed.khalid/Mma/"]
Get["FileReadingModules/FileReadingModules.m"];

$Pre=Function[Null,MemoryConstrained[#,32000000000-MemoryInUse[]],HoldAll];

en = 15;

ang = 0;

window = 40;

geom = "SiO2"

filename = "/users/waleed.khalid/SiO2_Track/Energy_"<>ToString[en]<>"/Angle_"<>ToString[ang]<>"/"

runNo = StringSplit[FileNames[filename <> "*"], "/"][[All, -1]];


exyzFile=Table[exyzReadingMod2[filename <> runNo[[i]] <> "/EXYZ.txt"], {i, 
   Length[runNo]}];
   
exyzFile = Flatten[DeleteDuplicates[exyzFile], 1];

lenOfFile=Length[exyzFile];

effList = Table[mainMod[exyzFile,i,window],{i,lenOfFile}];

ClearAll[exyzFile];

Export["/users/waleed.khalid/SRIM_Results/"<>geom<>"_"<>ToString[window]<>"effListEn"<>ToString[en]<>"Ang"<>ToString[ang]<>".mx",effList]

Print["File written: "<>geom<>"_"<>ToString[window]<>"effListEn"<>ToString[en]<>"Ang"<>ToString[ang]<>".mx"]

cceFileBest = Table[cceWeightingModulePerIon[effList[[i]],100,gammaVals[[1]],tauVals[[1]]], {i, lenOfFile}];

cceFileWorst = Table[cceWeightingModulePerIon[effList[[i]],100,gammaVals[[2]],tauVals[[2]]], {i, lenOfFile}];

Export["/users/waleed.khalid/SRIM_Results/"<>geom<>"_"<>ToString[window]<>"cceListEn"<>ToString[en]<>"Ang"<>ToString[ang]<>"Best.mx",cceFileBest]

Print["File written: "<>geom<>"_"<>ToString[window]<>"cceListEn"<>ToString[en]<>"Ang"<>ToString[ang]<>"Best.mx"]

Export["/users/waleed.khalid/SRIM_Results/"<>geom<>"_"<>ToString[window]<>"cceListEn"<>ToString[en]<>"Ang"<>ToString[ang]<>"Worse.mx",cceFileWorst]

Print["File written: "<>geom<>"_"<>ToString[window]<>"cceListEn"<>ToString[en]<>"Ang"<>ToString[ang]<>"Worse.mx"]

ClearAll[effList];

histoListBest=histogramMod2[cceFileBest,10];

histoListWorst=histogramMod2[cceFileWorst,10];

ClearAll[cceFileBest,cceFileWorst];

Export["/users/waleed.khalid/SRIM_Results/"<>geom<>"_"<>ToString[window]<>"histoListEn"<>ToString[en]<>"Ang"<>ToString[ang]<>"Best.mx",histoListBest]

Export["/users/waleed.khalid/SRIM_Results/"<>geom<>"_"<>ToString[window]<>"histoListEn"<>ToString[en]<>"Ang"<>ToString[ang]<>"Worse.mx",histoListWorst]

effListBest=efficiencyMod[histoListBest,lenOfFile];
effListWorst=efficiencyMod[histoListWorst,lenOfFile];



Export["/users/waleed.khalid/SRIM_Results/"<>geom<>"_"<>ToString[window]<>"FinalEffListEn"<>ToString[en]<>"Ang"<>ToString[ang]<>"Best.mx",effListBest]

Export["/users/waleed.khalid/SRIM_Results/"<>geom<>"_"<>ToString[window]<>"FinalEffListEn"<>ToString[en]<>"Ang"<>ToString[ang]<>"Worse.mx",effListWorst]


Print["Time consumption: ",(AbsoluteTime[]-t0)/3600.];
Print["MaxMemUse :",MaxMemoryUsed[]/1024/1024/1024.," GB"];
CloseKernels[];
Quit[];


