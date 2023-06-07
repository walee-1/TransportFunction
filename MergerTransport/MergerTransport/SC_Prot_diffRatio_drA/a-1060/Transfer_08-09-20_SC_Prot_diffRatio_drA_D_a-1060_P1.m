(* Wolfram Language package *)
t0=AbsoluteTime[];
SetDirectory["/home/waleed/Documents/Wolfram_Mathematica/TransportFunction"];

Print["Configured Kernels: ",$ConfiguredKernels];
Print["Processor Count: ",$ProcessorCount];

(*BackScatterBoole = True;*)

a=-0.106;
Bins = {50, 98};

Get["MergerTransport/Merger2DProton.m"];
Print["Modules loaded ..."];

kernels=4
LaunchKernels[kernels];
Print[kernels," Kernels launched ..."];


{{MeanX, MeanY}, {ApertX, ApertY}, {BDV}, {BF}, {BA}, {BRxB}, {BDet}}={{0.0301277, 1.00623}, {0.01, 
  0.035}, {0.220899}, {0.451107}, {0.219769}, {0.219722}, {0.207816}};
  
rF = BF/BDV;
rA = BA/BDV;
rRxB = BRxB/BDV;
rDet = BDet/BRxB;


rFSC=rF;
rASC=rA;
rRxBSC=rRxB;
bRxBSC=0.976;
g1SC=1.29;
g2SC=1.792;
alphaSC=180.03*Pi/180;
rDSC=rDet;
xAShiftSC=0;
yAShiftSC=2.873/1000;
yRxBShiftSC=1.518/1000;
xDShiftSC=0;
yDShiftSC=6.234/1000;
xAASC=10/1000;
yAASC=35/1000;
xAOffSC=30/1000;
yAOffSC=0;
wnxSC=10/100;
wnySC=6/100;
pnxSC=9/100;
pnySC=5/100;
kx1SC=0.9;
ky1SC=0.9;
kx2SC=0;
ky2SC=0;
kx3SC=-0.9;
ky3SC=-0.9;

xStart=0.02;
xEnd=0.065;
yStart=-0.03;
yEnd=0.04;
xyBins=bin2DGen[xStart,xEnd,yStart,yEnd,24,11];
 
IntMethod = {"GlobalAdaptive", Method -> "MultidimensionalRule"};

Print["Run started at ",DateString[],", Bins: ",Bins];
  
BinwYShiftPrec44OriginalBinsb0NewBGrad = 
 ParallelMap[
  AbsoluteTiming[
    BinIntShift[
     a, #, 
     {alphaSC, bRxBSC, rFSC, rRxBSC, rASC+rASC*10^-4, rDSC, g1SC, g2SC}, {xAShiftSC, yAShiftSC, yRxBShiftSC, xDShiftSC, yDShiftSC}, 
     {wnxSC, pnxSC, kx1SC, kx2SC, kx3SC, wnySC, pnySC, ky1SC, ky2SC, ky3SC}, {xAASC, yAASC, xAOffSC, yAOffSC},
     {IntMethod, 4, 2, 0, 2}, {IntMethod, 4, 2, 0, 2}, {IntMethod, 3, 6, 0, 1}
     ]] &, xyBins[[Bins[[1]] ;; Bins[[2]]]], 
  Method -> "FinestGrained"];

Export["MergerTransport/SC_Prot_diffRatio_drA/a"<>ToString[IntegerPart[a*10000]]<>"/TransferResult_08-09-20_SC_diffRatio_D_drA_a"<>ToString[IntegerPart[a*10000]]<>"_"<>ToString[Bins[[1]]]<>"-"<>ToString[Bins[[2]]]<>".txt",BinwYShiftPrec44OriginalBinsb0NewBGrad,"Table"];
 
Print["a = ",IntegerPart[a*10000],"e-4"];
Print["Time consumption: ",(AbsoluteTime[]-t0)/3600.];
Print["MaxMemUse :",MaxMemoryUsed[]/1024/1024.," MB"];
CloseKernels[];
Quit[];

