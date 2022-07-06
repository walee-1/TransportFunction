(* Wolfram Language package *)
t0=AbsoluteTime[];
SetDirectory["/users/waleed.khalid/Mma/TransferProject/TransportFunction/"];

Print["Configured Kernels: ",$ConfiguredKernels];
Print["Processor Count: ",$ProcessorCount];

(*BackScatterBoole = True;*)

a=-0.00075;
Bins = {129 , 128+32+16};

Get["MergerTransport/Merger2DProtonShifted.m"];
Print["Modules loaded ..."];

kernels=6;
LaunchKernels[kernels];
Print[kernels," Kernels launched ..."];

rFSC=2.04214;
rASC=0.994885;
rRxBSC=0.994672;
bRxBSC=0.976;
g1SC=1.29;
g2SC=1.792;
alphaSC=180.03*Pi/180;
rDSC=0.88;
xAShiftSC=0;
yAShiftSC=2.873/1000;
yRxBShiftSC=1.518/1000;
xDShiftSC=0;
yDShiftSC=6.234/1000;
xAASC=3.5/1000;
yAASC=35/1000;
xAOffSC=-10/1000;
yAOffSC=0;
{wnxSC, pnxSC, kx1SC, kx2SC, kx3SC, wnySC, pnySC, ky1SC, ky2SC, ky3SC}={0.001, 0.0009, 0.9, 0., -0.9, 0.001, 0.0009, 0.9, 0., -0.9}

xStart=-0.012;
xEnd=0.012;
yStart=-0.022;
yEnd=0.025;
xyBins=bin2DGen[xStart,xEnd,yStart,yEnd,64,1];
 
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

Export["MergerTransport/NC_opt_drF/b"<>ToString[IntegerPart[a*10000]]<>"/TransferResult_08-09-20_NC_opt_C_BSOn_drF_b"<>ToString[IntegerPart[a*10000]]<>"_"<>ToString[Bins[[1]]]<>"-"<>ToString[Bins[[2]]]<>".txt",BinwYShiftPrec44OriginalBinsb0NewBGrad,"Table"];
 
Print["a = ",IntegerPart[a*10000],"e-4"];
Print["Time consumption: ",(AbsoluteTime[]-t0)/3600.];
Print["MaxMemUse :",MaxMemoryUsed[]/1024/1024.," MB"];
CloseKernels[];
Quit[];
