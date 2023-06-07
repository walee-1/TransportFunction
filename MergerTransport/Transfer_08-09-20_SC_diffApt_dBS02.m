(* Wolfram Language package *)
t0=AbsoluteTime[];
SetDirectory["/home/dmoser/eclipse-workspace/TransportProject/"];

Print["Configured Kernels: ",$ConfiguredKernels];
Print["Processor Count: ",$ProcessorCount];




BackScatterBoole = True;

a=-0.00075;
Bins = {129 , 128+32+16};

Get["MergerTransport/Merger2D_SimpleBS.m"];
Print["Modules loaded ..."];


(*things changed here*)

yAOffSC=0;

(*end of changes changed here*)
kernels=3;
LaunchKernels[kernels];
Print[kernels," Kernels launched ..."];


rFSC=2.036;
rASC=1;
rRxBSC=1;
rDSC=1;
bRxBSC=0.976;
g1SC=1.29;
g2SC=1.792;
alphaSC=180.03*Pi/180;

xAShiftSC=0;
yAShiftSC=2.873/1000;
yRxBShiftSC=1.518/1000;
xDShiftSC=0;
yDShiftSC=6.234/1000;
xAASC=2.5/1000;
yAASC=35/1000;
xAOffSC=6.5/1000;

{wnxSC, pnxSC, kx1SC, kx2SC, kx3SC, wnySC, pnySC, ky1SC, ky2SC, ky3SC}={1, 1, 1, 1., -1, 1, 1, 1, 1, -1}

xStart=-0.011;
xEnd=0.01;
yStart=-0.025;
yEnd=0.025;
xyBins=bin2DGen[xStart,xEnd,yStart,yEnd,64,1];
 
IntMethod = {"GlobalAdaptive", Method -> "MultidimensionalRule"};

Print["Run started at ",DateString[],", Bins: ",Bins];
  
BinwYShiftPrec44OriginalBinsb0NewBGrad = 
 ParallelMap[
  AbsoluteTiming[
    BinIntShift[
     a, #, 
     {alphaSC, bRxBSC, rFSC, rRxBSC, rASC, rDSC, g1SC, g2SC}, {xAShiftSC, yAShiftSC, yRxBShiftSC, xDShiftSC, yDShiftSC}, 
     {wnxSC, pnxSC, kx1SC, kx2SC, kx3SC, wnySC, pnySC, ky1SC, ky2SC, ky3SC}, {xAASC, yAASC, xAOffSC, yAOffSC},
     {IntMethod, 4, 2, 0, 2}, {IntMethod, 4, 2, 0, 2}, {IntMethod, 3, 6, 0, 1}
     ]] &, xyBins[[Bins[[1]] ;; Bins[[2]]]], 
  Method -> "FinestGrained"];

Export["MergerTransport/NC_opt_drA/b"<>ToString[IntegerPart[b*10000]]<>"/TransferResult_08-09-20_NC_opt_D_BSOn_drA_b"<>ToString[IntegerPart[b*10000]]<>"_"<>ToString[Bins[[1]]]<>"-"<>ToString[Bins[[2]]]<>".txt",BinwYShiftPrec44OriginalBinsb0NewBGrad,"Table"];
 
Print["b = ",IntegerPart[b*10000],"e-4"];
Print["Time consumption: ",(AbsoluteTime[]-t0)/3600.];
Print["MaxMemUse :",MaxMemoryUsed[]/1024/1024.," MB"];
CloseKernels[];
Quit[];


(* Wolfram Language package *)