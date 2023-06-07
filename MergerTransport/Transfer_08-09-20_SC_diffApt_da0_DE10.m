(* Wolfram Language package *)
t0=AbsoluteTime[];
SetDirectory["/home/dmoser/eclipse-workspace/TransportProject/"];

Print["Configured Kernels: ",$ConfiguredKernels];
Print["Processor Count: ",$ProcessorCount];

BackScatterBoole = True;

a=-0.00075;
Bins = {129 , 128+32+16};

Get["MergerTransport/Merger2D_Simple.m"];
Print["Modules loaded ..."];

kernels=3;
LaunchKernels[kernels];
Print[kernels," Kernels launched ..."];

{{MeanX, MeanY}, {ApertX, ApertY}, {BDV}, {BF}, {BA}, {BRxB}, {BDet}}={{0.0301277, 1.00623}, {0.01, 
  0.035}, {0.220899}, {0.451107}, {0.219769}, {0.219722}, {0.207816}};
  
ApertXOff = 0.03; 
ApertYOff = 0.;

(*CenterfitG1=0.96907;
CenterfitG2=0.973888;*)
BlinefitG1=0.965547;
BlinefitG2=0.974009;

rF = BF/BDV;
rA = BA/BDV;
rRxB = BRxB/BDV;
rDet = BDet/BRxB;

(*AlphaFromBlines = 178.20869281586448`/180*Pi;*)

AlphaFromBlinesCorr = 180.11/180.*Pi;

(*BfromCenter=0.221032;*)

BfromCentralBline=0.219679;

XAShift = 0.0;YAShift = 0.007188045691239971`;YRxBShift = 0.00677284727516161`;YDetShift = 0.015928000000000164`;
XRxBShift = 0.;XDetShift = 0.;

rFSC=2.036;
rASC=0.937;
rRxBSC=0.902;
rDSC=0.880;
bRxBSC=0.976;
g1SC=1.29;
g2SC=1.792;
alphaSC=180.03*Pi/180;

xAShiftSC=0;
yAShiftSC=2.873/1000;
yRxBShiftSC=1.518/1000;
xDShiftSC=0;
yDShiftSC=6.234/1000;
xAASC=10/1000;
yAASC=35/1000;
xAOffSC=6.5/1000;
yAOffSC=0;
{wnxSC, pnxSC, kx1SC, kx2SC, kx3SC, wnySC, pnySC, ky1SC, ky2SC, ky3SC}={0.001, 0.0009, 0.9, 0., -0.9, 0.001, 0.0009, 0.9, 0., -0.9}

xStart=-0.013;
xEnd=0.013;
yStart=-0.041;
yEnd=0.041;
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


