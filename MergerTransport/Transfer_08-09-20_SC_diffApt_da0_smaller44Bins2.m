(* ::Package:: *)

(* Wolfram Language package *)
t0=AbsoluteTime[];
SetDirectory["/users/waleed.khalid/Mma/TransferProject/TransportFunction/"];

Print["Configured Kernels: ",$ConfiguredKernels];
Print["Processor Count: ",$ProcessorCount];

(*BackScatterBoole = True;*)

a=-0.00075;
Bins = {129 , 128+32+16};

Get["MergerTransport/Merger2DProton.m"];
Print["Modules loaded ..."];

kernels=6;
LaunchKernels[kernels];
Print[kernels," Kernels launched ..."];

DetHisto = 
  Flatten[{{{0.025`, 0.025568182000000002`, 0.026136364000000002`, 
      0.026704546000000003`, 0.027272728000000003`, 0.02784091`, 
      0.028409092`, 0.028977274`, 0.029545456`, 0.030113638`, 
      0.030681820000000002`, 0.031250002`, 0.031818184`, 0.032386366`,
       0.032954548`, 0.03352273`, 0.034090912`, 0.034659094`, 
      0.035227276`, 0.035795458`, 0.03636364`, 0.036931822`, 
      0.037500004`, 0.038068186000000004`, 0.038636368000000004`, 
      0.039204550000000005`, 0.039772732000000005`, 0.040340914`, 
      0.040909096`, 0.041477278`, 0.04204546`, 0.042613642`, 
      0.043181824`, 0.043750006`, 0.044318188`, 0.04488637`, 
      0.045454552`, 0.046022733999999996`, 0.046590915999999996`, 
      0.047159097999999997`, 0.04772728`, 0.048295462`, 0.048863644`, 
      0.049431826`}, {-0.02`, -0.01167`, -0.0033399999999999992`, 
      0.004990000000000001`, 0.013320000000000002`, 
      0.021650000000000006`, 0.029980000000000003`}}, {{0, 0, 0, 0, 0,
       0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 
      0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 
      0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 
      0}, {0, 14, 17, 13, 15, 0}, {6, 153, 149, 149, 126, 0}, {29, 
      505, 448, 459, 468, 1}, {75, 1060, 990, 1066, 1007, 1}, {149, 
      1872, 1820, 1766, 1710, 7}, {268, 2899, 2847, 2774, 2493, 
      19}, {413, 4015, 4040, 3853, 3535, 34}, {671, 5385, 5416, 5195, 
      4738, 68}, {876, 6841, 6841, 6728, 6099, 88}, {1205, 8336, 8352,
       8184, 7310, 149}, {1591, 9770, 9965, 9696, 8660, 202}, {1842, 
      11593, 11593, 11060, 10154, 280}, {2104, 12992, 13105, 12722, 
      11552, 358}, {2458, 14438, 14468, 14097, 12881, 478}, {2846, 
      15331, 15512, 15426, 13621, 528}, {2844, 15875, 16210, 16042, 
      14405, 611}, {2952, 15838, 16266, 16264, 14562, 669}, {2872, 
      15530, 15985, 15910, 14360, 709}, {2825, 14741, 15240, 15285, 
      13782, 678}, {2587, 13282, 14023, 14408, 12835, 644}, {2274, 
      11469, 12213, 12511, 11275, 637}, {1878, 9238, 10027, 10163, 
      9431, 474}, {1405, 6881, 7571, 7837, 7289, 432}, {1045, 4764, 
      5145, 5656, 5195, 341}, {730, 3332, 3694, 3937, 3521, 
      250}, {470, 2255, 2597, 2767, 2511, 182}, {346, 1486, 1699, 
      1692, 1693, 137}, {212, 974, 1062, 1125, 1024, 84}, {145, 524, 
      633, 689, 659, 33}, {60, 266, 357, 370, 343, 25}, {22, 128, 123,
       169, 197, 13}, {4, 42, 50, 60, 69, 1}, {1, 13, 14, 15, 14, 
      0}}}, 1];


xyBins= Flatten[Table[
    {
     {DetHisto[[1, xi]], DetHisto[[1, xi + 1]]},
     {DetHisto[[2, yi + 1]] , DetHisto[[2, yi]] }
     }, 
     {xi, 1, Length[DetHisto[[1]]] - 1}, 
     {yi, 1, Length[DetHisto[[2]]] - 1, 1}], 1];

rFSC=2.036;
rASC=0.937;
rRxBSC=0.902;
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
xAASC=1/1000;
yAASC=25/1000;
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

Export["MergerTransport/NC_opt_drF/b"<>ToString[IntegerPart[a*10000]]<>"/TransferResult_08-09-20_NC_opt_C_BSOn_drF_b"<>ToString[IntegerPart[a*10000]]<>"_"<>ToString[Bins[[1]]]<>"-"<>ToString[Bins[[2]]]<>".txt",BinwYShiftPrec44OriginalBinsb0NewBGrad,"Table"];
 
Print["a = ",IntegerPart[a*10000],"e-4"];
Print["Time consumption: ",(AbsoluteTime[]-t0)/3600.];
Print["MaxMemUse :",MaxMemoryUsed[]/1024/1024.," MB"];
CloseKernels[];
Quit[];

