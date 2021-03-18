(* ::Package:: *)

(* Wolfram Language package *)
t0=AbsoluteTime[];
SetDirectory["C:/Users/Synchrotronus/Desktop/Transport_24-09-20"];

Print["Configured Kernels: ",$ConfiguredKernels];
Print["Processor Count: ",$ProcessorCount];

BackScatterBoole = True;

b=0.0015;
Bins = {123, 138};

Get["MergerTransport/Merger2D.m"];
Print["Modules loaded ..."];

kernels=2;
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

R = 0.999002;

DetHisto= {{-0.07, -0.065, -0.06, -0.055, -0.05, -0.045, -0.04, -0.035, -0.03, \
-0.025, -0.02, -0.015, -0.01, -0.005, 1.38778*10^-17, 0.005, 0.01, 
  0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 
  0.045}, {-1.054, -1.049, -1.044, -1.039, -1.034, -1.029, -1.024, \
-1.019, -1.014, -1.009, -1.004, -0.999002, -0.994002, -0.989002, \
-0.984002, -0.979002, -0.974002, -0.969002, -0.964002, -0.959002, \
-0.954002, -0.949002, -0.944002}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
   0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 4, 5, 5, 8, 10, 8, 8, 
  6, 7, 1, 2, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 2, 10, 21, 53, 53, 54, 
  70, 54, 58, 46, 21, 8, 3, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 25, 53, 139, 
  172, 225, 286, 295, 273, 228, 180, 97, 55, 19, 2, 0, 0, 0, 0, 
  0}, {0, 0, 17, 66, 190, 414, 617, 832, 888, 901, 832, 724, 556, 375,
   198, 49, 14, 1, 0, 0, 0, 0}, {0, 1, 34, 175, 526, 986, 1611, 2023, 
  2292, 2277, 2214, 1873, 1477, 947, 504, 192, 41, 1, 0, 0, 0, 0}, {0,
   3, 88, 369, 1087, 2197, 3813, 4904, 5250, 5381, 5155, 4475, 3545, 
  2138, 1186, 436, 106, 11, 0, 0, 0, 0}, {0, 1, 110, 696, 1886, 4112, 
  7968, 10086, 11143, 11228, 10855, 9519, 7598, 4390, 2004, 769, 181, 
  9, 0, 0, 0, 0}, {0, 1, 155, 916, 2905, 6920, 14467, 18151, 19991, 
  20436, 19831, 18072, 14297, 7701, 3307, 1174, 253, 12, 0, 0, 0, 
  0}, {0, 2, 119, 1170, 3997, 10330, 22455, 28438, 30617, 31268, 
  30614, 28053, 23025, 11808, 4634, 1437, 214, 3, 0, 0, 0, 0}, {0, 0, 
  61, 1027, 4717, 13435, 30336, 38438, 41907, 42599, 41926, 39417, 
  31810, 16031, 5781, 1657, 145, 0, 0, 0, 0, 0}, {0, 0, 11, 674, 4489,
   15048, 36916, 46817, 50520, 51226, 51068, 48296, 39189, 18713, 
  6111, 1268, 47, 0, 0, 0, 0, 0}, {0, 0, 0, 280, 3692, 15404, 40139, 
  51198, 53790, 54718, 55258, 52828, 43694, 19608, 5208, 594, 0, 0, 0,
   0, 0, 0}, {0, 0, 0, 34, 2219, 13116, 38833, 48711, 51342, 51399, 
  52013, 51454, 43194, 18067, 3567, 184, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 
  1, 804, 9211, 32451, 41018, 41357, 42074, 43057, 43175, 36185, 
  13879, 1648, 15, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 155, 5049, 22938, 
  27689, 28143, 28324, 29185, 29578, 26021, 8580, 399, 0, 0, 0, 0, 0, 
  0, 0}, {0, 0, 0, 0, 7, 1882, 12109, 14079, 14010, 14468, 14745, 
  15047, 13718, 3798, 44, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 368, 
  3929, 4520, 4297, 4409, 4603, 4715, 4490, 1028, 0, 0, 0, 0, 0, 0, 0,
   0}, {0, 0, 0, 0, 0, 18, 475, 491, 468, 537, 499, 557, 553, 111, 0, 
  0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 4, 0, 0, 1, 0, 1, 4, 1, 0, 
  0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0}};
  
(*XYBinBoundariesOriginal = Flatten[Table[
    {
     {DetHisto[[1, xi]], DetHisto[[1, xi + 1]]},
     {-DetHisto[[2, yi + 1]] - R, -DetHisto[[2, yi]] - R}
     }, {xi, 1, Length[DetHisto[[1]]] - 1}, {yi, 1, 
     Length[DetHisto[[2]]] - 1}], 1];*)
     
XYBinBoundariesYHalf = Flatten[Table[
    {
     {DetHisto[[1, xi]], DetHisto[[1, xi + 1]]},
     {-DetHisto[[2, yi + 2]] - R, -DetHisto[[2, yi]] - R}
     }, 
     {xi, 1, Length[DetHisto[[1]]] - 1}, 
     {yi, 1, Length[DetHisto[[2]]] - 2, 2}], 1];

(*Print["XYBinLength = ",Length[XYBinBoundariesYHalf]];*)
 
IntMethod = {"GlobalAdaptive", Method -> "MultidimensionalRule"};

Print["Run started at ",DateString[],", Bins: ",Bins];
  
BinwYShiftPrec44OriginalBinsb0NewBGrad = 
 ParallelMap[
  AbsoluteTiming[
    BinIntShift[
     b, #, 
     {AlphaFromBlinesCorr, BfromCentralBline, rF, rRxB, rA, rDet, BlinefitG1, BlinefitG2}, {XAShift, YAShift, YRxBShift, XDetShift, YDetShift}, 
     {0.06, 0.059, 0.9, 0., -0.9, 0.1, 0.099, 0.9, 0., -0.9}, {ApertX, ApertY, ApertXOff+0.0002, ApertYOff},
     {IntMethod, 4, 2, 0, 2}, {IntMethod, 4, 2, 0, 2}, {IntMethod, 3, 6, 0, 1}
     ]] &, XYBinBoundariesYHalf[[Bins[[1]] ;; Bins[[2]]]], 
  Method -> "FinestGrained"];

Export["MergerTransport/NC_opt_dxAOff/b"<>ToString[IntegerPart[b*10000]]<>"/TransferResult_08-09-20_NC_opt_H_BSOn_dxAOff_b"<>ToString[IntegerPart[b*10000]]<>"_"<>ToString[Bins[[1]]]<>"-"<>ToString[Bins[[2]]]<>".txt",BinwYShiftPrec44OriginalBinsb0NewBGrad,"Table"];
 
Print["b = ",IntegerPart[b*10000],"e-4"];
Print["Time consumption: ",(AbsoluteTime[]-t0)/3600.];
Print["MaxMemUse :",MaxMemoryUsed[]/1024/1024.," MB"];
CloseKernels[];
Quit[];



