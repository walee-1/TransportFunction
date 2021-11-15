(* Wolfram Language package *)
SetOptions[ListPlot, ImageSize -> 800, Frame -> True, 
  GridLines -> Automatic, 
  LabelStyle -> Directive[FontSize -> 20, Bold, Black], 
  FrameStyle -> Thick, FrameTicksStyle -> Thick];

SetOptions[ListLinePlot, ImageSize -> 800, Frame -> True, 
  GridLines -> Automatic, 
  LabelStyle -> Directive[FontSize -> 20, Bold, Black], 
  FrameStyle -> Thick, FrameTicksStyle -> Thick];

SetOptions[Plot, ImageSize -> 800, Frame -> True, 
  GridLines -> Automatic, 
  LabelStyle -> Directive[FontSize -> 20, Bold, Black], 
  FrameStyle -> Thick, FrameTicksStyle -> Thick];

SetOptions[ListLogPlot, ImageSize -> 800, Frame -> True, 
  GridLines -> Automatic, 
  LabelStyle -> Directive[FontSize -> 20, Bold, Black], 
  FrameStyle -> Thick, FrameTicksStyle -> Thick];
  
  
scTicks[start_, end_, step_, power_] := 
 Table[{i*10^power, 
   Style[ToString[i] <> "\[Times]\!\(\*SuperscriptBox[\(10\), \(" <> 
     ToString[power] <> "\)]\)"]}, {i, start, end, step}]
     
 fixLegend[size_][legend_] := 
 ReplaceAll[
  legend, {g_Graphics :> 
    Show[g, ImageSize -> Automatic -> {Automatic, size}, 
     BaselinePosition -> Axis], Pane -> Function@#}]
     
MakeBoxes[stripGridAlignemnt[expr_], form_] ^:= 
 ReplaceAll[MakeBoxes[expr, form], 
  Rule[GridBoxAlignment, _] -> Sequence[]]
  
 plotter[list_, NumberOfTicks_, interOrder_: Automatic, 
  opts : OptionsPattern[]] := 
 Module[{min, max, plot}, min = Min[list[[All, 3]]]; 
  max = Max[list[[All, 3]]]; 
  plot = ListDensityPlot[list, PlotRange -> Full, 
    InterpolationOrder -> interOrder, opts, 
    ColorFunctionScaling -> False, 
    ColorFunction -> (MPLColorMap["Viridis"][
        LogarithmicScaling[#, min, max]] &), 
    PlotLegends -> 
     BarLegend[{MPLColorMap["Viridis"], {0, 1}}, 
      LegendMarkerSize -> 370, 
      Ticks -> ({LogarithmicScaling[#, min, max], 
           ScientificForm[#, 2]} & /@ (min (max/min)^
            Range[0, 1, 1/NumberOfTicks]))]]]

LogScaleLegend[min_, max_, colorfunction_, height_: 400] := 
 Module[{bareTicksList, numberedTicks, m, M, ml, Ml, 
   minInArbitraryScale, maxInArbitraryScale, linearScaling}, 
  bareTicksList = 
   First[Ticks /. AbsoluteOptions[LogLogPlot[x, {x, min, max}]]];
  numberedTicks = (Select[bareTicksList /. {Superscript -> Power}, 
      NumberQ[#[[2]]] &])[[All, {1, 2}]];
  m = Min[numberedTicks[[All, 2]]];
  M = Max[numberedTicks[[All, 2]]];
  ml = Min[numberedTicks[[All, 1]]];
  Ml = Max[numberedTicks[[All, 1]]];
  {minInArbitraryScale, maxInArbitraryScale} = 
   ml + (Ml - ml) Log[{min, max}/m]/Log[M/m];
  linearScaling[
    x_] := (x - minInArbitraryScale)/(maxInArbitraryScale - 
      minInArbitraryScale);
  DensityPlot[y, {x, 0, 0.04}, {y, 0, 1}, AspectRatio -> Automatic, 
   PlotRangePadding -> 0, ImageSize -> {Automatic, height}, 
   ColorFunction -> colorfunction, 
   FrameTicks -> {{None, 
      Select[Table[{linearScaling[r[[1]]], 
         r[[2]] /. {Superscript[10., n_] -> Superscript[10, n]}, {0, 
          If[r[[2]] === "", 0.15, 0.3]}, {If[r[[2]] === "", 
           Thickness[0.03], Thickness[0.06]]}}, {r, 
         bareTicksList}], (#[[1]] (1 - #[[1]]) >= 0 &)]}, {None, 
      None}}]]

plotterCustomContour[list_, min_, max_, NumberOfTicks_, 
  interOrder_: Automatic, opts : OptionsPattern[]] := 
 Module[{plot}, 
  plot = ListContourPlot[list, PlotRange -> Full, 
    InterpolationOrder -> interOrder, opts, 
    ColorFunctionScaling -> False, 
    ColorFunction -> (MPLColorMap["Viridis"][
        LogarithmicScaling[#, min, max]] &), 
    PlotLegends -> 
     LogScaleLegend[min, max, MPLColorMap["Viridis"], 350]]]

plotter2[list_, NumberOfTicks_, interOrder_: Automatic, 
  opts : OptionsPattern[]] := 
 Module[{min, max, plot}, min = Min[list[[All, 3]]]; 
  max = Max[list[[All, 3]]]; 
  plot = ListDensityPlot[list, PlotRange -> Full, 
    InterpolationOrder -> interOrder, opts, 
    ColorFunctionScaling -> False, 
    ColorFunction -> (MPLColorMap["Viridis"][
        LogarithmicScaling[#, min, max]] &), 
    PlotLegends -> 
     LogScaleLegend[min, max, ColorData["DeepSeaColors"], 350]]]

plotterContourTest[list_, min_, max_, NumberOfTicks_, 
  interOrder_: Automatic, opts : OptionsPattern[]] := 
 Module[{plot}, 
  plot = ListContourPlot[list, PlotRange -> {min, max}, 
    InterpolationOrder -> interOrder, opts, 
    ColorFunctionScaling -> False, 
    ColorFunction -> (MPLColorMap["Viridis"][
        LogarithmicScaling[#, min, max]] &), 
    PlotLegends -> 
     BarLegend[{MPLColorMap["Viridis"], {0, 1}}, 
      LegendMarkerSize -> 370, 
      Ticks -> ({LogarithmicScaling[#, min, max], 
           ScientificForm[#, 2]} & /@ (min (max/min)^
            Range[0, 1, 1/NumberOfTicks])), TicksStyle -> Thick]]]

myPlotter[list_, interOrder_: Automatic, opts : OptionsPattern[]] := 
 Module[{min, max, plot}, min = Min[list[[All, 3]]]; 
  max = Max[list[[All, 3]]]; 
  plot = ListDensityPlot[list, 
    ScalingFunctions -> {"Linear", "Linear", "Log10"}, opts, 
    ColorFunction -> JetCM, 
    PlotLegends -> 
     BarLegend[{JetCM, {min, max}}, LabelStyle -> Directive[Black], 
      LegendMarkerSize -> 350], InterpolationOrder -> interOrder]]

myPlotter[list_, interOrder_: Automatic, min_, max_, 
  opts : OptionsPattern[]] := 
 Module[{plot}, 
  plot = ListDensityPlot[list, 
    ScalingFunctions -> {"Linear", "Linear", "Log10"}, opts, 
    ColorFunction -> JetCM, 
    PlotLegends -> 
     BarLegend[{JetCM, {min, max}}, LabelStyle -> Directive[Black], 
      LegendMarkerSize -> 350], InterpolationOrder -> interOrder]]

stereoPlotterMod[list_, colorScheme_, radian_: 32, 
  padding_: {{1.7, 5}, {1, 5}}] := 
 Module[{data, func, rmin, rmax, phimin, phimax, dphi, plot},
  data = Transpose[{list[[All, 1]], list[[All, 2]] Degree, 
     list[[All, 3]]}];
  func = Interpolation[data];
  {{rmin, rmax}, {phimin, phimax}} = func["Domain"];
  dphi = 15 Degree;
  plot = Legended[
    Show[ParametricPlot[
      r {Cos[phi], Sin[phi]}, {r, rmin, rmax}, {phi, phimin, phimax}, 
      FrameLabel -> {"Polar Angle \[Theta](\[Degree])", 
        "Polar Angle \[Theta](\[Degree])"}, 
      ColorFunction -> 
       Function[{x, y, r, phi}, 
        colorScheme@Rescale[func[r, phi], MinMax@data[[All, 3]]]], 
      ColorFunctionScaling -> False, PlotPoints -> 50, 
      LabelStyle -> Directive[FontSize -> 20, Bold, Black], 
      Frame -> {{True, True}, {True, True}}, FrameTicksStyle -> Thick,
       FrameStyle -> Thick, PlotRangePadding -> padding, 
      PlotRange -> {{0, rmax}, {0, rmax}}
      ], PolarPlot[radian, {phi, phimin, phimax}, 
      PolarAxes -> Automatic, AxesStyle -> {Black, Thick}, 
      PolarTicks -> {"Degrees", None}, 
      TicksStyle -> Directive[Black, Bold, FontSize -> 20], 
      PolarGridLines -> {Range[phimin, phimax, dphi], Automatic}, 
      LabelStyle -> Directive[FontSize -> 20, Bold, Black], 
      GridLinesStyle -> Directive[{Black, Thick, Dashed}], 
      PlotRange -> All, PlotStyle -> Transparent
      ], ImageSize -> 600, 
     LabelStyle -> Directive[FontSize -> 20, Bold, Black]], 
    BarLegend[{colorScheme@Rescale[#, MinMax@data[[All, 3]]] &, 
      MinMax@data[[All, 3]]}, 
     LabelStyle -> Directive[FontSize -> 20, Bold, Black], 
     LegendMarkerSize -> 500, 
     LegendLabel -> 
      Placed["Backscattering Rate (%)", After, Rotate[#, Pi/2] &], 
     TicksStyle -> Directive[Black, Thick]]
    ];
  Return[plot]
  ]
     
     
     
stereoPlotterModEqualLegend[list_, colorScheme_,minmax_:-1, radian_: 32, 
  padding_: {{1.7, 5}, {1, 5}}] := 
 Module[{data, func, rmin, rmax, phimin, phimax, dphi, plot},
  data = Transpose[{list[[All, 1]], list[[All, 2]] Degree, 
     list[[All, 3]]}];
  if[minmax==-1,minmax=MinMax@data[[All,3]]];
  func = Interpolation[data];
  {{rmin, rmax}, {phimin, phimax}} = func["Domain"];
  dphi = 15 Degree;
  plot = Legended[
    Show[ParametricPlot[
      r {Cos[phi], Sin[phi]}, {r, rmin, rmax}, {phi, phimin, phimax}, 
      FrameLabel -> {"Polar Angle \[Theta](\[Degree])", 
        "Polar Angle \[Theta](\[Degree])"}, 
      ColorFunction -> 
       Function[{x, y, r, phi}, 
        colorScheme@Rescale[func[r, phi], minmax]], 
      ColorFunctionScaling -> False, PlotPoints -> 50, 
      LabelStyle -> Directive[FontSize -> 20, Bold, Black], 
      Frame -> {{True, True}, {True, True}}, FrameTicksStyle -> Thick,
       FrameStyle -> Thick, PlotRangePadding -> padding, 
      PlotRange -> {{0, rmax}, {0, rmax}}
      ], PolarPlot[radian, {phi, phimin, phimax}, 
      PolarAxes -> Automatic, AxesStyle -> {Black, Thick}, 
      PolarTicks -> {"Degrees", None}, 
      TicksStyle -> Directive[Black, Bold, FontSize -> 20], 
      PolarGridLines -> {Range[phimin, phimax, dphi], Automatic}, 
      LabelStyle -> Directive[FontSize -> 20, Bold, Black], 
      GridLinesStyle -> Directive[{Black, Thick, Dashed}], 
      PlotRange -> All, PlotStyle -> Transparent
      ], ImageSize -> 600, 
     LabelStyle -> Directive[FontSize -> 20, Bold, Black]], 
    BarLegend[{colorScheme@Rescale[#, minmax] &, 
      minmax}, 
     LabelStyle -> Directive[FontSize -> 20, Bold, Black], 
     LegendMarkerSize -> 500, 
     LegendLabel -> 
      Placed["Backscattering Rate (%)", After, Rotate[#, Pi/2] &], 
     TicksStyle -> Directive[Black, Thick]]
    ];
  Return[plot]
  ]
     
     