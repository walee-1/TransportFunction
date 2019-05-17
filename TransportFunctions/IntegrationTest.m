(* Wolfram Language package *)

Get["TransportFunctions/Integrands.m"];



IntegrationTest[opts : OptionsPattern[NIntegrate]] := Module[
  {t0 = AbsoluteTime[], t1, intresult},
  intresult = Reap@NIntegrate[
     IntegrandApert[
      0., {-0.05, p, th, Pi, 0.2, 1., 1.}, {0.01, 0.035, 0., 0., x0, 0., 1., 3}],
     {th,0,thetamax[2.]},
     {x0, -0.01/2 - rG[pmax, th, 0.2], 0.01/2 + rG[pmax, th, 0.2]},
     {
      p,
      pminCases[x0, -0.05, th, 0.2, Pi, th, 1., 0.2],
      pmaxCases[x0, -0.05, th, 0.2, Pi, th, 1., 0.2]
      },
     opts, EvaluationMonitor :> Sow[{x0, p}]
     ];
  t1 = AbsoluteTime[];
  {intresult[[1]],
   Histogram3D[
    Transpose[{intresult[[2, 1, All, 1]], 
      intresult[[2, 1, All, 2]]/1000}], 20, ImageSize -> 700, 
    AxesLabel -> {"x0 [m]", "p [keV/c", "SampleP"}, 
    ViewPoint -> {-1, -1, 1}],
   t1 - t0,
   Length[intresult[[2, 1]]]
   }
  ]


