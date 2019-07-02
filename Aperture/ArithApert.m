(* Wolfram Language package *)


Get["Common/CommonFunctions.m"];


IntersectionAngles[x_] := Piecewise[
   {
    (*{{0,2*Pi},x===ComplexInfinity},*)
    {{0, 0}, x <= -1},
    {{0, 2 \[Pi]}, x >= 1},
    {{ArcCos[x], 2 \[Pi] - ArcCos[x]}, True}
    }
   ];

IntervalMod[{min_, max_}] := Interval @@ Piecewise[
    {
     {{{0, max - 2 \[Pi]}, {min, 2 \[Pi]}}, min <= 2 \[Pi] <= max},
     {{Mod[{min, max}, 2 \[Pi]]}, True}
     }];

ArithApert[xA_, yA_, xOff_, yOff_, xGC_, yGC_, r_] :=
 
 Piecewise[(*we need a piecewise for r=0*)
  {
   {1, r <10^-12 && -xA/2 + xOff <= xGC <= xA/2 + xOff && -yA/2 + yOff <=
       yGC <= yA/2 + yOff},
   {0, r <10^-12(*&&(-xA/2+xOff>xGC||xA/2+xOff<xGC)&&(-yA/2+yOff>yGC||yA/
    2+yOff<yGC)*)}(*this is sufficient as it will first check first \
case, then second, so it has to be outside in this case*)
   },
  Module[{arr = {xOff - xGC + xA/2, 
       yOff - yGC + yA/2, -(xOff - xGC - xA/2), -(yOff - yGC - yA/2)}/
      r, int},
   (*{arr,
   Table[IntersectionAngles[arr[[i]]]+(i-1)\[Pi]/2,{i,4}],
   Table[IntervalMod[IntersectionAngles[arr[[i]]]+(i-1)\[Pi]/2],{i,
   4}],*)
   
   int = IntervalIntersection @@ 
     Table[IntervalMod[
       IntersectionAngles[arr[[i]]] + (i - 1) \[Pi]/2], {i, 4}];
    RegionMeasure[int, 1]/(2 Pi)(*,
   r RegionMeasure[int]
   }*)
   ]
  ]