(* Wolfram Language package *)

Get["CommonFunctions.m"];

CircleInRectangle[xA_, yA_, xOff_, yOff_, x0_, y0_, r_] =
  ImplicitRegion[
   -xA/2 + xOff <= x <= xA/2 + xOff && -yA/2 + yOff <= y <= 
     yA/2 + yOff &&
    {x, y} \[Element] 
     Circle[{x0, y0}, r](*&&
   {x,y,xA,yA,xOff,yOff,x0,y0,
   r}\[Element]Reals*),
   {x, y}];



ApertureFunc[xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, 
  yOff_?NumericQ, x0_?NumericQ, y0_?NumericQ, th2_?NumericQ, 
  p_?NumericQ, BA_?NumericQ,Prec_?IntegerQ] :=
 Piecewise[
  {
   {
    Piecewise[
     {
      {1, -xA/2 + xOff <= x0 <= xA/2 + xOff && -yA/2 + yOff <= y0 <= 
         yA/2 + yOff}
      },
     0.],
    rG[p, th2, BA] < 10^-6
    }(*if rG small, 
   just one ore zero dependent on inside or outside*)
   },
  (*val zero*)
  Piecewise[
   {
    {(*completely outide of region*)
     0., 
     x0 + rG[p, th2, BA] < -xA/2 || x0 - rG[p, th2, BA] > xA/2 || 
      y0 + rG[p, th2, BA] < -yA/2 || y0 - rG[p, th2, BA] > yA/2
     },
    {(*completely inside of region*)
     1., 
     x0 - rG[p, th2, BA] >= -xA/2 && x0 + rG[p, th2, BA] <= xA/2 && 
      y0 - rG[p, th2, BA] >= -yA/2 && y0 + rG[p, th2, BA] <= yA/2
     }
    },
   NIntegrate[
     1, {x, y} \[Element] 
      CircleInRectangle[xA, yA, xOff, yOff, x0, y0, rG[p, th2, BA]], 
     PrecisionGoal -> Prec]/(2*Pi*rG[p, th2, BA])
   (*RegionMeasure[CircleInRectangle[xA,yA,xOff,yOff,x0,y0,rG[p,th2,
   BA]]]/(2*Pi*rG[p,th2,BA])*)
   ]
  ]
  
  
  
  MonteCarloAperture[xA_?NumericQ, yA_?NumericQ, xOff_?NumericQ, 
  yOff_?NumericQ, x0_?NumericQ, y0_?NumericQ, th2_?NumericQ, 
  p_?NumericQ, BA_?NumericQ, steps_?NumericQ] :=
 Module[
  {rGlocal = rG[p, th2, BA], step = 2*Pi/steps, XYTable},
  XYTable = 
   Table[{rGlocal*Cos[phi] + x0, rGlocal*Sin[phi] + y0}, {phi, 0, 
     2*Pi - step, step}];
  Piecewise[
   {
    {
     Piecewise[
      {
       {1, -xA/2 + xOff <= x0 <= xA/2 + xOff && -yA/2 + yOff <= y0 <= 
          yA/2 + yOff}
       },
      0.],
     rGlocal < 10^-6
     }(*if rG small, 
    just one ore zero dependent on inside or outside*)
    },
   (*val zero*)
   Piecewise[
    {
     {(*completely outide of region*)
      0., 
      x0 + rGlocal < -xA/2 || x0 - rGlocal > xA/2 || 
       y0 + rGlocal < -yA/2 || y0 - rGlocal > yA/2
      },
     {(*completely inside of region*)
      1., 
      x0 - rGlocal >= -xA/2 && x0 + rGlocal <= xA/2 && 
       y0 - rGlocal >= -yA/2 && y0 + rGlocal <= yA/2
      }
     },
    N[Count[
       XYTable, {_?(-xA/2 + xOff <= # <= 
            xA/2 + xOff &), _?(-yA/2 + yOff <= # <= yA/2 + yOff &)}]/
      steps]
    (*RegionMeasure[CircleInRectangle[xA,yA,xOff,yOff,x0,y0,rG[p,th2,
    BA]]]/(2*Pi*rG[p,th2,BA])*)
    ]
   ]
  ]
  
  
  
  