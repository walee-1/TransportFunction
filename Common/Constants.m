(* Wolfram Language package *)

\[Lambda]0 = -1.2732;
\[Kappa]0=1.85;
a0 = -0.106;
me = 510998.95; 
c = 299792458;
hq = 1.054572*10^-34;
mp = 938.27208 * 10^6;
mn = 939.56541 * 10^6;
del = 1.293332*10^6;
xr =  me/del;
epMax = mp + (del^2 - me^2)/(2*mn);
E0 = del - (del^2 - me^2)/(2*mn);
tpMax = epMax - mp;
\[Alpha] = 1/137;


SetOptions[ListPlot,GridLines->Automatic,Frame->True,ImageSize->700,FrameStyle->20];
SetOptions[Plot,GridLines->Automatic,Frame->True,ImageSize->700,FrameStyle->20];