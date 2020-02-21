(* Wolfram Language package *)

\[Lambda]0 = -1.2724;
\[Kappa]0=1.85;
me = 510998.95; 
c = 299792458;
hq = 1.055*10^-34;
mp = 938.27208 * 10^6;
mn = 939.56541 * 10^6;
del = mn - mp;
xr =  me/del;
epMax = mp + (del^2 - me^2)/(2*mn);
E0 = del - (del^2 - me^2)/(2*mn);
tpMax = epMax - mp;
\[Alpha] = 1/137;

<<ErrorBarPlots`;
SetOptions[ListPlot,GridLines->Automatic,Frame->True,ImageSize->700,FrameStyle->20];
SetOptions[Plot,GridLines->Automatic,Frame->True,ImageSize->700,FrameStyle->20];