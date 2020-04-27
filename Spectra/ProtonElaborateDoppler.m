(* Wolfram Language package *)

Get["Common/Constants.m"];

thLab[TCMS_, thCMS_, Tn_, mp_, mn_] := 
 ArcTan[Cos[thCMS]/(Sin[thCMS] + Sqrt[mp/mn*Tn/TCMS])]
 
 thLAB[TCMS_, thCMS_, Tn_, mp_, mn_] := 
 ArcTan[(Sin[thCMS] + Sqrt[mp/mn*Tn/TCMS])/Cos[thCMS]]
 
 
 TCMS[tlab_, thlab_] := 
 Solve[
 	tlab == TLAB[TCMS, thCMSs, 4*10^-3, mp, mn] && thlab == thLAB[TCMS, thCMSs, 4*10^-3, mp, mn] && thCMSs < Pi,
 	{TCMS, thCMSs}, Reals][[1, 1, 2, 1]]
    
    
    thCMS[tlab_, thlab_] := 
 Solve[tlab == TLAB[TCMS, thCMSs, 4*10^-3, mp, mn] && 
    thlab == thLAB[TCMS, thCMSs, 4*10^-3, mp, mn] && 
    thCMSs < Pi, {TCMS, thCMSs}, Reals][[1, 2, 2, 1, 1]]