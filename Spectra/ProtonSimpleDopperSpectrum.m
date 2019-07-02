(* Wolfram Language package *)

Get["ProtonSpectrumNachtmann.m"];
Needs["Common/Constants.m"];

Tn0 = 4*10^-3;

TLAB[TCMS_, thCMS_, Tn_, mp_, mn_] := 
 TCMS + mp/mn*Tn + 2*Sqrt[mp*TCMS*Tn/mn]*Sin[thCMS]
 
tcms[tlab_, thCMS_, Tn_, mpp_, mnn_] = 
 FullSimplify[Solve[tlab == TLAB[TCMS, thCMS, Tn, mpp, mnn], TCMS]][[
  1, 1, 2]]
  
  
  
  dTCMSoverdTLAB[tlab_, thCMS_, Tn_, mpp_, mnn_] = 
 FullSimplify[D[tcms[tlab, thCMS, Tn, mpp, mnn], tlab], 
  Assumptions -> 
   0 <= thCMS <= Pi && Tn > 0 && mpp > 0 && mnn > 0 && tlab > 0]
   
   dTLABoverdpLAB[pLAB_] = D[TofP[pLAB], pLAB]
   
   ProtonMomSimpleDoppler[plab_, a_, thCMS_, Tn_] := 
 dwdt[tcms[TofP[plab], thCMS, Tn, mp, mn], a]*
  dTCMSoverdTLAB[TofP[plab], thCMS, Tn, mp, mn]*dTLABoverdpLAB[plab]
  
  (* ::Section:: *)
  (* limits of integration of momentum *)
  
  PLAB[TCMS_, thCMS_, Tn_] := 
 pofTClassic[TLAB[TCMS, thCMS, Tn, mp, mn], mp]

pmaxSimpleDoppler[thCMS_, Tn_] := 
 pofTClassic[TLAB[tpMax, thCMS, Tn, mp, mn], mp]

pminSimpleDoppler[thCMS_, Tn_] := 
 pofTClassic[TLAB[0., thCMS, Tn, mp, mn], mp]
 
 ProtonMomSimpleDopplerNorm[a_, thCMS_, Tn_] := 
 ProtonMomSimpleDopplerNorm[a, thCMS, Tn] = 
  NIntegrate[
   ProtonMomSimpleDoppler[plab, a, thCMS, Tn], {plab, 
    pminSimpleDoppler[thCMS, Tn], pmaxSimpleDoppler[thCMS, Tn]}]

ProtonMomSimpleDopplerNormed[plab_, a_, thCMS_, Tn_] := 
 ProtonMomSimpleDoppler[plab, a, thCMS, Tn]/
  ProtonMomSimpleDopplerNorm[a, thCMS, Tn]
  
  
  
  
  