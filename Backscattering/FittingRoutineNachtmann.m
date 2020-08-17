(* Wolfram Language package *)
Get["Backscattering.m"]
Get["Spectra/ProtonSpectrumNachtmann.m"]
Get["bins.m"]


bsCorrSin815FuncNacht[T_?NumericQ, a_?NumericQ, c1_?NumericQ, 
  c2_?NumericQ, c3_?NumericQ, c4_?NumericQ, c5_?NumericQ, 
  rd_?NumericQ, postAcc_?NumericQ, thMax_?NumericQ] := 
 NIntegrate[
  dwdt[T, a]*
   Sin[theta]*(1 - 
     fitFuncSin815ForFit[TofPClassical[pacc[pofTClassic[T, mp], postAcc]]/1000, 
       thDetAcc[pofTClassic[T, mp], theta, rd, postAcc]*180/Pi, c1, c2, c3, c4, 
       c5]/100), {theta, 0, thMax}]

FitFuncTableNacht[{a_, c1_, c2_, c3_, c4_, c5_, rd_, postAcc_, 
   thMax_}] := 
 FitFuncTableNacht[{a, c1, c2, c3, c4, c5, rd, postAcc, thMax}] = 
  ParallelTable[
   NIntegrate[
    Re[bsCorrSin815FuncNacht[T, a, c1, c2, c3, c4, c5, rd, postAcc, 
      thMax]], {T,
     If[xbins[[i]] < 0, 0, xbins[[i]]]
     , If[xbins[[i + 1]] >= tpMax, tpMax, xbins[[i + 1]]]}
    ], {i, Length[xbins] - 1}]
    
fitFuncBinNacht[
  bin_?NumericQ, {a_?NumericQ, c1_?NumericQ, c2_?NumericQ, 
   c3_?NumericQ, c4_?NumericQ, c5_?NumericQ, rd_?NumericQ, 
   postAcc_?NumericQ, thMax_?NumericQ}] := 
 FitFuncTableNacht[{a, c1, c2, c3, c4, c5, rd, postAcc, 
     thMax}][[bin]]/
  Total[FitFuncTableNacht[{a, c1, c2, c3, c4, c5, rd, postAcc, 
     thMax}]]
     
     
fitFuncNacht[
  bin_?NumericQ, {a_?NumericQ, c1_?NumericQ, c2_?NumericQ, 
   c3_?NumericQ, c4_?NumericQ, c5_?NumericQ, rd_?NumericQ, 
   postAcc_?NumericQ, thMax_?NumericQ}] := Piecewise[{
   {c1, bin == -1},
   {c2, bin == -2},
   {c3, bin == -3},
   {c4, bin == -4},
   {c5, bin == -5},
   {fitFuncBinNacht[bin, {a, c1, c2, c3, c4, c5, rd, postAcc, thMax}],
     bin > 0}}
  ]
  
FitFuncNachtBSOFF[T_?NumericQ, a_?NumericQ, thMax_?NumericQ] := 
 NIntegrate[dwdt[T, a]*Sin[theta], {theta, 0, thMax}]
 
 
FitFuncTableNachtBSOFF[{a_, thMax_}] := 
 FitFuncTableNachtBSOFF[{a, thMax}] = ParallelTable[
   NIntegrate[Re[FitFuncNachtBSOFF[T, a, thMax]], {T,
     If[xbins[[i]] < 0, 0, xbins[[i]]]
     , If[xbins[[i + 1]] >= tpMax, tpMax, xbins[[i + 1]]]}
    ], {i, Length[xbins] - 1}]  
    
    
fitFuncBinNachtBSOFF[bin_?NumericQ, {a_?NumericQ,thMax_?NumericQ}] := 
 FitFuncTableNachtBSOFF[{a, thMax}][[bin]]/
  Total[FitFuncTableNachtBSOFF[{a, thMax}]]
