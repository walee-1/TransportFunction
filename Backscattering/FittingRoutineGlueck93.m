(* Wolfram Language package *)
Get["Backscattering.m"]
Get["Spectra/Glueck93ProtonSpectrum.m"]
Get["bins.m"]

bsCorrSin815FitFuncGlueck93[T_?NumericQ, a_?NumericQ, c1_?NumericQ, 
  c2_?NumericQ, c3_?NumericQ, c4_?NumericQ, c5_?NumericQ, rd_?NumericQ, postAcc_?NumericQ, thMax_?NumericQ] := 
 NIntegrate[
  wpNormedGlueck93[T, lambdaOfa[a], \[Kappa]0]*
   Sin[theta]*(1 - 
     fitFuncSin815ForFit[
       TofPClassical[pacc[pofTClassic[T, mp], postAcc]]/1000, 
       thDetAcc[pofTClassic[T, mp], theta, rd, postAcc]*180/Pi, c1, c2, 
       c3, c4, c5]/100), {theta, 0, thetaMax}]
       

FitFuncTableGlueck93[{a_, c1_, c2_, c3_, c4_, c5_,rd_,postAcc_,thMax_}] := 
 FitFuncTableGlueck93[{a, c1, c2, c3, c4, c5, rd, postAcc, thMax}] = ParallelTable[
   NIntegrate[
    Re[bsCorrSin815FitFuncGlueck93[T, a, c1, c2, c3, c4, c5, rd, postAcc, thMax]], {T,
     If[xbins[[i]] < 0, 0, xbins[[i]]]
     , If[xbins[[i + 1]] >= tpMax, tpMax, xbins[[i + 1]]]}, 
    PrecisionGoal -> 4, MinRecursion -> 2, MaxRecursion -> 10
    ], {i, Length[xbins] - 1}]
    
 fitFuncBinGlueck93[
  bin_?NumericQ, {a_?NumericQ, c1_?NumericQ, c2_?NumericQ, 
   c3_?NumericQ, c4_?NumericQ, c5_?NumericQ, rd_?NumericQ, postAcc_?NumericQ, thMax_?NumericQ}] := 
 FitFuncTableGlueck93[{a, c1, c2, c3, c4, c5, rd, postAcc, thMax}][[bin]]/
  Total[FitFuncTableGlueck93[{a, c1, c2, c3, c4, c5, rd, postAcc, thMax}]]
  
  
  fitFuncGlueck93BS[
  bin_?NumericQ, {a_?NumericQ, c1_?NumericQ, c2_?NumericQ, 
   c3_?NumericQ, c4_?NumericQ, c5_?NumericQ, rd_?NumericQ, postAcc_?NumericQ, thMax_?NumericQ}] := Piecewise[{
   {c1, bin == -1},
   {c2, bin == -2},
   {c3, bin == -3},
   {c4, bin == -4},
   {c5, bin == -5},
   {fitFuncBinGlueck93[bin, {a, c1, c2, c3, c4, c5, rd, postAcc, thMax}], bin > 0}}
  ]
  
FitFuncGlueck93BSOFF[T_?NumericQ, a_?NumericQ, thMax_?NumericQ] := 
 NIntegrate[
  wpNormedGlueck93[T, lambdaOfa[a], \[Kappa]0]*Sin[theta], {theta, 0, 
   thMax}]
  
FitFuncTableGlueck93BSOFF[{a_,thMax_}] := 
 FitFuncTableGlueck93BSOFF[{a,thMax}] = ParallelTable[
   NIntegrate[Re[FitFuncGlueck93BSOFF[T, a, thMax]], {T,
     If[xbins[[i]] < 0, 0, xbins[[i]]]
     , If[xbins[[i + 1]] >= tpMax, tpMax, xbins[[i + 1]]]}
    ], {i, Length[xbins] - 1}]  
  
  
fitFuncBinGlueck93BSOFF[bin_?NumericQ, {a_?NumericQ, thMax_?NumericQ}] := 
 FitFuncTableGlueck93BSOFF[{a, thMax}][[bin]]/Total[FitFuncTableGlueck93BSOFF[{a, thMax}]] 
  
  
  