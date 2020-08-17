binNo=64;

xbins = Table[
   bin, {bin, 0, Ceiling[tpMax, Ceiling[tpMax/binNo]], 
    Ceiling[tpMax/binNo]}];(* Wolfram Language package *)
    
 (*for Nomos*)
rd = 1;

postAcc = -15000;

thMax = Pi/4