(* Wolfram Language Raw Program *)
Get["Common/Constants.m"]
gausHalf[th_, sigma_] := 
 Sqrt[2]/(sigma*Sqrt[Pi]) E^(-th^2/(2*sigma^2))
 
pacc[p_, U_] := Sqrt[2*mp*(p^2/(2*mp) - U)]
thDetAcc[p_, th0_, rD_, U_] := 
 ArcSin[Sqrt[rD*p^2*Sin[th0]^2/(2*mp*(p^2/(2*mp) - U))]]
 	
FuncBackScat[en_?NumericQ, th_?NumericQ, c1_?NumericQ, c2_?NumericQ, c3_?NumericQ, 
  c4_?NumericQ, 
  c5_?NumericQ] := (c1*en + c2)*(1 - c3*gausHalf[th, c5] + (th/c4)^2)
  
