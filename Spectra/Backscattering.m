(* Wolfram Language package *)

Get["Common/CommonFunctions.m"];

bsFuncElectron[thDet_, Te_] := 0.10596090314315686` - 0.00010437408424908431` Te + 9.635959369457236`*^-7 thDet^3

bsFuncElectron2[p_, th0_, rD_] = bsFuncElectron[theta2[th0,rD], Tofpe[p]]