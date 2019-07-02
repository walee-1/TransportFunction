(* Wolfram Language package *)

Get["Common/Constants.m"];

Eofp[p_] := Sqrt[me^2 + p^2]

TofPClassical[p_] := p^2/(2*mp);

pofT[T_, m_] := Sqrt[T^2 + 2*T*m]
pofTClassic[T_, m_] := Sqrt[2*T*m]

thetamax[rB1_] := ArcSin[Sqrt[1/rB1]]

rG[p_, th2_, B_] := p*Sin[th2]/c/B

theta2[th0_, rB2_] := ArcSin[Sin[th0]*Sqrt[rB2]]