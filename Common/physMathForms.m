(* Wolfram Language package *)

Get["Common/Constants.m"];

critAngLindhard[Z1_,Z2_,e_,d_]:=Sqrt[2*Z1*Z2*1.44/(e*1000*d)]*180/Pi

critAngLindhard::usage = "critAngLindhard[Z1, Z2, e, d] gives the lindhard critical angle in degrees for channeling.
Arguments: 
	Z1, Z2 = Atomic numbers of the ion and material. 
	e = Kinetic energy of the impinging ion in keV. 
	d = interplanar distance between crystal planes.
For interplanar distance use: interplanarDist"

interplanarDist[h_,k_,l_,a_:0.543]:=a/Sqrt[h^2+k^2+l^2]

interplanarDist::usage = "interplanarDist[h, k, l ,a] gives interplanar distance between two crystalline planes for a specific crystal orientation 
Arguments: 
	h, k, l = miller indices of the crystal
	a = distance between two crystal planes, a is set to the default value of 0.543 nm (Si)"

betaCalc[mass_,ke_]:=Sqrt[1-mass^2/(mass+ke)^2]

betaCalc::usage = "betaCalc[mass_,ke_] calculates the beta factor for a given mass and ke."

errorBarMod[prob_, totalN_] := 
 100/totalN*
  Sqrt[totalN*
    prob*(1 - 
      prob)] (*because we convert the mean to a percentage by \
dividing by totalN and multiplying by 100, 
  so similar thing should be done for the errorbars... basic \
propogation of errors is a thing dumbass*)
errorBarMod::usage = "errorBarMod[prob_, totalN_] gives the binomial error bar for percentage of N particles in a bin"
    statErrorCalc[percent_, TotalNumber_] := 
 percent/Sqrt[(percent*TotalNumber)]
 
 statErrorCalc::usage = "statErrorCalc[percent_, TotalNumber_] gives the statistical error bar for percentage of N particles in a bin"
 
 relErrCalc[valEXP_, valOrig_] := (valEXP - valOrig)/valOrig
 
 relErrCalc::usage = "relErrCalc[valEXP_, valOrig_] gives the relative error between two numbers"
 
     PercCalc[val_, per_, round_: 0.001] := 
 Block[{err, return}, err = val*per; 
  return = Sort[{val + err, val - err}]; Return[Round[return, round]]]
  
  PercCalc::usage = "PercCalc[val_, per_, round_: 0.001] gives the error bar values of a number when the error is given in a percent around the number."
  
  
  