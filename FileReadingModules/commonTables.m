(* Wolfram Language package *)
energyTable = {5, 10, 15, 20, 25, 30};
   
anglesTable = {0, 5, 10, 15, 25, 35, 45};

angleFunc[it_] := it*5 - 5;

energyIteratorFunc[n_] := 
  1/120*(-n^5 + 20 n^4 - 155 n^3 + 580 n^2 - 444 n + 120);