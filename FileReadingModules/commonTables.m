(* Wolfram Language package *)
energyTable = {5, 10, 15, 20, 25, 30};
   
anglesTable = {0, 5, 10, 15, 25, 35, 45};

angleFunc[it_] := it*5 - 5;

energyIteratorFunc[n_] := 
  1/120*(-n^5 + 20 n^4 - 155 n^3 + 580 n^2 - 444 n + 120);

energyTableFine15 = {15}~Join~Table[i, {i, 15.1, 15.8, 0.1}]



energyTableFine30 = {30}~Join~Table[i, {i, 30.1, 30.8, 0.1}]

{30, 30.1, 30.2, 30.3, 30.4, 30.5, 30.6, 30.7, 30.8}

angleTableForFit2 = Table[i, {i, 0, 10}];

energyTableFiner15 = {15}~Join~
  Table[i, {i, 15.05, 15.8, 0.05}]

energyTableFiner30 = {30}~Join~
  Table[i, {i, 30.05, 30.8, 0.05}]

angTab10 = Table[i, {i, 0, 10}];

angTab13 = Table[i, {i, 0, 13}];

angTab10Gap2 = Table[i, {i, 0, 10, 2}];

angTab13Gap2 = {0, 2, 4, 6, 8, 10, 12, 13};

angTab12 = Table[i, {i, 0, 12}];

angTab12Gap2 = {0, 2, 4, 6, 8, 10, 12};

angTab10Gap2w5 = {0, 2, 4, 5, 6, 8, 10};

angTab18Gap2 = Table[i, {i, 0, 18, 2}];

tiltAngles = Table[i, {i, 0, 20, 1}]

anglesTable = {0, 5, 10, 15, 25, 35, 45};

phiTable = Table[i, {i, 0, 180, 30}];

anglesTable2 = {0, 2, 4, 5, 6, 8, 10, 15, 25, 35, 45};

energyTable = {5, 10, 15, 20, 25, 30};


