(* Wolfram Language package *)

ApertBoole[xtraj_, ytraj_, {xA_, yA_, xOff_, yOff_}] := 
 Boole[xtraj < xOff + xA/2 && xtraj > xOff - xA/2 && 
   ytraj < yOff + yA/2 && ytraj > yOff - yA/2]
   
ApertBooleCompiled = Compile[
  {{xtraj, _Real}, {ytraj, _Real}, {xA, _Real}, {yA, _Real}, {xOff, _Real}, {yOff, _Real}},
  Boole[xtraj < xOff + xA/2 && xtraj > xOff - xA/2 && ytraj < yOff + yA/2 && ytraj > yOff - yA/2], 
  CompilationTarget -> "C", RuntimeAttributes -> {Listable}, RuntimeOptions -> "Speed"
  ];   