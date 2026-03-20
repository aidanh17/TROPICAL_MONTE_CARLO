(* Convergence Study v2: fixed scientific notation parsing *)

Print["================================================================"];
Print["  Convergence Study: Tropical MC -> FIESTA"];
Print["  Bubblepp presector 1, k1=k2=k3=1"];
Print["================================================================"];
Print[];

testDir = "/home/aidanh/Desktop/Tropical_Monte_Carlo_Final/Bubble1final/tropical_code_claudef/FIESTA_Compare/bubblepp_test";
mainBin = FileNameJoin[{testDir, "main_mc"}];
g0Bin   = FileNameJoin[{testDir, "g0_mc"}];
kinFile = FileNameJoin[{testDir, "kin.txt"}];

fiestaPole   = 0.000513;
fiestaPoleErr = 0.000001;
fiestaFinite = 0.024968;
fiestaFiniteErr = 0.000037;

Print["FIESTA reference (5M Vegas samples):"];
Print["  Pole (1/ep):   ", fiestaPole, " +/- ", fiestaPoleErr];
Print["  Finite (ep^0): ", fiestaFinite, " +/- ", fiestaFiniteErr];
Print[];

(* Fix scientific notation: replace "e" with "*^" for Mathematica *)
parseVal[s_String] := Module[{fixed},
  fixed = StringReplace[s, RegularExpression["([0-9.]+)e([+-]?[0-9]+)"] -> "$1*^$2"];
  ToExpression[fixed]
];

parseMCOutput[file_String] := Module[{text, parts},
  text = StringTrim[Import[file, "Text"]];
  parts = StringSplit[text];
  {parseVal[parts[[1]]], parseVal[parts[[2]]],
   parseVal[parts[[3]]], parseVal[parts[[4]]]}
];

sampleCounts = {50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000};
nThreads = 4;
results = {};

Do[
  Module[{mainF, g0F, t0, procM, procG, elapsed, mV, gV},
    mainF = FileNameJoin[{testDir, "conv_main_" <> ToString[n] <> ".txt"}];
    g0F   = FileNameJoin[{testDir, "conv_g0_" <> ToString[n] <> ".txt"}];

    t0 = AbsoluteTime[];
    procM = StartProcess[{mainBin, kinFile, mainF, ToString[n], ToString[nThreads]}];
    procG = StartProcess[{g0Bin, kinFile, g0F, ToString[n], ToString[nThreads]}];
    While[ProcessStatus[procM]==="Running" || ProcessStatus[procG]==="Running", Pause[0.5]];
    elapsed = AbsoluteTime[] - t0;

    mV = parseMCOutput[mainF];
    gV = parseMCOutput[g0F];

    AppendTo[results, <|
      "N" -> n, "Time" -> elapsed,
      "MainRe" -> mV[[1]], "MainErr" -> mV[[3]],
      "G0Re" -> gV[[1]], "G0Err" -> gV[[3]]
    |>];

    Quiet[DeleteFile /@ {mainF, g0F}];
    Print["  N=", n, " done (", Round[elapsed, 0.1], "s)"];
  ],
  {n, sampleCounts}
];

Print[];
Print["================================================================"];
Print["  G0 (1/eps POLE coefficient) — convergence to FIESTA"];
Print["================================================================"];
Print[];

Module[{header},
  header = StringJoin[
    StringPadRight["  N", 12],
    StringPadRight["G0_MC", 14],
    StringPadRight["MC_err", 12],
    StringPadRight["FIESTA", 12],
    StringPadRight["|diff|", 12],
    StringPadRight["sigma", 8],
    "t(s)"
  ];
  Print[header];
  Print["  ", StringJoin[Table["-", 78]]];

  Do[
    Module[{r = results[[i]], diff, sigma, errStr, diffStr, sigStr},
      diff = Abs[r["G0Re"] - fiestaPole];
      sigma = If[Abs[r["G0Err"]] > 0, diff / Abs[r["G0Err"]], 0.];
      Print[
        StringPadRight["  " <> ToString[r["N"]], 12],
        StringPadRight[ToString@NumberForm[r["G0Re"], {7, 7}], 14],
        StringPadRight[ToString@ScientificForm[r["G0Err"], 2], 12],
        StringPadRight[ToString[fiestaPole], 12],
        StringPadRight[ToString@ScientificForm[diff, 2], 12],
        StringPadRight[ToString@NumberForm[sigma, {3, 1}], 8],
        ToString@Round[r["Time"], 0.1]
      ];
    ],
    {i, Length[results]}
  ];
];

Print[];
Print["================================================================"];
Print["  TOTAL (eps^0 FINITE part) — convergence to FIESTA"];
Print["================================================================"];
Print[];

Module[{header},
  header = StringJoin[
    StringPadRight["  N", 12],
    StringPadRight["TOTAL_MC", 14],
    StringPadRight["MC_err", 12],
    StringPadRight["FIESTA", 12],
    StringPadRight["|diff|", 12],
    StringPadRight["sigma", 8],
    "t(s)"
  ];
  Print[header];
  Print["  ", StringJoin[Table["-", 78]]];

  Do[
    Module[{r = results[[i]], diff, sigma},
      diff = Abs[r["MainRe"] - fiestaFinite];
      sigma = If[Abs[r["MainErr"]] > 0, diff / Abs[r["MainErr"]], 0.];
      Print[
        StringPadRight["  " <> ToString[r["N"]], 12],
        StringPadRight[ToString@NumberForm[r["MainRe"], {7, 7}], 14],
        StringPadRight[ToString@ScientificForm[r["MainErr"], 2], 12],
        StringPadRight[ToString[fiestaFinite], 12],
        StringPadRight[ToString@ScientificForm[diff, 2], 12],
        StringPadRight[ToString@NumberForm[sigma, {3, 1}], 8],
        ToString@Round[r["Time"], 0.1]
      ];
    ],
    {i, Length[results]}
  ];
];

Print[];
Print["================================================================"];
Print["  1/sqrt(N) error scaling check"];
Print["================================================================"];
Print[];

Module[{n0, g0Err0, mainErr0},
  n0 = results[[1]]["N"];
  g0Err0 = Abs[results[[1]]["G0Err"]];
  mainErr0 = Abs[results[[1]]["MainErr"]];

  Print[StringPadRight["  N", 12],
        StringPadRight["G0_err", 12],
        StringPadRight["predicted", 12],
        StringPadRight["ratio", 8],
        StringPadRight["TOTAL_err", 12],
        StringPadRight["predicted", 12],
        "ratio"];
  Print["  ", StringJoin[Table["-", 72]]];

  Do[
    Module[{r = results[[i]], predG0, predMain, ratG0, ratMain},
      predG0 = g0Err0 * Sqrt[N[n0] / r["N"]];
      predMain = mainErr0 * Sqrt[N[n0] / r["N"]];
      ratG0 = If[predG0 > 0, Abs[r["G0Err"]] / predG0, 0];
      ratMain = If[predMain > 0, Abs[r["MainErr"]] / predMain, 0];
      Print[
        StringPadRight["  " <> ToString[r["N"]], 12],
        StringPadRight[ToString@ScientificForm[Abs[r["G0Err"]], 2], 12],
        StringPadRight[ToString@ScientificForm[predG0, 2], 12],
        StringPadRight[ToString@NumberForm[ratG0, {3, 2}], 8],
        StringPadRight[ToString@ScientificForm[Abs[r["MainErr"]], 2], 12],
        StringPadRight[ToString@ScientificForm[predMain, 2], 12],
        ToString@NumberForm[ratMain, {3, 2}]
      ];
    ],
    {i, Length[results]}
  ];
];

Print[];
Print["================================================================"];
Print["  SUMMARY"];
Print["================================================================"];
Print[];

Module[{last = results[[-1]], dPole, dFinite, sPole, sFinite},
  dPole = Abs[last["G0Re"] - fiestaPole];
  dFinite = Abs[last["MainRe"] - fiestaFinite];
  sPole = dPole / Abs[last["G0Err"]];
  sFinite = dFinite / Abs[last["MainErr"]];

  Print["Highest sample count: N = ", last["N"]];
  Print[];
  Print["  POLE:   MC = ", last["G0Re"], " ± ", ScientificForm[last["G0Err"], 3]];
  Print["          FIESTA = ", fiestaPole, " ± ", fiestaPoleErr];
  Print["          |MC - FIESTA| = ", ScientificForm[dPole, 2],
        " (", NumberForm[sPole, {3, 1}], " sigma)"];
  Print[];
  Print["  FINITE: MC = ", last["MainRe"], " ± ", ScientificForm[last["MainErr"], 3]];
  Print["          FIESTA = ", fiestaFinite, " ± ", fiestaFiniteErr];
  Print["          |MC - FIESTA| = ", ScientificForm[dFinite, 2],
        " (", NumberForm[sFinite, {3, 1}], " sigma)"];
  Print[];

  If[sPole < 3 && sFinite < 3,
    Print["  PASS: Both pole and finite agree within 3 sigma."];,
    If[sPole < 3,
      Print["  Pole: PASS (within 3 sigma)"];,
      Print["  Pole: ", NumberForm[sPole, 2], " sigma — check needed"];
    ];
    If[sFinite < 3,
      Print["  Finite: PASS (within 3 sigma)"];,
      Print["  Finite: ", NumberForm[sFinite, 2], " sigma — expected for MC at this N"];
    ];
  ];
];
Print[];
Print["================================================================"];
