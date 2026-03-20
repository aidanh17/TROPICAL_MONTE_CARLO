(* ============================================================================
   Convergence Study: Tropical MC -> FIESTA as N_samples increases

   Uses the already-compiled C++ binaries for bubblepp presector 1
   at the equilateral point k1=k2=k3=1.

   FIESTA reference values (from fiesta_bubblepp_standalone.wl):
     1/ep pole:  0.000513 ± 1e-6
     ep^0 finite: 0.024968 ± 3.7e-5

   We run the MC at N = 50K, 100K, 200K, 500K, 1M, 2M, 5M, 10M
   and track convergence of both the pole (G0) and finite (TOTAL).
   ============================================================================ *)

Print["================================================================"];
Print["  Convergence Study: Tropical MC → FIESTA"];
Print["  Bubblepp presector 1, k1=k2=k3=1"];
Print["================================================================"];
Print[];

testDir = "/home/aidanh/Desktop/Tropical_Monte_Carlo_Final/Bubble1final/tropical_code_claudef/FIESTA_Compare/bubblepp_test";

mainBin = FileNameJoin[{testDir, "main_mc"}];
g0Bin   = FileNameJoin[{testDir, "g0_mc"}];
kinFile = FileNameJoin[{testDir, "kin.txt"}];

(* FIESTA reference *)
fiestaPole   = 0.000513;
fiestaFinite = 0.024968;

Print["FIESTA reference:"];
Print["  Pole (1/ep):  ", fiestaPole, " ± 1e-6"];
Print["  Finite (ep^0): ", fiestaFinite, " ± 3.7e-5"];
Print[];

(* Sample counts to test *)
sampleCounts = {50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000};

nThreads = 4;

(* Collect results *)
results = {};

Do[
  Module[{mainResultFile, g0ResultFile, runMain, runG0,
          mainVal, g0Val, mainErr, g0Err,
          mainOutput, g0Output, mainVals, g0Vals,
          t0},

    mainResultFile = FileNameJoin[{testDir, "conv_main_" <> ToString[nSamp] <> ".txt"}];
    g0ResultFile   = FileNameJoin[{testDir, "conv_g0_" <> ToString[nSamp] <> ".txt"}];

    t0 = AbsoluteTime[];

    (* Run main and G0 in parallel *)
    Module[{procMain, procG0},
      procMain = StartProcess[{mainBin, kinFile, mainResultFile,
                               ToString[nSamp], ToString[nThreads]}];
      procG0   = StartProcess[{g0Bin, kinFile, g0ResultFile,
                               ToString[nSamp], ToString[nThreads]}];

      While[ProcessStatus[procMain] === "Running" ||
            ProcessStatus[procG0] === "Running",
        Pause[0.5]];
    ];

    Module[{elapsed = AbsoluteTime[] - t0},
      (* Parse main result *)
      mainOutput = Import[mainResultFile, "Text"];
      mainVals = ToExpression /@ StringSplit[StringTrim[mainOutput]];
      mainVal = mainVals[[1]];
      mainErr = mainVals[[3]];

      (* Parse G0 result *)
      g0Output = Import[g0ResultFile, "Text"];
      g0Vals = ToExpression /@ StringSplit[StringTrim[g0Output]];
      g0Val = g0Vals[[1]];
      g0Err = g0Vals[[3]];

      AppendTo[results, <|
        "N" -> nSamp,
        "MainVal" -> mainVal, "MainErr" -> mainErr,
        "G0Val" -> g0Val, "G0Err" -> g0Err,
        "Time" -> elapsed
      |>];

      (* Clean up *)
      Quiet[DeleteFile /@ {mainResultFile, g0ResultFile}];
    ];
  ],
  {nSamp, sampleCounts}
];

(* ============================================================ *)
(* Print results table                                           *)
(* ============================================================ *)

Print["================================================================"];
Print["  G0 (1/eps pole coefficient) Convergence"];
Print["================================================================"];
Print[];
Print[StringPadRight["  N_samples", 14],
      StringPadRight["G0 (MC)", 16],
      StringPadRight["MC error", 14],
      StringPadRight["FIESTA ref", 14],
      StringPadRight["|diff|/err", 12],
      "Time(s)"];
Print["  ", StringJoin[Table["-", 80]]];

Do[
  Module[{r = results[[i]], diff, sigmas},
    diff = Abs[r["G0Val"] - fiestaPole];
    sigmas = If[r["G0Err"] > 0, diff / r["G0Err"], Infinity];
    Print[
      StringPadRight["  " <> ToString[r["N"]], 14],
      StringPadRight[ToString@NumberForm[r["G0Val"], {8, 7}], 16],
      StringPadRight[ToString@ScientificForm[r["G0Err"], 3], 14],
      StringPadRight[ToString[fiestaPole], 14],
      StringPadRight[ToString@NumberForm[sigmas, {4, 1}], 12],
      ToString@NumberForm[r["Time"], {4, 1}]
    ];
  ],
  {i, Length[results]}
];
Print[];

Print["================================================================"];
Print["  TOTAL (finite part, eps^0) Convergence"];
Print["================================================================"];
Print[];
Print[StringPadRight["  N_samples", 14],
      StringPadRight["TOTAL (MC)", 16],
      StringPadRight["MC error", 14],
      StringPadRight["FIESTA ref", 14],
      StringPadRight["|diff|/err", 12],
      "Time(s)"];
Print["  ", StringJoin[Table["-", 80]]];

Do[
  Module[{r = results[[i]], diff, sigmas},
    diff = Abs[r["MainVal"] - fiestaFinite];
    sigmas = If[r["MainErr"] > 0, diff / r["MainErr"], Infinity];
    Print[
      StringPadRight["  " <> ToString[r["N"]], 14],
      StringPadRight[ToString@NumberForm[r["MainVal"], {8, 7}], 16],
      StringPadRight[ToString@ScientificForm[r["MainErr"], 3], 14],
      StringPadRight[ToString[fiestaFinite], 14],
      StringPadRight[ToString@NumberForm[sigmas, {4, 1}], 12],
      ToString@NumberForm[r["Time"], {4, 1}]
    ];
  ],
  {i, Length[results]}
];
Print[];

(* ============================================================ *)
(* Verify 1/sqrt(N) scaling                                      *)
(* ============================================================ *)

Print["================================================================"];
Print["  Error scaling verification (should be ~ 1/sqrt(N))"];
Print["================================================================"];
Print[];

Module[{firstG0Err, firstN, firstMainErr},
  firstN = results[[1]]["N"];
  firstG0Err = Abs[results[[1]]["G0Err"]];
  firstMainErr = Abs[results[[1]]["MainErr"]];

  Print[StringPadRight["  N_samples", 14],
        StringPadRight["G0 err", 14],
        StringPadRight["predicted", 14],
        StringPadRight["TOTAL err", 14],
        "predicted"];
  Print["  ", StringJoin[Table["-", 70]]];

  Do[
    Module[{r = results[[i]], predicted, predictedMain, ratio},
      predicted = firstG0Err * Sqrt[firstN / r["N"]];
      predictedMain = firstMainErr * Sqrt[firstN / r["N"]];
      Print[
        StringPadRight["  " <> ToString[r["N"]], 14],
        StringPadRight[ToString@ScientificForm[Abs[r["G0Err"]], 3], 14],
        StringPadRight[ToString@ScientificForm[predicted, 3], 14],
        StringPadRight[ToString@ScientificForm[Abs[r["MainErr"]], 3], 14],
        ToString@ScientificForm[predictedMain, 3]
      ];
    ],
    {i, Length[results]}
  ];
];
Print[];

(* ============================================================ *)
(* Summary                                                       *)
(* ============================================================ *)

Print["================================================================"];
Print["  SUMMARY"];
Print["================================================================"];
Print[];

Module[{lastR = results[[-1]], lastDiffPole, lastDiffFinite,
        lastSigmaPole, lastSigmaFinite},
  lastDiffPole = Abs[lastR["G0Val"] - fiestaPole];
  lastDiffFinite = Abs[lastR["MainVal"] - fiestaFinite];
  lastSigmaPole = If[lastR["G0Err"] > 0, lastDiffPole / lastR["G0Err"], 0];
  lastSigmaFinite = If[lastR["MainErr"] > 0, lastDiffFinite / lastR["MainErr"], 0];

  Print["At ", lastR["N"], " samples:"];
  Print["  G0 pole:  MC = ", lastR["G0Val"], " ± ", lastR["G0Err"]];
  Print["            FIESTA = ", fiestaPole];
  Print["            |diff| = ", lastDiffPole, " = ", NumberForm[lastSigmaPole, 2], " sigma"];
  Print[];
  Print["  Finite:   MC = ", lastR["MainVal"], " ± ", lastR["MainErr"]];
  Print["            FIESTA = ", fiestaFinite];
  Print["            |diff| = ", lastDiffFinite, " = ", NumberForm[lastSigmaFinite, 2], " sigma"];
  Print[];

  If[lastSigmaPole < 3 && lastSigmaFinite < 3,
    Print["  ✓ Both pole and finite agree within 3 sigma."];,
    Print["  Note: larger-than-3-sigma difference (may need more samples)."];
  ];
];
Print[];
Print["================================================================"];
