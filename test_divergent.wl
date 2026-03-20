(* Run only the 1/epsilon divergent pipeline tests: Tests 3, 4, 8-18 *)

SetDirectory[DirectoryName[$InputFileName]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];

Print["Package loaded."];
Print[];

(* Load the example file to get the test functions *)
Get[FileNameJoin[{Directory[], "tropical_eval_examples.wl"}]];

Print[];
Print["================================================================"];
Print["  1/epsilon Divergent Pipeline Tests"];
Print["================================================================"];
Print[];

results = {};
nPass = 0;
nFail = 0;

tests = {
  {3,  RunTest3},
  {4,  RunTest4},
  {8,  RunTest8},
  {9,  RunTest9},
  {10, RunTest10},
  {11, RunTest11},
  {12, RunTest12},
  {13, RunTest13},
  {14, RunTest14},
  {15, RunTest15},
  {16, RunTest16},
  {17, RunTest17},
  {18, RunTest18}
};

Do[
  Module[{num, fn, pass},
    {num, fn} = t;
    Print[];
    pass = fn[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test " <> ToString[num], pass}];
    Print[];
  ],
  {t, tests}
];

Print[];
Print["================================================================"];
Print["  DIVERGENT PIPELINE SUMMARY"];
Print["================================================================"];
Print[];
Do[
  Print["  ", r[[1]], ": ", If[r[[2]], "PASS", "FAIL"]],
  {r, results}
];
Print[];
Print["  Passed: ", nPass, " / ", nPass + nFail];
If[nFail > 0,
  Print["  FAILED: ", nFail, " test(s)"];,
  Print["  All tests passed!"];
];
Print[];
