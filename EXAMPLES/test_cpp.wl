(* Test C++ compilation and MC execution on Linux *)
(* Based on Example 6: C++ code generation *)

SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];

Print["=== Test: C++ code generation and execution ==="];
Print["Integral[0,Inf] dx1 dx2 / (1 + x1^2 + x2^2)^3 = Pi/8"];
Print[];

Module[{poly, vars, spec, verts, fanData, dualVerts, simplices,
        allSectorData, cppFile, binary, cppResult, kinFile, resultFile,
        nSamples = 100000},

  poly = 1 + x[1]^2 + x[2]^2;
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {-3},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  verts    = PolytopeVertices[poly^(-3), vars];
  fanData  = ComputeDecomposition[verts, "ShowProgress" -> False];
  {dualVerts, simplices} = fanData;

  Print["Fan: ", Length[dualVerts], " rays, ", Length[simplices], " sectors"];

  allSectorData = Table[
    ProcessSector[spec, dualVerts, simplices[[s]], s],
    {s, Length[simplices]}
  ];
  allSectorData = Select[allSectorData, AssociationQ];

  (* Generate C++ *)
  Quiet[CreateDirectory[FileNameJoin[{Directory[], "INTERFILES"}]]];
  cppFile = FileNameJoin[{Directory[], "INTERFILES", "test_mc_generated.cpp"}];
  binary  = FileNameJoin[{Directory[], "INTERFILES", "test_mc"}];

  cppResult = GenerateCppMonteCarlo[
    allSectorData, {},
    spec, cppFile,
    "NSamples" -> nSamples
  ];

  Print["C++ generated: ", AssociationQ[cppResult]];

  (* Compile *)
  If[CompileCpp[cppFile, binary, False] =!= $Failed,
    Print["Compilation successful!"];

    (* Write dummy kinematic data (no params, just need one line) *)
    kinFile = FileNameJoin[{Directory[], "INTERFILES", "test_kin_data.txt"}];
    resultFile = FileNameJoin[{Directory[], "INTERFILES", "test_mc_results.txt"}];
    Export[kinFile, "0\n", "Text"];

    (* Run *)
    Module[{proc, result},
      proc = RunProcess[{binary, kinFile, resultFile,
                        ToString[nSamples], "4"}];
      Print["Exit code: ", proc["ExitCode"]];
      If[proc["ExitCode"] == 0,
        result = Import[resultFile, "Table"];
        Print["MC result: ", result];
        Print["Exact: Pi/8 = ", N[Pi/8]];
        If[Length[result] > 0 && Length[result[[1]]] >= 4,
          Print["MC value: ", result[[1, 1]], " +/- ", result[[1, 3]]];
          Print["Relative error: ",
            Abs[result[[1, 1]] - N[Pi/8]] / N[Pi/8] * 100, " %"];
        ];,
        Print["STDERR: ", proc["StandardError"]];
      ];
    ];,
    Print["Compilation FAILED"];
  ];
];
Print[];
Print["=== C++ test complete ==="];
