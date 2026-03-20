(* ::Package:: *)

(* ============================================================================
   tropical_eval_examples.wl

   Examples demonstrating the tropical_eval evaluation pipeline.
   Requires: tropical_fan.wl, tropical_eval.wl, Polymake, g++ with OpenMP.

   Load the package first:
     SetDirectory[NotebookDirectory[]];
     Get["tropical_eval.wl"];

   Or run as a script:
     wolframscript -file tropical_eval_examples.wl
   ============================================================================ *)

(* --- Load package --- *)
SetDirectory[NotebookDirectory[]];
Get["tropical_fan.wl"];
Get["tropical_eval.wl"];


(* ============================================================================
   Example 1: ProcessSector \[LongDash] basic 2D convergent integral

   Integral[0,Inf] dx1 dx2 / (1 + x1^2 + x2^2)^3

   This is the simplest case: one polynomial with positive coefficients,
   real negative exponent, no kinematic parameters, no epsilon regulator.
   ============================================================================ *)

Print["=== Example 1: Basic 2D convergent integral ==="];
Print["Integral[0,Inf] dx1 dx2 / (1 + x1^2 + x2^2)^3"];
Print[];

Module[{poly, vars, spec, verts, fanData, dualVerts, simplices},

  poly = 1 + x[1]^2 + x[2]^2;
  vars = {x[1], x[2]};

  (* IntegrandSpec: the standard input format for all pipeline functions *)
  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0},          (* no x_i^{A_i} prefactor *)
    "PolynomialExponents" -> {-3},            (* P^{-3} *)
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},               (* no parameters *)
    "RegulatorSymbol"    -> None              (* no epsilon *)
  |>;

  (* Step 1: Compute the tropical fan *)
  verts    = PolytopeVertices[poly^(-3), vars];
  fanData  = ComputeDecomposition[verts, "ShowProgress" -> False];
  {dualVerts, simplices} = fanData;

  Print["Fan: ", Length[dualVerts], " rays, ",
        Length[simplices], " sectors"];

  (* Step 2: Process each sector *)
  Do[
    Module[{sd},
      sd = ProcessSector[spec, dualVerts, simplices[[s]], s,
                         "Verbose" -> True];
      If[AssociationQ[sd],
        Print["  Effective exponents: ", sd["NewExponents"]];
        Print["  Prefactor: ", sd["Prefactor"]];
        Print["  # monomials in Q: ",
              Length /@ sd["FlattenedPolys"]];
        Print["  Dominant monomial (should have exps=0): ",
              SelectFirst[sd["FlattenedPolys"][[1]],
                AllTrue[#[[2]], TrueQ[# == 0] &] &]];
        Print[];
      ];
    ],
    {s, Length[simplices]}
  ];

  (* Step 3: Validate the decomposition against NIntegrate *)
  Print["Validating against NIntegrate..."];
  Module[{vr},
    vr = Quiet@ValidateDecomposition[spec, fanData, {}, 4];
    Print["  Direct NIntegrate:  ", vr["DirectResult"]];
    Print["  Sector sum:         ", vr["SectorSum"]];
    Print["  Relative error:     ", vr["RelativeError"]];
    Print["  Per-sector results: ", vr["SectorResults"]];
  ];
];
Print[];


(* ============================================================================
   Example 2: Complex exponents and the flattening check

   Integral[0,Inf] dx1 dx2 / (1 + 2 x1^2 + x2^2 + x1 x2^2 + 3 x1^2 x2)^{2+i}

   Demonstrates complex polynomial exponents. The polynomial always
   evaluates to a positive real (all coefficients positive), so
   P^{2+i} = exp((2+i) * log(P)) is well-defined.
   ============================================================================ *)

Print["=== Example 2: Complex exponents ==="];
Print["Integral[0,Inf] dx1 dx2 / (1+2x1^2+x2^2+x1*x2^2+3x1^2*x2)^{2+I}"];
Print[];

Module[{poly, vars, A, spec, verts, fanData, sectorData},

  poly = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2];
  vars = {x[1], x[2]};
  A    = 2 + I;

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {-A},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  verts   = PolytopeVertices[poly^(-Re[A]), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  (* Process sector 1 and spot-check the integrand magnitude *)
  sectorData = ProcessSector[spec, fanData[[1]], fanData[[2, 1]], 1];

  Print["Sector 1 details:"];
  Print["  Effective exponents (complex): ", sectorData["NewExponents"]];
  Print["  Prefactor (complex): ", sectorData["Prefactor"]];

  Print[];
  Print["Flattening magnitude check (should be O(1)):"];
  Module[{check},
    check = CheckFlatteningMagnitude[sectorData, 10];
    Print["  Mean |integrand|: ", check["Mean"]];
    Print["  Min:  ", check["Min"]];
    Print["  Max:  ", check["Max"]];
  ];

  Print[];
  Print["Full validation:"];
  Module[{vr},
    vr = Quiet@ValidateDecomposition[spec, fanData, {}, 3];
    Print["  Direct NIntegrate: ", vr["DirectResult"]];
    Print["  Sector sum:        ", vr["SectorSum"]];
    Print["  Relative error:    ", vr["RelativeError"]];
  ];
];
Print[];


(* ============================================================================
   Example 3: Kinematic-dependent integral

   Integral[0,Inf] dx1 dx2 / (1 + lam*x1^2 + x2^2 + x1*x2^2)^2

   The polynomial coefficient 'lam' is a kinematic parameter that varies.
   This demonstrates how to set up kinematic symbols and evaluate at
   multiple parameter values.
   ============================================================================ *)

Print["=== Example 3: Kinematic-dependent integral ==="];
Print["Integral[0,Inf] dx1 dx2 / (1 + lam*x1^2 + x2^2 + x1*x2^2)^2"];
Print[];

Module[{lam, poly, vars, spec, verts, fanData, lamValues},

  lam  = Symbol["lam"];
  poly = 1 + lam x[1]^2 + x[2]^2 + x[1] x[2]^2;
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {lam},
    "RegulatorSymbol"    -> None
  |>;

  (* Compute fan at lam=1 (topology doesn't change with lam) *)
  verts   = PolytopeVertices[(poly /. lam -> 1)^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Fan: ", Length[fanData[[2]]], " sectors"];

  (* Process one sector with symbolic lam to see the structure *)
  Module[{sd},
    sd = ProcessSector[spec, fanData[[1]], fanData[[2, 1]], 1];
    Print["Sector 1 (symbolic):"];
    Print["  Effective exponents: ", sd["NewExponents"]];
    Print["  Prefactor: ", sd["Prefactor"]];
    Print["  First monomial coeff: ", sd["FlattenedPolys"][[1, 1, 1]]];
    Print[];
  ];

  (* Evaluate at several values of lam *)
  lamValues = {0.5, 1.0, 2.0, 5.0};

  Print["Validation at multiple kinematic points:"];
  Do[
    Module[{kinRules, vr},
      kinRules = {lam -> lamVal};
      vr = Quiet@ValidateDecomposition[spec, fanData, kinRules, 3];
      If[AssociationQ[vr],
        Print["  lam = ", lamVal,
              ":  direct = ", vr["DirectResult"],
              "  sector sum = ", vr["SectorSum"],
              "  rel err = ", vr["RelativeError"]];
      ];
    ],
    {lamVal, lamValues}
  ];
];
Print[];


(* ============================================================================
   Example 4: Inspecting the tropical factoring

   This example shows what happens inside ProcessSector step by step:
   the monomial substitution, the tropical factoring (extracting the
   dominant monomial), and the resulting cleared polynomials.

   Polynomial: P = 1 + 2 x1^2 + x2^2 + x1 x2^2 + 3 x1^2 x2
   Rays: rho_1 = (0,1), rho_2 = (1,1)  (sector 5 of the fan)
   ============================================================================ *)

Print["=== Example 4: Tropical factoring step by step ==="];
Print[];

Module[{poly, vars, spec, dualVerts, simplex, sd},

  poly = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2];
  vars = {x[1], x[2]};

  dualVerts = {{0, 1}, {1, 1}};
  simplex   = {1, 2};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {-Symbol["A"]},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  sd = ProcessSector[spec, dualVerts, simplex, 1, "Verbose" -> True];

  Print[];
  Print["Step-by-step breakdown:"];

  Print["  Rays: rho1 = ", dualVerts[[1]], ", rho2 = ", dualVerts[[2]]];
  Print["  M = -Transpose[rays] = ", sd["RayMatrix"]];
  Print["  det(M) = ", sd["DetM"]];
  Print[];

  Print["  Raw a (from monomial exponents only):"];
  Print["    a = (monoExps+1) . M = ", sd["RawExponents"]];
  Print[];

  Print["  Transformed polynomial monomials (in y-space):"];
  Do[
    Print["    ", mono[[1]], " * y^", mono[[2]]],
    {mono, sd["TransformedPolys"][[1]]}
  ];
  Print[];

  Print["  Min exponents (dominant monomial): ", sd["MinExponents"]];
  Print[];

  Print["  Cleared polynomial Q (non-negative exponents):"];
  Do[
    Print["    ", mono[[1]], " * y^", mono[[2]]],
    {mono, sd["ClearedPolys"][[1]]}
  ];
  Print[];

  Print["  Effective exponents:"];
  Print["    a_eff = rawA + B*minExp = ", sd["NewExponents"]];
  Print["    (For A=2: a_eff = ", sd["NewExponents"] /. Symbol["A"] -> 2, ")"];
  Print["    (For A=3: a_eff = ", sd["NewExponents"] /. Symbol["A"] -> 3, ")"];
  Print[];

  Print["  Flattened polynomial (exponents / a_eff):"];
  Do[
    Print["    ", mono[[1]], " * y'^", mono[[2]]],
    {mono, sd["FlattenedPolys"][[1]]}
  ];
];
Print[];


(* ============================================================================
   Example 5: MmaToC \[LongDash] Mathematica to C++ conversion

   Demonstrates the expression converter used by the C++ code generator.
   ============================================================================ *)

Print["=== Example 5: MmaToC expression converter ==="];
Print[];

Module[{lam, paramMap},

  lam = Symbol["lam"];
  paramMap = <|lam -> "params[0]"|>;

  Print["Integer:      ", MmaToC[42]];
  Print["Rational:     ", MmaToC[3/7]];
  Print["Real:         ", MmaToC[3.14159]];
  Print["Complex:      ", MmaToC[2 + 3 I]];
  Print["Power:        ", MmaToC[Power[x, 3]]];
  Print["Sqrt:         ", MmaToC[Sqrt[x]]];
  Print["Log:          ", MmaToC[Log[x]]];
  Print["Exp:          ", MmaToC[Exp[x]]];
  Print["Sum:          ", MmaToC[a + b + c]];
  Print["Product:      ", MmaToC[a * b * c]];
  Print["With params:  ", MmaToC[1 + 2 lam, paramMap]];
  Print["Complex expr: ", MmaToC[Exp[(2 + I) * Log[1 + lam]], paramMap]];
];
Print[];


(* ============================================================================
   Example 6: C++ code generation (convergent, no kinematics)

   Generates a C++ Monte Carlo source file for a simple convergent integral,
   compiles it, and runs it.

   Integral[0,Inf] dx1 dx2 / (1 + x1^2 + x2^2)^3  =  Pi/8
   ============================================================================ *)

Print["=== Example 6: C++ code generation and execution ==="];
Print["Integral[0,Inf] dx1 dx2 / (1 + x1^2 + x2^2)^3 = Pi/8"];
Print[];

Module[{poly, vars, spec, verts, fanData, allSectors,
        convergent, cppResult, cppFile, exact},

  poly = 1 + x[1]^2 + x[2]^2;
  vars = {x[1], x[2]};
  exact = Pi/8;

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {-3},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  verts   = PolytopeVertices[poly^(-3), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  (* Process all sectors *)
  allSectors = Table[
    ProcessSector[spec, fanData[[1]], fanData[[2, s]], s],
    {s, Length[fanData[[2]]]}
  ];

  convergent = Select[allSectors,
    AssociationQ[#] && !#["IsDivergent"] &];

  Print["Processed ", Length[convergent], " convergent sectors"];

  (* Generate C++ *)
  cppFile = FileNameJoin[{Directory[], "example6_mc.cpp"}];
  cppResult = GenerateCppMonteCarlo[
    convergent, {},   (* no divergent sectors *)
    spec, cppFile,
    "NSamples" -> 500000
  ];

  If[AssociationQ[cppResult],
    Print[];

    (* Compile *)
    Module[{binary, kinFile, resultFile, runResult},
      binary     = FileNameJoin[{Directory[], "example6_mc"}];
      kinFile    = FileNameJoin[{Directory[], "example6_kin.txt"}];
      resultFile = FileNameJoin[{Directory[], "example6_results.txt"}];

      (* No kinematics: write count = 1 *)
      Export[kinFile, "1\n", "Text"];

      If[CompileCpp[cppFile, binary, False] =!= $Failed,
        Print[];
        runResult = RunProcess[{binary, kinFile, resultFile,
                                "500000", "2"}];
        If[runResult["ExitCode"] == 0,
          Module[{output, vals, mcVal, mcErr, relErr},
            output = Import[resultFile, "Text"];
            vals   = ToExpression /@ StringSplit[StringTrim[output]];
            mcVal  = vals[[1]] + I vals[[2]];
            mcErr  = vals[[3]] + I vals[[4]];

            Print["Results:"];
            Print["  Exact   = ", N[exact]];
            Print["  MC      = ", mcVal, " +/- ", Abs[mcErr]];
            Print["  Rel err = ",
                  Abs[(Re[mcVal] - N[exact]) / N[exact]]];
          ];,
          Print["Run failed: ", runResult["StandardError"]];
        ];

        (* Clean up *)
        Quiet[DeleteFile /@ {cppFile, binary, kinFile, resultFile}];
      ];
    ];
  ];
];
Print[];


(* ============================================================================
   Example 7: C++ code generation with kinematic parameters

   Generates C++ code for a kinematic-dependent integral and evaluates
   at 10 values of lambda.

   Integral[0,Inf] dx1 dx2 / (1 + lam*x1^2 + x2^2 + x1*x2^2)^2
   ============================================================================ *)

Print["=== Example 7: C++ with kinematic scan ==="];
Print["Integral[0,Inf] dx1 dx2 / (1+lam*x1^2+x2^2+x1*x2^2)^2"];
Print["10 values of lam in [0.5, 5]"];
Print[];

Module[{lam, poly, vars, spec, verts, fanData, allSectors,
        convergent, cppResult, cppFile,
        lamValues, kinPoints},

  lam  = Symbol["lam"];
  poly = 1 + lam x[1]^2 + x[2]^2 + x[1] x[2]^2;
  vars = {x[1], x[2]};

  lamValues = Table[0.5 + 4.5 (i - 1)/9, {i, 10}];
  kinPoints = List /@ lamValues;

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {lam},
    "RegulatorSymbol"    -> None
  |>;

  verts   = PolytopeVertices[(poly /. lam -> 1)^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  allSectors = Table[
    ProcessSector[spec, fanData[[1]], fanData[[2, s]], s],
    {s, Length[fanData[[2]]]}
  ];

  convergent = Select[allSectors,
    AssociationQ[#] && !#["IsDivergent"] &];

  cppFile = FileNameJoin[{Directory[], "example7_mc.cpp"}];
  cppResult = GenerateCppMonteCarlo[
    convergent, {}, spec, cppFile,
    "NSamples" -> 200000
  ];

  If[AssociationQ[cppResult],
    Module[{binary, kinFile, resultFile, runResult},
      binary     = FileNameJoin[{Directory[], "example7_mc"}];
      kinFile    = FileNameJoin[{Directory[], "example7_kin.txt"}];
      resultFile = FileNameJoin[{Directory[], "example7_results.txt"}];

      (* Write kinematic data: one lam value per line *)
      Export[kinFile,
        StringRiffle[ToString[CForm[#]] & /@ lamValues, "\n"] <> "\n",
        "Text"
      ];

      If[CompileCpp[cppFile, binary, False] =!= $Failed,
        Print[];
        runResult = RunProcess[{binary, kinFile, resultFile,
                                "200000", "2"}];
        If[runResult["ExitCode"] == 0,
          Module[{lines, parsed},
            lines  = Select[
              StringSplit[Import[resultFile, "Text"], "\n"],
              StringLength[StringTrim[#]] > 0 &
            ];
            parsed = (ToExpression /@ StringSplit[#]) & /@ lines;

            Print["Results (MC vs NIntegrate):"];
            Print[StringPadRight["  lam", 10],
                  StringPadRight["MC Re", 16],
                  StringPadRight["MC err", 14],
                  "NIntegrate"];
            Do[
              Module[{mcRe, mcErr, niResult},
                mcRe  = parsed[[i, 1]];
                mcErr = parsed[[i, 3]];
                niResult = Quiet@NIntegrate[
                  1 / (1 + lamValues[[i]] t1^2 + t2^2 + t1 t2^2)^2,
                  {t1, 0, Infinity}, {t2, 0, Infinity},
                  PrecisionGoal -> 5
                ];
                Print["  ",
                  StringPadRight[ToString@NumberForm[lamValues[[i]], {4,2}], 8],
                  StringPadRight[ToString@NumberForm[mcRe, {8,6}], 16],
                  StringPadRight[ToString[mcErr], 14],
                  ToString@NumberForm[niResult, {8,6}]];
              ],
              {i, Length[lamValues]}
            ];
          ];,
          Print["Run failed: ", runResult["StandardError"]];
        ];

        Quiet[DeleteFile /@ {cppFile, binary, kinFile, resultFile}];
      ];
    ];
  ];
];
Print[];


(* ============================================================================
   Example 8: ParsePolynomial utility

   Shows how ParsePolynomial extracts {coefficient, exponentVector} pairs
   from a symbolic polynomial.
   ============================================================================ *)

Print["=== Example 8: ParsePolynomial utility ==="];
Print[];

Module[{lam, poly, vars, parsed},

  lam  = Symbol["lam"];
  poly = 1 + 3 lam x[1]^2 + x[2]^3 + 2 x[1] x[2];
  vars = {x[1], x[2]};

  parsed = ParsePolynomial[poly, vars];

  Print["Polynomial: ", poly];
  Print["Variables:  ", vars];
  Print[];
  Print["Parsed monomials {coefficient, exponents}:"];
  Do[
    Print["  ", mono[[1]], " * x^", mono[[2]]],
    {mono, parsed}
  ];
];
Print[];


(* ============================================================================
   Example 9: Multiple polynomials (product of two factors)

   Integral[0,Inf] dx1 dx2 x1^{1/2} / [(1+x1+x2)^2 * (1+x1*x2)]

   Shows how to handle a product of polynomials with separate exponents,
   plus a non-trivial monomial prefactor x1^{1/2}.
   The two factors have different Newton polytopes, ensuring a
   full-dimensional combined polytope.
   ============================================================================ *)

Print["=== Example 9: Multiple polynomials ==="];
Print["Integral[0,Inf] dx1 dx2 x1^{1/2} / [(1+x1+x2)^2 * (1+x1*x2)]"];
Print[];

Module[{p1, p2, vars, spec, integrand, verts, fanData, vr},

  p1   = 1 + x[1] + x[2];
  p2   = 1 + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {p1, p2},
    "MonomialExponents"  -> {1/2, 0},        (* x1^{1/2} prefactor *)
    "PolynomialExponents" -> {-2, -1},        (* P1^{-2} * P2^{-1} *)
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  (* Combined integrand for fan computation *)
  integrand = p1^(-2) * p2^(-1);
  verts     = PolytopeVertices[integrand, vars];
  fanData   = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Fan: ", Length[fanData[[2]]], " sectors"];

  (* Process sectors *)
  Do[
    Module[{sd},
      sd = ProcessSector[spec, fanData[[1]], fanData[[2, s]], s];
      If[AssociationQ[sd],
        Print["Sector ", s, ": effA = ", sd["NewExponents"],
              ", div = ", sd["IsDivergent"]];
      ];
    ],
    {s, Length[fanData[[2]]]}
  ];

  Print[];
  Print["Validation:"];
  vr = Quiet@ValidateDecomposition[spec, fanData, {}, 3];
  If[AssociationQ[vr],
    Print["  Direct NIntegrate: ", vr["DirectResult"]];
    Print["  Sector sum:        ", vr["SectorSum"]];
    Print["  Relative error:    ", vr["RelativeError"]];
  ];
];
Print[];


Print["=== All examples complete ==="];
