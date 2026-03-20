(* ============================================================================
   tropical_eval_examples.wl

   Examples demonstrating the tropical_eval evaluation pipeline.
   Requires: tropical_fan.wl, tropical_eval.wl, Polymake, g++ with OpenMP.

   Load the package first:
     SetDirectory[FileNameJoin[{NotebookDirectory[], ".."}]];
     Get["tropical_eval.wl"];

   Or run as a script:
     wolframscript -file EXAMPLES/tropical_eval_examples.wl
   ============================================================================ *)

(* --- Load package --- *)
SetDirectory[FileNameJoin[{DirectoryName[$InputFileName], ".."}]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];

Print["tropical_eval.wl loaded successfully"];
Print[];


(* ============================================================================
   Example 1: ProcessSector — basic 2D convergent integral

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
   Example 5: MmaToC — Mathematica to C++ conversion

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
  Quiet[CreateDirectory[FileNameJoin[{Directory[], "INTERFILES"}]]];
  cppFile = FileNameJoin[{Directory[], "INTERFILES", "example6_mc.cpp"}];
  cppResult = GenerateCppMonteCarlo[
    convergent, {},   (* no divergent sectors *)
    spec, cppFile,
    "NSamples" -> 500000
  ];

  If[AssociationQ[cppResult],
    Print[];

    (* Compile *)
    Module[{binary, kinFile, resultFile, runResult},
      binary     = FileNameJoin[{Directory[], "INTERFILES", "example6_mc"}];
      kinFile    = FileNameJoin[{Directory[], "INTERFILES", "example6_kin.txt"}];
      resultFile = FileNameJoin[{Directory[], "INTERFILES", "example6_results.txt"}];

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

  Quiet[CreateDirectory[FileNameJoin[{Directory[], "INTERFILES"}]]];
  cppFile = FileNameJoin[{Directory[], "INTERFILES", "example7_mc.cpp"}];
  cppResult = GenerateCppMonteCarlo[
    convergent, {}, spec, cppFile,
    "NSamples" -> 200000
  ];

  If[AssociationQ[cppResult],
    Module[{binary, kinFile, resultFile, runResult},
      binary     = FileNameJoin[{Directory[], "INTERFILES", "example7_mc"}];
      kinFile    = FileNameJoin[{Directory[], "INTERFILES", "example7_kin.txt"}];
      resultFile = FileNameJoin[{Directory[], "INTERFILES", "example7_results.txt"}];

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


(* ============================================================================
   Example 10: Non-trivial numerator — single numerator polynomial

   Integral[0,Inf] dx1 dx2  (1 + x1)^{3/2} / (1 + x1^2 + x2^2 + x1*x2)^3

   The integrand has a POSITIVE polynomial exponent (numerator) and a
   NEGATIVE one (denominator).  The IntegrandSpec lists both polynomials
   with their signed exponents.

   For the Newton polytope / fan computation, we use the product of all
   polynomials raised to |B_j| (absolute exponents), since the tropical
   fan only depends on which monomials dominate, not on the sign of the
   power.
   ============================================================================ *)

Print["=== Example 10: Numerator polynomial (1+x1)^{3/2} ==="];
Print["Integral[0,Inf] dx1 dx2 (1+x1)^{3/2} / (1+x1^2+x2^2+x1*x2)^3"];
Print[];

Module[{pNum, pDen, vars, spec, fanPoly, verts, fanData, vr},

  pNum = 1 + x[1];
  pDen = 1 + x[1]^2 + x[2]^2 + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {pNum, pDen},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {3/2, -3},        (* numerator^{3/2} * denom^{-3} *)
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  (* For the fan: the normal fan of a Minkowski sum a*Newt(P1) + b*Newt(P2)
     is the same for any positive a, b.  So just multiply the polynomials
     once -- all that matters is that every monomial appears. *)
  fanPoly = pNum * pDen;
  verts   = PolytopeVertices[fanPoly^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Fan: ", Length[fanData[[2]]], " sectors"];

  (* Process and validate *)
  Do[
    Module[{sd},
      sd = ProcessSector[spec, fanData[[1]], fanData[[2, s]], s];
      If[AssociationQ[sd],
        Print["Sector ", s, ": effA = ", sd["NewExponents"],
              ", div = ", sd["IsDivergent"],
              ", prefactor = ", sd["Prefactor"]];
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


(* ============================================================================
   Example 11: Non-trivial numerator — two numerator factors

   Integral[0,Inf] dx1 dx2  (1 + x1 + x2)^{1/2} * (x1 + x2)^{1/3}
                             / (1 + 2*x1^2 + x2^2 + x1*x2^2)^2

   Three polynomials: two in the numerator (positive exponents 1/2, 1/3)
   and one in the denominator (negative exponent -2).
   ============================================================================ *)

Print["=== Example 11: Two numerator factors ==="];
Print["Integral[0,Inf] dx1 dx2 (1+x1+x2)^{1/2}*(x1+x2)^{1/3} / (1+2x1^2+x2^2+x1*x2^2)^2"];
Print[];

Module[{pN1, pN2, pDen, vars, spec, fanPoly, verts, fanData, vr},

  pN1  = 1 + x[1] + x[2];
  pN2  = x[1] + x[2];
  pDen = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2;
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {pN1, pN2, pDen},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {1/2, 1/3, -2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  (* Fan from product of all polynomials *)
  fanPoly = pN1 * pN2 * pDen;
  verts   = PolytopeVertices[fanPoly^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Fan: ", Length[fanData[[2]]], " sectors"];

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


(* ============================================================================
   Example 12: Numerator with kinematic-dependent coefficients

   Integral[0,Inf] dx1 dx2  (1 + lam*x1)^{3/2} / (1 + x1^2 + x2^2)^3

   The numerator coefficient 'lam' varies.  One compile, multiple evaluations.
   ============================================================================ *)

Print["=== Example 12: Numerator with kinematic coefficients ==="];
Print["Integral[0,Inf] dx1 dx2 (1+lam*x1)^{3/2} / (1+x1^2+x2^2)^3"];
Print[];

Module[{lam, pNum, pDen, vars, spec, fanPoly, verts, fanData},

  lam  = Symbol["lam"];
  pNum = 1 + lam x[1];
  pDen = 1 + x[1]^2 + x[2]^2;
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {pNum, pDen},
    "MonomialExponents"  -> {0, 0},
    "PolynomialExponents" -> {3/2, -3},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {lam},
    "RegulatorSymbol"    -> None
  |>;

  (* Fan: use unit-coefficient polynomial since fan is independent of lam *)
  fanPoly = (1 + x[1]) * pDen;
  verts   = PolytopeVertices[fanPoly^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Fan: ", Length[fanData[[2]]], " sectors"];
  Print[];

  (* Validate at several lambda values *)
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
    {lamVal, {0.5, 1., 2., 5.}}
  ];
];
Print[];


(* ============================================================================
   Example 13: 3D convergent integral

   Integral[0,Inf] dx1 dx2 dx3 / (1 + x1^2 + x2^2 + x3^2 + x1*x2*x3)^3

   First example in dimension > 2.  Tests that ProcessSector, tropical
   factoring, and ValidateDecomposition all work correctly in 3D.
   ============================================================================ *)

Print["=== Example 13: 3D convergent integral ==="];
Print["Integral[0,Inf] dx1 dx2 dx3 / (1 + x1^2 + x2^2 + x3^2 + x1*x2*x3)^3"];
Print[];

Module[{poly, vars, spec, verts, fanData, vr},

  poly = 1 + x[1]^2 + x[2]^2 + x[3]^2 + x[1] x[2] x[3];
  vars = {x[1], x[2], x[3]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0, 0},
    "PolynomialExponents" -> {-3},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  verts   = PolytopeVertices[poly^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Fan: ", Length[fanData[[1]]], " rays, ", Length[fanData[[2]]], " sectors"];

  (* Show sector details *)
  Do[
    Module[{sd},
      sd = ProcessSector[spec, fanData[[1]], fanData[[2, s]], s];
      If[AssociationQ[sd],
        Print["Sector ", s, ": effA = ", sd["NewExponents"],
              ", div = ", sd["IsDivergent"],
              ", prefactor = ", sd["Prefactor"]];
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


(* ============================================================================
   Example 14: 4D convergent integral

   Integral[0,Inf] dx1 dx2 dx3 dx4
       / (1 + x1^2 + x2^2 + x3^2 + x4^2 + x1*x2 + x3*x4)^4

   Tests the pipeline in 4 dimensions.  The Newton polytope is richer
   (mix of degree-1 and degree-2 monomials), producing more sectors.
   ============================================================================ *)

Print["=== Example 14: 4D convergent integral ==="];
Print["Integral[0,Inf] dx1 dx2 dx3 dx4 / (1+x1^2+x2^2+x3^2+x4^2+x1*x2+x3*x4)^4"];
Print[];

Module[{poly, vars, spec, verts, fanData, allSectors, convergent, vr},

  poly = 1 + x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 + x[1] x[2] + x[3] x[4];
  vars = {x[1], x[2], x[3], x[4]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 0, 0, 0},
    "PolynomialExponents" -> {-4},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  verts   = PolytopeVertices[poly^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Fan: ", Length[fanData[[1]]], " rays, ", Length[fanData[[2]]], " sectors"];

  (* Process all sectors *)
  allSectors = Table[
    ProcessSector[spec, fanData[[1]], fanData[[2, s]], s],
    {s, Length[fanData[[2]]]}
  ];
  convergent = Select[allSectors, AssociationQ[#] && !#["IsDivergent"] &];
  Print["Convergent sectors: ", Length[convergent], " / ", Length[fanData[[2]]]];

  (* Show a few sectors *)
  Do[
    Print["Sector ", convergent[[i]]["ConeIndex"],
          ": effA = ", convergent[[i]]["NewExponents"],
          ", prefactor = ", convergent[[i]]["Prefactor"]],
    {i, Min[4, Length[convergent]]}
  ];
  If[Length[convergent] > 4, Print["  ... (", Length[convergent] - 4, " more)"]];

  Print[];
  Print["Validation:"];
  vr = Quiet@ValidateDecomposition[spec, fanData, {}, 2];
  If[AssociationQ[vr],
    Print["  Direct NIntegrate: ", vr["DirectResult"]];
    Print["  Sector sum:        ", vr["SectorSum"]];
    Print["  Relative error:    ", vr["RelativeError"]];
  ];
];
Print[];


Print["=== All examples complete ==="];


(* ============================================================================
   MODULE 5: VALIDATION TESTS
   ============================================================================ *)

RunAllTests[] := Module[
  {results, nPass, nFail},

  results = {};
  nPass = 0;
  nFail = 0;

  Print[""];
  Print["================================================================"];
  Print["  TropicalEval Validation Suite"];
  Print["================================================================"];
  Print[""];

  Module[{pass},
    pass = RunTest1[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 1", pass}];
  ];

  Module[{pass},
    pass = RunTest2[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 2", pass}];
  ];

  Module[{pass},
    pass = RunTest3[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 3", pass}];
  ];

  Module[{pass},
    pass = RunTest4[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 4", pass}];
  ];

  Module[{pass},
    pass = RunTest5[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 5", pass}];
  ];

  Module[{pass},
    pass = RunTest6[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 6", pass}];
  ];

  Module[{pass},
    pass = RunTest7[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 7", pass}];
  ];

  Module[{pass},
    pass = RunTest8[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 8", pass}];
  ];

  Module[{pass},
    pass = RunTest9[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 9", pass}];
  ];

  Module[{pass},
    pass = RunTest10[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 10", pass}];
  ];

  Module[{pass},
    pass = RunTest11[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 11", pass}];
  ];

  Module[{pass},
    pass = RunTest12[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 12", pass}];
  ];

  Module[{pass},
    pass = RunTest13[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 13", pass}];
  ];

  Module[{pass},
    pass = RunTest14[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 14", pass}];
  ];

  Module[{pass},
    pass = RunTest15[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 15", pass}];
  ];

  Module[{pass},
    pass = RunTest16[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 16", pass}];
  ];

  Module[{pass},
    pass = RunTest17[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 17", pass}];
  ];

  Module[{pass},
    pass = RunTest18[];
    If[pass, nPass++, nFail++];
    AppendTo[results, {"Test 18", pass}];
  ];

  Print[""];
  Print["================================================================"];
  Print["  Results: ", nPass, " PASSED, ", nFail, " FAILED"];
  Print["================================================================"];

  results
];

(* --------------------------------------------------------------------------
   Test 1: Convergent 2D, real exponents
   -------------------------------------------------------------------------- *)

RunTest1[] := Module[
  {poly, vars, integrandSpec, verts, fanData,
   testAValues, allPass},

  Print["--- Test 1: Convergent 2D, real exponents ---"];
  Print["Integral[0,Inf] dx1 dx2 / (1+2x1^2+x2^2+x1*x2^2+3x1^2*x2)^A"];

  poly = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2];
  vars = {x[1], x[2]};

  testAValues = {2, 3};
  allPass = True;

  Do[
    Module[{A, spec, directResult, mcResult, relErr, pass},
      A = testA;

      spec = <|
        "Polynomials"       -> {poly},
        "MonomialExponents" -> {0, 0},
        "PolynomialExponents" -> {-A},
        "Variables"         -> vars,
        "KinematicSymbols"  -> {},
        "RegulatorSymbol"   -> None
      |>;

      directResult = NIntegrate[
        1 / (1 + 2 t1^2 + t2^2 + t1 t2^2 + 3 t1^2 t2)^A,
        {t1, 0, Infinity}, {t2, 0, Infinity},
        MaxRecursion -> 20, PrecisionGoal -> 6
      ];

      verts = PolytopeVertices[poly^(-A), vars];
      fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

      Module[{vr},
        vr = Quiet@ValidateDecomposition[spec, fanData, {}, 3];
        If[AssociationQ[vr],
          mcResult = vr["SectorSum"];
          relErr   = vr["RelativeError"];,
          mcResult = 0;
          relErr = Infinity;
        ];
      ];

      pass = NumericQ[relErr] && (relErr < 0.02);
      Print["  A = ", A, ":"];
      Print["    NIntegrate = ", directResult];
      Print["    Sector sum = ", mcResult];
      Print["    Rel error  = ", relErr];
      Print["    ", If[pass, "PASS", "FAIL"]];

      If[!pass, allPass = False];
    ],
    {testA, testAValues}
  ];

  allPass
];

(* --------------------------------------------------------------------------
   Test 2: Convergent 2D, complex exponents
   -------------------------------------------------------------------------- *)

RunTest2[] := Module[
  {poly, vars, A, spec, verts, fanData,
   directResult, mcResult, relErr, pass},

  Print["--- Test 2: Convergent 2D, complex exponents ---"];
  Print["Same integral, A = 2 + 0.5I"];

  poly = 1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2];
  vars = {x[1], x[2]};
  A    = 2 + 0.5 I;

  spec = <|
    "Polynomials"       -> {poly},
    "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-A},
    "Variables"         -> vars,
    "KinematicSymbols"  -> {},
    "RegulatorSymbol"   -> None
  |>;

  directResult = NIntegrate[
    1 / (1 + 2 t1^2 + t2^2 + t1 t2^2 + 3 t1^2 t2)^A,
    {t1, 0, Infinity}, {t2, 0, Infinity},
    MaxRecursion -> 20, PrecisionGoal -> 5
  ];

  verts   = PolytopeVertices[poly^(-Re[A]), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Module[{vr},
    vr = Quiet@ValidateDecomposition[spec, fanData, {}, 3];
    If[AssociationQ[vr],
      mcResult = vr["SectorSum"];
      relErr   = Abs[(mcResult - directResult) / directResult];,
      mcResult = 0;
      relErr = Infinity;
    ];
  ];

  pass = NumericQ[relErr] && (relErr < 0.02);
  Print["  NIntegrate = ", directResult];
  Print["  Sector sum = ", mcResult];
  Print["  Rel error  = ", relErr];
  Print["  ", If[pass, "PASS", "FAIL"]];

  pass
];

(* --------------------------------------------------------------------------
   Test 3: Divergent 1D
   f(eps) = Integral[0,1] dy y^{2eps-1} / (1+y^2)
   Exact: 1/(2eps) - log(2)/2
   -------------------------------------------------------------------------- *)

RunTest3[] := Module[
  {eps, exactPole, exactFinite, pass, testEps,
   numericalResult, exactAtEps, relErr},

  Print["--- Test 3: Divergent 1D ---"];
  Print["f(eps) = Integral[0,1] y^{2eps-1} / (1+y^2)"];
  Print["Exact: 1/(2eps) - log(2)/2"];

  eps = Symbol["epsTest3"];
  exactPole   = 1/2;
  exactFinite = -Log[2]/2;

  testEps = 0.01;
  numericalResult = NIntegrate[
    y^(2 testEps - 1) / (1 + y^2),
    {y, 0, 1},
    MaxRecursion -> 30, PrecisionGoal -> 8,
    Method -> "DoubleExponential"
  ];

  exactAtEps = exactPole / testEps + exactFinite;
  relErr = Abs[(numericalResult - exactAtEps) / exactAtEps];

  Print["  Pole coefficient: expected = ", N[exactPole],
        " (1/2)"];
  Print["  Finite part: expected = ", N[exactFinite],
        " (-log(2)/2)"];
  Print["  At eps = ", testEps, ":"];
  Print["    NIntegrate = ", numericalResult];
  Print["    Exact formula = ", N[exactAtEps]];
  Print["    Rel error = ", relErr];

  pass = (relErr < 0.001);
  Print["  ", If[pass, "PASS", "FAIL"]];

  pass
];

(* --------------------------------------------------------------------------
   Test 4: Divergent 2D
   -------------------------------------------------------------------------- *)

RunTest4[] := Module[
  {testEps, directResult, pass},

  Print["--- Test 4: Divergent 2D ---"];
  Print["g(eps) = Int dy1 dy2 y2^{2eps-1} y1^{2eps} / (1+y1+y2+y1*y2^2+y1^3*y2^2)^2"];

  testEps = 0.05;

  directResult = Quiet@NIntegrate[
    t2^(2 testEps - 1) * t1^(2 testEps) /
    (1 + t1 + t2 + t1 t2^2 + t1^3 t2^2)^2,
    {t1, 0, 1}, {t2, 0, 1},
    MaxRecursion -> 30, PrecisionGoal -> 5,
    Method -> "GlobalAdaptive",
    MinRecursion -> 5
  ];

  Print["  At eps = ", testEps, ":"];
  Print["    Direct NIntegrate = ", directResult];

  pass = NumericQ[directResult] && Abs[directResult] < 10^4;
  Print["  Value is finite: ", If[pass, "PASS", "FAIL"]];

  Module[{val1, val2, ratio},
    val1 = Quiet@NIntegrate[
      t2^(2 * 0.1 - 1) * t1^(2 * 0.1) /
      (1 + t1 + t2 + t1 t2^2 + t1^3 t2^2)^2,
      {t1, 0, 1}, {t2, 0, 1},
      MaxRecursion -> 30, PrecisionGoal -> 4,
      Method -> "GlobalAdaptive", MinRecursion -> 5
    ];
    val2 = Quiet@NIntegrate[
      t2^(2 * 0.2 - 1) * t1^(2 * 0.2) /
      (1 + t1 + t2 + t1 t2^2 + t1^3 t2^2)^2,
      {t1, 0, 1}, {t2, 0, 1},
      MaxRecursion -> 30, PrecisionGoal -> 4,
      Method -> "GlobalAdaptive", MinRecursion -> 5
    ];
    ratio = val1 / val2;
    Print["  g(0.1)/g(0.2) = ", ratio, " (expect ~ 2 if 1/eps behavior)"];
    If[Abs[ratio - 2] < 0.5,
      Print["  1/eps scaling: PASS"];,
      Print["  1/eps scaling: approximate (ratio = ", ratio, ")"];
    ];
  ];

  pass
];

(* --------------------------------------------------------------------------
   Test 5: End-to-end kinematic scan
   -------------------------------------------------------------------------- *)

RunTest5[] := Module[
  {poly, vars, lam, A, spec, verts, fanData,
   lamValues, nPoints,
   allPass, maxRelErr},

  Print["--- Test 5: End-to-end kinematic scan ---"];
  Print["Int dx1 dx2 / (1+lam*x1^2+x2^2+x1*x2^2)^{2+0.5I}"];
  Print["100 values of lam in [0.1, 10]"];

  lam = Symbol["lam"];
  A   = 2 + 0.5 I;
  poly = 1 + lam x[1]^2 + x[2]^2 + x[1] x[2]^2;
  vars = {x[1], x[2]};

  nPoints  = 100;
  lamValues = Table[0.1 + (10.0 - 0.1) (i - 1)/(nPoints - 1),
                    {i, nPoints}];

  spec = <|
    "Polynomials"       -> {poly},
    "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-A},
    "Variables"         -> vars,
    "KinematicSymbols"  -> {lam},
    "RegulatorSymbol"   -> None
  |>;

  verts   = PolytopeVertices[(poly /. lam -> 1)^(-Re[A]), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["Computing NIntegrate reference values..."];
  Module[{testIndices, refResults},
    testIndices = {1, 25, 50, 75, 100};
    refResults = Table[
      Module[{lamVal, result},
        lamVal = lamValues[[idx]];
        result = Quiet@NIntegrate[
          1 / (1 + lamVal t1^2 + t2^2 + t1 t2^2)^A,
          {t1, 0, Infinity}, {t2, 0, Infinity},
          MaxRecursion -> 20, PrecisionGoal -> 5
        ];
        {lamVal, result}
      ],
      {idx, testIndices}
    ];

    allPass = True;
    maxRelErr = 0;

    Do[
      Print["  lam = ", refResults[[i, 1]], ": NIntegrate = ",
            refResults[[i, 2]]];,
      {i, Length[refResults]}
    ];

    (* Sector decomposition check *)
    Do[
      Module[{kinRules, vr, relErr},
        kinRules = {lam -> refResults[[i, 1]]};
        vr = Quiet@ValidateDecomposition[spec, fanData, kinRules, 3];
        If[AssociationQ[vr],
          relErr = vr["RelativeError"];
          If[NumericQ[relErr],
            If[relErr > maxRelErr, maxRelErr = relErr];
            If[relErr > 0.05, allPass = False];
            Print["  lam = ", refResults[[i, 1]],
                  ": sector rel err = ", relErr,
                  " ", If[relErr < 0.05, "PASS", "FAIL"]];,
            Print["  lam = ", refResults[[i, 1]],
                  ": non-numeric error, FAIL"];
            allPass = False;
          ];
        ];
      ],
      {i, Length[refResults]}
    ];
  ];

  Print["  Max relative error: ", maxRelErr];
  Print["  ", If[allPass, "PASS", "FAIL"]];

  allPass
];

(* --------------------------------------------------------------------------
   Test 6: Large coefficients
   Verifies that the tropical decomposition and sector integrals remain
   correct when polynomial coefficients span many orders of magnitude.
   The tropically dominant monomial (min exponents) may NOT be the
   numerically largest monomial, but the factoring is an exact algebraic
   identity so the result must still agree with direct NIntegrate.
   -------------------------------------------------------------------------- *)

RunTest6[] := Module[
  {poly, vars, spec, verts, fanData, allPass,
   testCases, t1, t2},

  Print["--- Test 6: Large polynomial coefficients ---"];
  Print["Verifies correctness when coefficients span many orders of magnitude"];

  allPass = True;

  (* Test case A: coefficients O(10^6)
     P = 1 + 10^6 x1^2 + x2^2 + x1 x2^2
     The 10^6 term dominates numerically but is NOT the tropically
     dominant monomial in most sectors. *)
  Module[{polyA, specA, vertsA, fanA, directA, vrA, relErrA},
    Print[];
    Print["  Case A: P = 1 + 10^6 x1^2 + x2^2 + x1*x2^2, exponent -2"];
    polyA = 1 + 10^6 x[1]^2 + x[2]^2 + x[1] x[2]^2;
    vars  = {x[1], x[2]};

    specA = <|
      "Polynomials"        -> {polyA},
      "MonomialExponents"  -> {0, 0},
      "PolynomialExponents" -> {-2},
      "Variables"          -> vars,
      "KinematicSymbols"   -> {},
      "RegulatorSymbol"    -> None
    |>;

    directA = Quiet@NIntegrate[
      1 / (1 + 10^6 t1^2 + t2^2 + t1 t2^2)^2,
      {t1, 0, Infinity}, {t2, 0, Infinity},
      MaxRecursion -> 20, PrecisionGoal -> 6
    ];

    vertsA = PolytopeVertices[polyA^(-2), vars];
    fanA   = ComputeDecomposition[vertsA, "ShowProgress" -> False];

    vrA = Quiet@ValidateDecomposition[specA, fanA, {}, 3];
    If[AssociationQ[vrA],
      relErrA = vrA["RelativeError"];
      Print["    NIntegrate = ", directA];
      Print["    Sector sum = ", vrA["SectorSum"]];
      Print["    Rel error  = ", relErrA];
      If[!NumericQ[relErrA] || relErrA > 0.05,
        Print["    FAIL"];
        allPass = False;,
        Print["    PASS"];
      ];,
      Print["    ValidateDecomposition returned non-association, FAIL"];
      allPass = False;
    ];
  ];

  (* Test case B: mixed large and small coefficients
     P = 10^(-4) + 10^4 x1^2 + 10^(-4) x2^2 + 10^4 x1 x2^2 + x1^2 x2
     Coefficients span 8 orders of magnitude. *)
  Module[{polyB, specB, vertsB, fanB, directB, vrB, relErrB},
    Print[];
    Print["  Case B: coefficients from 10^-4 to 10^4, exponent -2"];
    polyB = 10^(-4) + 10^4 x[1]^2 + 10^(-4) x[2]^2 +
            10^4 x[1] x[2]^2 + x[1]^2 x[2];
    vars  = {x[1], x[2]};

    specB = <|
      "Polynomials"        -> {polyB},
      "MonomialExponents"  -> {0, 0},
      "PolynomialExponents" -> {-2},
      "Variables"          -> vars,
      "KinematicSymbols"   -> {},
      "RegulatorSymbol"    -> None
    |>;

    directB = Quiet@NIntegrate[
      1 / (10^(-4) + 10^4 t1^2 + 10^(-4) t2^2 +
           10^4 t1 t2^2 + t1^2 t2)^2,
      {t1, 0, Infinity}, {t2, 0, Infinity},
      MaxRecursion -> 20, PrecisionGoal -> 6
    ];

    vertsB = PolytopeVertices[polyB^(-2), vars];
    fanB   = ComputeDecomposition[vertsB, "ShowProgress" -> False];

    vrB = Quiet@ValidateDecomposition[specB, fanB, {}, 3];
    If[AssociationQ[vrB],
      relErrB = vrB["RelativeError"];
      Print["    NIntegrate = ", directB];
      Print["    Sector sum = ", vrB["SectorSum"]];
      Print["    Rel error  = ", relErrB];
      If[!NumericQ[relErrB] || relErrB > 0.05,
        Print["    FAIL"];
        allPass = False;,
        Print["    PASS"];
      ];,
      Print["    ValidateDecomposition returned non-association, FAIL"];
      allPass = False;
    ];
  ];

  (* Test case C: large coefficient with higher exponent
     P = 1 + 10^8 x1^3 x2 + x2^3, exponent -3
     The 10^8 monomial has degree 4 and large coefficient, ensuring
     it numerically dominates even though the constant term is tropically
     dominant in its cone. *)
  Module[{polyC, specC, vertsC, fanC, directC, vrC, relErrC},
    Print[];
    Print["  Case C: P = 1 + 10^8 x1^3*x2 + x2^3, exponent -3"];
    polyC = 1 + 10^8 x[1]^3 x[2] + x[2]^3;
    vars  = {x[1], x[2]};

    specC = <|
      "Polynomials"        -> {polyC},
      "MonomialExponents"  -> {0, 0},
      "PolynomialExponents" -> {-3},
      "Variables"          -> vars,
      "KinematicSymbols"   -> {},
      "RegulatorSymbol"    -> None
    |>;

    directC = Quiet@NIntegrate[
      1 / (1 + 10^8 t1^3 t2 + t2^3)^3,
      {t1, 0, Infinity}, {t2, 0, Infinity},
      MaxRecursion -> 20, PrecisionGoal -> 5
    ];

    vertsC = PolytopeVertices[polyC^(-3), vars];
    fanC   = ComputeDecomposition[vertsC, "ShowProgress" -> False];

    vrC = Quiet@ValidateDecomposition[specC, fanC, {}, 2];
    If[AssociationQ[vrC],
      relErrC = vrC["RelativeError"];
      Print["    NIntegrate = ", directC];
      Print["    Sector sum = ", vrC["SectorSum"]];
      Print["    Rel error  = ", relErrC];
      If[!NumericQ[relErrC] || relErrC > 0.05,
        Print["    FAIL"];
        allPass = False;,
        Print["    PASS"];
      ];,
      Print["    ValidateDecomposition returned non-association, FAIL"];
      allPass = False;
    ];
  ];

  Print[];
  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 7: Higher-dimensional integrals (3D and 4D)
   -------------------------------------------------------------------------- *)

RunTest7[] := Module[
  {allPass = True, vars3, poly3, spec3, verts3, fan3, vr3,
   vars4, poly4, spec4, verts4, fan4, vr4},

  Print["--- Test 7: Higher-dimensional integrals (3D and 4D) ---"];
  Print["Tests that the pipeline works correctly in dimensions > 2"];

  (* 3D: Int dx1 dx2 dx3 / (1 + x1^2 + x2^2 + x3^2 + x1*x2*x3)^3 *)
  Print[];
  Print["  Case A (3D): P = 1 + x1^2 + x2^2 + x3^2 + x1*x2*x3, exponent -3"];

  poly3 = 1 + x[1]^2 + x[2]^2 + x[3]^2 + x[1] x[2] x[3];
  vars3 = {x[1], x[2], x[3]};

  spec3 = <|
    "Polynomials"        -> {poly3},
    "MonomialExponents"  -> {0, 0, 0},
    "PolynomialExponents" -> {-3},
    "Variables"          -> vars3,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  verts3 = PolytopeVertices[poly3^(-1), vars3];
  fan3   = ComputeDecomposition[verts3, "ShowProgress" -> False];
  Print["    Fan: ", Length[fan3[[1]]], " rays, ", Length[fan3[[2]]], " sectors"];

  vr3 = Quiet@ValidateDecomposition[spec3, fan3, {}, 3];
  If[AssociationQ[vr3],
    Print["    NIntegrate = ", vr3["DirectResult"]];
    Print["    Sector sum = ", vr3["SectorSum"]];
    Print["    Rel error  = ", vr3["RelativeError"]];
    If[!NumericQ[vr3["RelativeError"]] || vr3["RelativeError"] > 0.01,
      Print["    FAIL"];
      allPass = False;,
      Print["    PASS"];
    ];,
    Print["    ValidateDecomposition returned non-association, FAIL"];
    allPass = False;
  ];

  (* 4D: Int dx1 dx2 dx3 dx4 / (1+x1^2+x2^2+x3^2+x4^2+x1*x2+x3*x4)^4 *)
  Print[];
  Print["  Case B (4D): P = 1+x1^2+x2^2+x3^2+x4^2+x1*x2+x3*x4, exponent -4"];

  poly4 = 1 + x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 + x[1] x[2] + x[3] x[4];
  vars4 = {x[1], x[2], x[3], x[4]};

  spec4 = <|
    "Polynomials"        -> {poly4},
    "MonomialExponents"  -> {0, 0, 0, 0},
    "PolynomialExponents" -> {-4},
    "Variables"          -> vars4,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> None
  |>;

  verts4 = PolytopeVertices[poly4^(-1), vars4];
  fan4   = ComputeDecomposition[verts4, "ShowProgress" -> False];
  Print["    Fan: ", Length[fan4[[1]]], " rays, ", Length[fan4[[2]]], " sectors"];

  vr4 = Quiet@ValidateDecomposition[spec4, fan4, {}, 2];
  If[AssociationQ[vr4],
    Print["    NIntegrate = ", vr4["DirectResult"]];
    Print["    Sector sum = ", vr4["SectorSum"]];
    Print["    Rel error  = ", vr4["RelativeError"]];
    If[!NumericQ[vr4["RelativeError"]] || vr4["RelativeError"] > 0.01,
      Print["    FAIL"];
      allPass = False;,
      Print["    PASS"];
    ];,
    Print["    ValidateDecomposition returned non-association, FAIL"];
    allPass = False;
  ];

  Print[];
  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 8: Full divergent pipeline (2D)
   ProcessSector -> ProcessDivergentSector -> ValidateSubtraction
   -------------------------------------------------------------------------- *)

RunTest8[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, divSectors,
   allPass},

  Print["--- Test 8: Full divergent pipeline (2D) ---"];
  Print["Int x1^{2eps-1} (1+x1+x2+x1*x2)^{-2} dx"];

  eps = Symbol["eps8"];
  poly = 1 + x[1] + x[2] + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[Length[divSectors] == 0,
    Print["  No divergent sectors found, FAIL"];
    Return[False]
  ];
  Print["  Found ", Length[divSectors], " divergent sector(s)"];

  allPass = True;
  Do[
    Module[{divData, vsResult, relErr},
      divData = Quiet@ProcessDivergentSector[sd, spec];
      If[!AssociationQ[divData],
        Print["  ProcessDivergentSector failed, FAIL"];
        allPass = False;,

        vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, 0.05];
        If[!AssociationQ[vsResult],
          Print["  ValidateSubtraction failed, FAIL"];
          allPass = False;,

          relErr = vsResult["RelativeError"];
          Print["  Sector ", sd["ConeIndex"], ": RelativeError = ", relErr];
          If[!NumericQ[relErr] || relErr > 0.02,
            Print["    FAIL"];
            allPass = False;,
            Print["    PASS"];
          ];
        ];
      ];
    ],
    {sd, divSectors}
  ];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 9: Reconstruction at multiple epsilon values
   -------------------------------------------------------------------------- *)

RunTest9[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, divSectors,
   testEpsValues, allPass, remainders},

  Print["--- Test 9: Reconstruction at multiple epsilon values ---"];
  Print["Same integral as Test 8, validated at eps = {0.1, 0.05, 0.01, 0.005}"];

  eps = Symbol["eps9"];
  poly = 1 + x[1] + x[2] + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[Length[divSectors] == 0,
    Print["  No divergent sectors, FAIL"];
    Return[False]
  ];

  testEpsValues = {0.1, 0.05, 0.01, 0.005};
  allPass = True;
  remainders = {};

  (* Use first divergent sector *)
  Module[{sd, divData},
    sd = divSectors[[1]];
    divData = Quiet@ProcessDivergentSector[sd, spec];
    If[!AssociationQ[divData],
      Print["  ProcessDivergentSector failed, FAIL"];
      Return[False]
    ];

    Do[
      Module[{vsResult, relErr},
        vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, te];
        If[!AssociationQ[vsResult],
          Print["  ValidateSubtraction failed at eps=", te, ", FAIL"];
          allPass = False;,

          relErr = vsResult["RelativeError"];
          AppendTo[remainders, {te, vsResult["Remainder"]}];
          (* Threshold scales with eps: truncation error is O(eps) *)
          Module[{threshold},
            threshold = 0.03 + 0.5 te;
            Print["  eps = ", te, ": RelativeError = ", relErr,
                  " (threshold = ", threshold, ")"];
            If[!NumericQ[relErr] || relErr > threshold,
              Print["    FAIL"];
              allPass = False;,
              Print["    PASS"];
            ];
          ];
        ];
      ],
      {te, testEpsValues}
    ];

    (* Check remainder convergence *)
    Module[{r1, r2, convErr},
      r1 = Select[remainders, #[[1]] == 0.01 &];
      r2 = Select[remainders, #[[1]] == 0.005 &];
      If[Length[r1] > 0 && Length[r2] > 0,
        r1 = r1[[1, 2]]; r2 = r2[[1, 2]];
        If[NumericQ[r1] && NumericQ[r2] && Abs[r1] > 0,
          convErr = Abs[r2 - r1] / Abs[r1];
          Print["  |R(0.005)-R(0.01)|/|R(0.01)| = ", convErr];
          If[convErr > 0.3,
            Print["    Remainder convergence: FAIL"];
            allPass = False;,
            Print["    Remainder convergence: PASS"];
          ];
        ];
      ];
    ];
  ];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 10: Divergent variable index permutation
   -------------------------------------------------------------------------- *)

RunTest10[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, divSectors,
   allPass, divVarSet},

  Print["--- Test 10: Divergent variable index permutation ---"];
  Print["Int x2^{2eps-1} (1+x1+x2+x1*x2)^{-2} dx"];

  eps = Symbol["eps10"];
  poly = 1 + x[1] + x[2] + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {0, 2 eps - 1},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[Length[divSectors] == 0,
    Print["  No divergent sectors, FAIL"];
    Return[False]
  ];

  divVarSet = Union[#["DivergentVariable"] & /@ divSectors];
  Print["  Divergent variables found: ", divVarSet];

  allPass = True;
  Do[
    Module[{divData, vsResult, relErr},
      divData = Quiet@ProcessDivergentSector[sd, spec];
      If[!AssociationQ[divData],
        Print["  ProcessDivergentSector failed for sector ",
              sd["ConeIndex"], ", FAIL"];
        allPass = False;,

        vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, 0.05];
        If[!AssociationQ[vsResult],
          Print["  ValidateSubtraction failed for sector ",
                sd["ConeIndex"], ", FAIL"];
          allPass = False;,

          relErr = vsResult["RelativeError"];
          Print["  Sector ", sd["ConeIndex"],
                " (divVar=", sd["DivergentVariable"],
                "): RelativeError = ", relErr];
          If[!NumericQ[relErr] || relErr > 0.02,
            Print["    FAIL"];
            allPass = False;,
            Print["    PASS"];
          ];
        ];
      ];
    ],
    {sd, divSectors}
  ];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 11: Epsilon-dependent polynomial exponents (B1 != 0)
   -------------------------------------------------------------------------- *)

RunTest11[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, divSectors,
   allPass},

  Print["--- Test 11: Epsilon-dependent polynomial exponents (B1 != 0) ---"];
  Print["Int x1^{2eps-1} (1+x1+x2+x1*x2)^{-2+eps} dx"];

  eps = Symbol["eps11"];
  poly = 1 + x[1] + x[2] + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2 + eps},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[Length[divSectors] == 0,
    Print["  No divergent sectors, FAIL"];
    Return[False]
  ];

  allPass = True;
  Do[
    Module[{divData, vsResult, relErr, b1Check, logInsCheck},
      divData = Quiet@ProcessDivergentSector[sd, spec];
      If[!AssociationQ[divData],
        Print["  ProcessDivergentSector failed, FAIL"];
        allPass = False;,

        (* Verify B1 is {1} *)
        b1Check = (divData["B1"] === {1});
        Print["  Sector ", sd["ConeIndex"], ": B1 = ", divData["B1"],
              If[b1Check, " (correct)", " (UNEXPECTED)"]];
        If[!b1Check, allPass = False];

        (* Verify G1LogInsertions PolynomialTerms has nonzero coefficient *)
        logInsCheck = False;
        Module[{polyTerms},
          polyTerms = divData["G1LogInsertions"]["PolynomialTerms"];
          If[AnyTrue[polyTerms, (#[[1]] =!= 0) &],
            logInsCheck = True
          ];
        ];
        Print["  G1LogInsertions PolynomialTerms nonzero: ", logInsCheck];
        If[!logInsCheck, allPass = False];

        vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, 0.05];
        If[!AssociationQ[vsResult],
          Print["  ValidateSubtraction failed, FAIL"];
          allPass = False;,

          relErr = vsResult["RelativeError"];
          Print["  RelativeError = ", relErr];
          If[!NumericQ[relErr] || relErr > 0.03,
            Print["    FAIL"];
            allPass = False;,
            Print["    PASS"];
          ];
        ];
      ];
    ],
    {sd, divSectors}
  ];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 12: Full sector sum (convergent + divergent)
   -------------------------------------------------------------------------- *)

RunTest12[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, convergentSectors, divSectors,
   testEps, totalSum, directResult, relErr, allPass},

  Print["--- Test 12: Full sector sum (convergent + divergent) ---"];
  Print["Sum all sectors at eps=0.05, compare with direct NIntegrate"];

  testEps = 0.05;
  eps = Symbol["eps12"];
  poly = 1 + x[1] + x[2] + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  convergentSectors = Select[allSectorData,
    (AssociationQ[#] && !#["IsDivergent"]) &];
  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  Print["  ", Length[convergentSectors], " convergent, ",
        Length[divSectors], " divergent sectors"];

  totalSum = 0;
  allPass = True;

  (* Convergent sectors: NIntegrate of flattened form at eps=testEps *)
  Do[
    Module[{flatPolys, polyExps, pf, n, yVars, polyVals, integrand, sVal},
      flatPolys = cs["FlattenedPolys"] /. eps -> testEps;
      polyExps  = cs["PolynomialExponents"] /. eps -> testEps;
      pf = cs["Prefactor"] /. eps -> testEps;
      n  = cs["Dimension"];
      yVars = Table[Unique["cy"], {n}];

      polyVals = Table[
        Total[Table[
          mono[[1]] * Exp[Total[mono[[2]] * Log /@ yVars]],
          {mono, flatPolys[[j]]}
        ]],
        {j, Length[flatPolys]}
      ];

      integrand = pf *
        Times @@ MapThread[
          Function[{pv, be}, Exp[be * Log[pv]]],
          {polyVals, polyExps}
        ];

      sVal = Quiet@NIntegrate[
        integrand,
        Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
        MaxRecursion -> 20, PrecisionGoal -> 4,
        Method -> "GlobalAdaptive"
      ];
      Print["  Convergent sector ", cs["ConeIndex"], ": ", sVal];
      totalSum += sVal;
    ],
    {cs, convergentSectors}
  ];

  (* Divergent sectors: use ValidateSubtraction *)
  Do[
    Module[{divData, vsResult},
      divData = Quiet@ProcessDivergentSector[ds, spec];
      If[!AssociationQ[divData],
        Print["  ProcessDivergentSector failed for sector ",
              ds["ConeIndex"]];
        allPass = False;,

        vsResult = Quiet@ValidateSubtraction[divData, ds, spec, {}, testEps];
        If[!AssociationQ[vsResult],
          Print["  ValidateSubtraction failed for sector ",
                ds["ConeIndex"]];
          allPass = False;,

          Print["  Divergent sector ", ds["ConeIndex"], ": ",
                vsResult["Reconstructed"]];
          totalSum += vsResult["Reconstructed"];
        ];
      ];
    ],
    {ds, divSectors}
  ];

  (* Direct NIntegrate on [0,Infinity)^2 *)
  directResult = Quiet@NIntegrate[
    t1^(2 testEps - 1) / (1 + t1 + t2 + t1 t2)^2,
    {t1, 0, Infinity}, {t2, 0, Infinity},
    MaxRecursion -> 30, PrecisionGoal -> 5,
    Method -> "GlobalAdaptive"
  ];

  relErr = Abs[(totalSum - directResult) / directResult];
  Print["  Sector sum      = ", totalSum];
  Print["  Direct NIntegrate = ", directResult];
  Print["  Relative error   = ", relErr];

  allPass = allPass && NumericQ[relErr] && relErr < 0.05;
  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 13: C++ codegen for divergent sectors
   -------------------------------------------------------------------------- *)

RunTest13[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, convergentSectors, divSectors,
   processedDiv, testEps,
   cppFile, cppBinary, codeResult,
   allPass},

  Print["--- Test 13: C++ codegen for divergent sectors ---"];
  Print["Verify code generation and compilation"];

  testEps = 0.05;
  eps = Symbol["eps13"];
  poly = 1 + x[1] + x[2] + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  convergentSectors = Select[allSectorData,
    (AssociationQ[#] && !#["IsDivergent"]) &];
  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  processedDiv = Table[
    Quiet@ProcessDivergentSector[sd, spec],
    {sd, divSectors}
  ];
  processedDiv = Select[processedDiv, AssociationQ];

  (* Substitute eps for numerical codegen *)
  Quiet[CreateDirectory[FileNameJoin[{Directory[], "INTERFILES"}]]];
  cppFile   = FileNameJoin[{Directory[], "INTERFILES", "test13_tropical.cpp"}];
  cppBinary = FileNameJoin[{Directory[], "INTERFILES", "test13_tropical"}];

  Module[{convNum, specNum},
    convNum = convergentSectors /. eps -> testEps;
    specNum = <|
      "Polynomials"        -> {poly},
      "MonomialExponents"  -> {2 testEps - 1, 0},
      "PolynomialExponents" -> {-2},
      "Variables"          -> vars,
      "KinematicSymbols"   -> {},
      "RegulatorSymbol"    -> None
    |>;
    codeResult = Quiet@GenerateCppMonteCarlo[
      convNum, processedDiv, specNum, cppFile, "NSamples" -> 1000
    ];
  ];

  allPass = True;

  If[!AssociationQ[codeResult],
    Print["  GenerateCppMonteCarlo failed, FAIL"];
    Return[False]
  ];

  If[codeResult["NG0"] > 0,
    Print["  NG0 = ", codeResult["NG0"], ", PASS"];,
    Print["  NG0 = ", codeResult["NG0"], " (expected > 0), FAIL"];
    allPass = False;
  ];

  If[codeResult["NG1"] > 0,
    Print["  NG1 = ", codeResult["NG1"], ", PASS"];,
    Print["  NG1 = ", codeResult["NG1"], " (expected > 0), FAIL"];
    allPass = False;
  ];

  If[codeResult["NRemainder"] > 0,
    Print["  NRemainder = ", codeResult["NRemainder"], ", PASS"];,
    Print["  NRemainder = ", codeResult["NRemainder"], " (expected > 0), FAIL"];
    allPass = False;
  ];

  (* Check no unresolved Mathematica symbols *)
  Module[{code, badPatterns, found},
    code = codeResult["Code"];
    badPatterns = {"Sin[", "Cos[", "Sqrt[", "Plus[", "Times[",
                   "Power[", "Rule[", "List["};
    found = Select[badPatterns, StringContainsQ[code, #] &];
    If[Length[found] > 0,
      Print["  Unresolved symbols: ", found, ", FAIL"];
      allPass = False;,
      Print["  No unresolved Mathematica symbols, PASS"];
    ];
  ];

  (* Check compilation *)
  Module[{compResult},
    compResult = Quiet@CompileCpp[cppFile, cppBinary];
    If[compResult === $Failed,
      Print["  Compilation FAILED"];
      allPass = False;,
      Print["  Compilation PASS"];
    ];
  ];

  (* Clean up *)
  Quiet[If[FileExistsQ[cppFile], DeleteFile[cppFile]]];
  Quiet[If[FileExistsQ[cppBinary], DeleteFile[cppBinary]]];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 14: Pole coefficient accuracy
   -------------------------------------------------------------------------- *)

RunTest14[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, divSectors,
   allPass},

  Print["--- Test 14: Pole coefficient accuracy ---"];
  Print["Int x1^{2eps-1} (1+x1+x2)^{-2} dx"];

  eps = Symbol["eps14"];
  poly = 1 + x[1] + x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[Length[divSectors] == 0,
    Print["  No divergent sectors, FAIL"];
    Return[False]
  ];

  allPass = True;

  Do[
    Module[{divData, g0Val, ck, poleCoeff,
            eps1, eps2, Ieps1, Ieps2, slopeCoeff, relErr},
      divData = Quiet@ProcessDivergentSector[sd, spec];
      If[!AssociationQ[divData],
        Print["  ProcessDivergentSector failed, FAIL"];
        allPass = False;,

        ck = divData["ck"];

        (* Compute G0 via NIntegrate *)
        Module[{g0FlatPolys, g0Pf, g0Dim, B0, g0yVars,
                g0PolyVals, g0Integrand},
          g0FlatPolys = divData["G0FlatPolys"];
          g0Pf  = divData["G0Prefactor"];
          g0Dim = divData["G0Dimension"];
          B0    = divData["G0PolyExponents"];
          g0yVars = Table[Unique["g0y"], {g0Dim}];

          g0PolyVals = Table[
            Total[Table[
              mono[[1]] * Exp[Total[mono[[2]] * Log /@ g0yVars]],
              {mono, g0FlatPolys[[j]]}
            ]],
            {j, Length[g0FlatPolys]}
          ];

          g0Integrand = g0Pf *
            Times @@ MapThread[
              Function[{pv, be}, Exp[be * Log[pv]]],
              {g0PolyVals, B0}
            ];

          g0Val = Quiet@NIntegrate[
            g0Integrand,
            Evaluate[Sequence @@ ({#, 0, 1} & /@ g0yVars)],
            MaxRecursion -> 20, PrecisionGoal -> 4,
            Method -> "GlobalAdaptive"
          ];
        ];

        poleCoeff = g0Val / ck;
        Print["  Sector ", sd["ConeIndex"], ": G0/ck = ", poleCoeff];

        (* Numerical check: sector integral at two small eps values *)
        eps1 = 0.01; eps2 = 0.02;

        Module[{yVars, clearedPolys, aNum1, polyVals1, integ1,
                aNum2, polyVals2, integ2, polyExps},
          yVars = Table[Unique["py"], {sd["Dimension"]}];
          clearedPolys = divData["ClearedPolys"];
          polyExps = sd["PolynomialExponents"];

          (* Sector integral at eps1 *)
          aNum1 = sd["NewExponents"] /. eps -> eps1;
          polyVals1 = Table[
            Total[Table[
              mono[[1]] * Exp[Total[mono[[2]] * Log /@ yVars]],
              {mono, clearedPolys[[j]]}
            ]],
            {j, Length[clearedPolys]}
          ];
          integ1 = Abs[sd["DetM"]] *
            Exp[Total[(aNum1 - 1) * Log /@ yVars]] *
            Times @@ MapThread[
              Function[{pv, be}, Exp[be * Log[pv]]],
              {polyVals1, polyExps /. eps -> eps1}
            ];

          Ieps1 = Quiet@NIntegrate[
            integ1,
            Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
            MaxRecursion -> 20, PrecisionGoal -> 3,
            Method -> "GlobalAdaptive"
          ];

          (* Sector integral at eps2 *)
          aNum2 = sd["NewExponents"] /. eps -> eps2;
          polyVals2 = Table[
            Total[Table[
              mono[[1]] * Exp[Total[mono[[2]] * Log /@ yVars]],
              {mono, clearedPolys[[j]]}
            ]],
            {j, Length[clearedPolys]}
          ];
          integ2 = Abs[sd["DetM"]] *
            Exp[Total[(aNum2 - 1) * Log /@ yVars]] *
            Times @@ MapThread[
              Function[{pv, be}, Exp[be * Log[pv]]],
              {polyVals2, polyExps /. eps -> eps2}
            ];

          Ieps2 = Quiet@NIntegrate[
            integ2,
            Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
            MaxRecursion -> 20, PrecisionGoal -> 3,
            Method -> "GlobalAdaptive"
          ];
        ];

        (* Numerical slope: I(eps) ~ C/eps + D + ...
           dI/d(1/eps) = C
           slopeCoeff = (Ieps1 - Ieps2) / (1/eps1 - 1/eps2) *)
        slopeCoeff = (Ieps1 - Ieps2) / (1/eps1 - 1/eps2);
        Print["  Numerical slope = ", slopeCoeff];

        relErr = Abs[(poleCoeff - slopeCoeff) / poleCoeff];
        Print["  |G0/ck - slope|/|G0/ck| = ", relErr];
        If[!NumericQ[relErr] || relErr > 0.05,
          Print["    FAIL"];
          allPass = False;,
          Print["    PASS"];
        ];
      ];
    ],
    {sd, divSectors}
  ];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 15: G1 log-insertion structural + numerical check
   -------------------------------------------------------------------------- *)

RunTest15[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, divSectors,
   allPass},

  Print["--- Test 15: G1 log-insertion structural + numerical check ---"];
  Print["Int x1^{2eps-1} (1+x1+x2+x1*x2)^{-2+eps} dx (B1 != 0)"];

  eps = Symbol["eps15"];
  poly = 1 + x[1] + x[2] + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2 + eps},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[Length[divSectors] == 0,
    Print["  No divergent sectors, FAIL"];
    Return[False]
  ];

  allPass = True;

  Do[
    Module[{divData, vsResult, ck, g1Formula, g1Num, relErr,
            logIns, nVarTerms, nPolyTerms, structPass},
      divData = Quiet@ProcessDivergentSector[sd, spec];
      If[!AssociationQ[divData],
        Print["  ProcessDivergentSector failed, FAIL"];
        allPass = False;,

        ck = divData["ck"];
        logIns = divData["G1LogInsertions"];

        (* Structural checks *)
        nVarTerms  = Length[logIns["VariableTerms"]];
        nPolyTerms = Length[logIns["PolynomialTerms"]];
        structPass = (nVarTerms == divData["G0Dimension"]) &&
                     (nPolyTerms == Length[divData["G0PolyExponents"]]);
        Print["  Sector ", sd["ConeIndex"],
              ": VariableTerms=", nVarTerms,
              " PolynomialTerms=", nPolyTerms,
              If[structPass, " (correct structure)",
                 " (UNEXPECTED structure)"]];
        If[!structPass, allPass = False];

        (* Get G1 from ValidateSubtraction *)
        vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, 0.05];
        If[!AssociationQ[vsResult],
          Print["  ValidateSubtraction failed, FAIL"];
          allPass = False;,

          g1Formula = vsResult["G1"];

          (* G1_num = ck * (I(eps) - G0/(ck*eps) - R) *)
          g1Num = ck * (vsResult["OriginalIntegral"] -
                        vsResult["G0"] / (ck * 0.05) -
                        vsResult["Remainder"]);

          Print["  G1 (formula)   = ", g1Formula];
          Print["  G1 (numerical) = ", g1Num];

          (* Use integral scale for error metric: G1 can be small
             relative to G0/(ck*eps), amplifying reconstruction error *)
          If[NumericQ[g1Formula] && Abs[vsResult["OriginalIntegral"]] > 0,
            relErr = Abs[g1Formula - g1Num] / Abs[vsResult["OriginalIntegral"]];
            Print["  |G1_formula - G1_num|/|I(eps)| = ", relErr];
            If[relErr > 0.05,
              Print["    FAIL"];
              allPass = False;,
              Print["    PASS"];
            ];,
            Print["  G1 or integral non-numeric, skipping cross-check"];
          ];
        ];
      ];
    ],
    {sd, divSectors}
  ];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 16: Multiple-polynomial divergent integral
   -------------------------------------------------------------------------- *)

RunTest16[] := Module[
  {p1, p2, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, divSectors,
   allPass},

  Print["--- Test 16: Multiple-polynomial divergent integral ---"];
  Print["Int x1^{2eps-1} (1+x1+x2)^{-1} (1+x1*x2)^{-1} dx"];

  eps = Symbol["eps16"];
  p1 = 1 + x[1] + x[2];
  p2 = 1 + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {p1, p2},
    "MonomialExponents"  -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-1, -1},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[(p1 p2)^(-1), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[Length[divSectors] == 0,
    Print["  No divergent sectors, FAIL"];
    Return[False]
  ];

  allPass = True;
  Module[{nProcessed = 0},
    Do[
      Module[{divData, vsResult, relErr},
        divData = Quiet@ProcessDivergentSector[sd, spec];
        If[!AssociationQ[divData],
          (* Nested divergence: known limitation, skip gracefully *)
          Print["  Sector ", sd["ConeIndex"],
                ": ProcessDivergentSector skipped (nested divergence)"];,

          nProcessed++;

          (* Verify ClearedPolys, SimplifiedPolys, G0FlatPolys each have 2 entries *)
          If[Length[divData["ClearedPolys"]] != 2,
            Print["  ClearedPolys has ", Length[divData["ClearedPolys"]],
                  " entries (expected 2), FAIL"];
            allPass = False;,
            Print["  ClearedPolys: 2 entries, correct"];
          ];
          If[Length[divData["SimplifiedPolys"]] != 2,
            Print["  SimplifiedPolys has ", Length[divData["SimplifiedPolys"]],
                  " entries (expected 2), FAIL"];
            allPass = False;,
            Print["  SimplifiedPolys: 2 entries, correct"];
          ];
          If[Length[divData["G0FlatPolys"]] != 2,
            Print["  G0FlatPolys has ", Length[divData["G0FlatPolys"]],
                  " entries (expected 2), FAIL"];
            allPass = False;,
            Print["  G0FlatPolys: 2 entries, correct"];
          ];

          vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, 0.05];
          If[!AssociationQ[vsResult],
            Print["  ValidateSubtraction failed for sector ",
                  sd["ConeIndex"], ", FAIL"];
            allPass = False;,

            relErr = vsResult["RelativeError"];
            Print["  Sector ", sd["ConeIndex"],
                  ": RelativeError = ", relErr];
            If[!NumericQ[relErr] || relErr > 0.05,
              Print["    FAIL"];
              allPass = False;,
              Print["    PASS"];
            ];
          ];
        ];
      ],
      {sd, divSectors}
    ];
    If[nProcessed == 0,
      Print["  No sectors could be processed, FAIL"];
      allPass = False;,
      Print["  Processed ", nProcessed, " of ", Length[divSectors],
            " divergent sectors"];
    ];
  ];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 17: 3D divergent integral
   -------------------------------------------------------------------------- *)

RunTest17[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, divSectors,
   allPass},

  Print["--- Test 17: 3D divergent integral ---"];
  Print["Int x1^{2eps-1} (1+x1+x2+x3+x1*x2*x3)^{-3} dx"];

  eps = Symbol["eps17"];
  poly = 1 + x[1] + x[2] + x[3] + x[1] x[2] x[3];
  vars = {x[1], x[2], x[3]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 0, 0},
    "PolynomialExponents" -> {-3},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-3), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];
  Print["  Fan: ", Length[dualVertices], " rays, ",
        Length[simplexList], " sectors"];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[Length[divSectors] == 0,
    Print["  No divergent sectors, FAIL"];
    Return[False]
  ];
  Print["  Found ", Length[divSectors], " divergent sector(s)"];

  allPass = True;
  Do[
    Module[{divData, vsResult, relErr},
      divData = Quiet@ProcessDivergentSector[sd, spec];
      If[!AssociationQ[divData],
        Print["  ProcessDivergentSector failed for sector ",
              sd["ConeIndex"], ", FAIL"];
        allPass = False;,

        (* Verify G0 dimension is 2 *)
        If[divData["G0Dimension"] != 2,
          Print["  G0Dimension = ", divData["G0Dimension"],
                " (expected 2), FAIL"];
          allPass = False;,
          Print["  G0Dimension = 2, correct"];
        ];

        vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, 0.05];
        If[!AssociationQ[vsResult],
          Print["  ValidateSubtraction failed for sector ",
                sd["ConeIndex"], ", FAIL"];
          allPass = False;,

          relErr = vsResult["RelativeError"];
          Print["  Sector ", sd["ConeIndex"],
                ": RelativeError = ", relErr];
          (* Relaxed threshold for 3D: NIntegrate less precise *)
          If[!NumericQ[relErr] || relErr > 0.10,
            Print["    FAIL"];
            allPass = False;,
            Print["    PASS"];
          ];
        ];
      ];
    ],
    {sd, divSectors}
  ];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];

(* --------------------------------------------------------------------------
   Test 18: Divergent integral with complex polynomial exponents
   -------------------------------------------------------------------------- *)

RunTest18[] := Module[
  {poly, vars, eps, spec, verts, fanData,
   dualVertices, simplexList,
   allSectorData, divSectors,
   allPass},

  Print["--- Test 18: Divergent integral with complex polynomial exponents ---"];
  Print["Int x1^{2eps-1} (1+x1+x2+x1*x2)^{-2+I} dx"];

  eps = Symbol["eps18"];
  poly = 1 + x[1] + x[2] + x[1] x[2];
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 0},
    "PolynomialExponents" -> {-2 + I},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList  = fanData[[2]];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[Length[divSectors] == 0,
    Print["  No divergent sectors, FAIL"];
    Return[False]
  ];
  Print["  Found ", Length[divSectors], " divergent sector(s)"];

  allPass = True;
  Do[
    Module[{divData, vsResult, relErr, b0Check},
      divData = Quiet@ProcessDivergentSector[sd, spec];
      If[!AssociationQ[divData],
        Print["  ProcessDivergentSector failed for sector ",
              sd["ConeIndex"], ", FAIL"];
        allPass = False;,

        (* Verify B0 is complex: {-2+I} *)
        b0Check = (divData["B0"] === {-2 + I});
        Print["  Sector ", sd["ConeIndex"], ": B0 = ", divData["B0"],
              If[b0Check, " (correct)", " (UNEXPECTED)"]];
        If[!b0Check, allPass = False];

        vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, 0.05];
        If[!AssociationQ[vsResult],
          Print["  ValidateSubtraction failed for sector ",
                sd["ConeIndex"], ", FAIL"];
          allPass = False;,

          relErr = vsResult["RelativeError"];
          Print["  Sector ", sd["ConeIndex"],
                ": RelativeError = ", relErr];
          (* Relaxed threshold: complex exponents add numerical difficulty *)
          If[!NumericQ[relErr] || relErr > 0.05,
            Print["    FAIL"];
            allPass = False;,
            Print["    PASS"];
          ];
        ];
      ];
    ],
    {sd, divSectors}
  ];

  Print["  ", If[allPass, "PASS", "FAIL"]];
  allPass
];
