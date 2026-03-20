(* ::Package:: *)

(* ============================================================================
   tropical_eval.wl

   Numerical evaluation of generalized Euler integrals via tropical
   decomposition.  Consumes the fan output of tropical_fan.wl and produces
   C++ Monte-Carlo code.

   Pipeline:
     Module 1 - ProcessSector        (symbolic coordinate transform + flattening)
     Module 2 - Divergence regulation (tropical subtraction for 1/eps poles)
     Module 3 - C++ code generation   (MmaToC + GenerateCppMonteCarlo)
     Module 4 - EvaluateTropicalMC    (driver: fan -> C++ -> results)
     Module 5 - RunAllTests           (validation suite)

   Dependencies: tropical_fan.wl (loaded automatically from same directory)
   ============================================================================ *)

(* Load tropical_fan.wl BEFORE BeginPackage so TropicalFan` is on $ContextPath.
   Set $SkipPolymakeLoad = True before loading this file to skip the polymake
   dependency (e.g. on cluster nodes without polymake installed). *)
If[!TrueQ[$SkipPolymakeLoad],
  Get[FileNameJoin[{DirectoryName[$InputFileName], "tropical_fan.wl"}]];
  BeginPackage["TropicalEval`", {"TropicalFan`"}],
  BeginPackage["TropicalEval`"]
];

(* ---- Public symbols ---- *)

ProcessSector::usage =
  "ProcessSector[integrandSpec, dualVertices, simplex, coneIndex] performs \
monomial change of variables and flattening for one simplicial cone.";

CheckFlatteningMagnitude::usage =
  "CheckFlatteningMagnitude[sectorData, nSamples] spot-checks the flattened \
integrand magnitude at random points.";

ValidateDecomposition::usage =
  "ValidateDecomposition[integrandSpec, fanData, testKinematics, precisionGoal] \
cross-checks sector sum against direct NIntegrate.";

IdentifyDivergences::usage =
  "IdentifyDivergences[sectorData, eps] identifies divergent variables.";

ProcessDivergentSector::usage =
  "ProcessDivergentSector[sectorData, integrandSpec] constructs G0, G1, \
remainder, and analytic pole for a divergent sector.";

ValidateSubtraction::usage =
  "ValidateSubtraction[sectorData, testKinematics, testEpsilon] checks \
self-consistency of the tropical subtraction.";

MmaToC::usage =
  "MmaToC[expr, paramMap] converts a Mathematica expression to a C++ string.";

GenerateCppMonteCarlo::usage =
  "GenerateCppMonteCarlo[convergentSectors, divergentSectors, integrandSpec, \
outputFile] generates self-contained C++ Monte Carlo source.";

EvaluateTropicalMC::usage =
  "EvaluateTropicalMC[integrandSpec, fanData, kinematicPoints, opts] runs \
the full tropical Monte Carlo pipeline.";

RunAllTests::usage =
  "RunAllTests[] runs the validation suite (17 tests) with structured reporting. \
Defined in tropical_eval_examples.wl.";

ParsePolynomial::usage =
  "ParsePolynomial[poly, vars] parses a polynomial into a list of \
{coefficient, exponentVector} pairs.";

CompileCpp::usage =
  "CompileCpp[srcFile, outputBinary, debug] compiles generated C++ \
Monte Carlo source code. debug=True adds -DTROPICAL_MC_DEBUG.";

(* ---- Error messages ---- *)

TropicalEval::degenerate = "Sector `1`: degenerate cone, det(M) = 0.";
TropicalEval::notsimplicial = "Sector `1`: `2` rays for `3` variables (not simplicial).";
TropicalEval::divergent = "Sector `1`: variable y_`2` divergent, a_`2` = `3`.";
TropicalEval::nested = "Sector `1`: multiple divergent variables (`2`). Nested subtraction not implemented.";
TropicalEval::badck = "Sector `1`: c_k = 0 for divergent variable y_`2`. Higher-order pole.";
TropicalEval::validate = "Validation `1`: relative error `2` exceeds tolerance `3`.";

(* ============================================================================
   PRIVATE IMPLEMENTATION
   ============================================================================ *)

Begin["`Private`"]

(* tropical_fan.wl is loaded before BeginPackage above *)

(* --------------------------------------------------------------------------
   Helper: ParsePolynomial
   Converts a polynomial into {coefficient, exponentVector} pairs.
   -------------------------------------------------------------------------- *)

ParsePolynomial[poly_, vars_List] := Module[
  {expanded, terms, result},
  expanded = Expand[poly];
  terms = If[Head[expanded] === Plus, List @@ expanded, {expanded}];
  result = Table[
    Module[{coeff, exps},
      exps = Exponent[term, #] & /@ vars;
      coeff = term / (Times @@ MapThread[Power, {vars, exps}]);
      {Simplify[coeff], exps}
    ],
    {term, terms}
  ];
  result
];

(* --------------------------------------------------------------------------
   Helper: TransformExponents
   Given original exponent vector and matrix M, compute new exponent vector.
   -------------------------------------------------------------------------- *)

TransformExponents[expVec_List, mMatrix_List] := expVec . mMatrix;

(* --------------------------------------------------------------------------
   MODULE 1: ProcessSector

   Key insight (tropical factoring):
   After the monomial substitution x_i = prod y_j^{M_ij}, each polynomial
   P_k becomes a sum of monomials in y with SIGNED exponents.  The dominant
   monomial (the one the tropical fan says dominates in this cone) has the
   minimum exponents.  We factor it out:
       P_k(y) = prod_j y_j^{d_{k,j}} * Q_k(y)
   where d_{k,j} = min_m (transformed exponent of y_j in monomial m),
   and Q_k has all non-negative y-exponents with a constant term.

   The effective monomial prefactor then becomes:
       a_j^eff = rawA_j + sum_k B_k * d_{k,j}

   For a properly constructed tropical fan and convergent integral,
   a_j^eff > 0 for all j, and we can flatten using these effective exponents.
   -------------------------------------------------------------------------- *)

Options[ProcessSector] = {"Verbose" -> False};

ProcessSector[integrandSpec_Association, dualVertices_List,
              simplex_List, coneIndex_Integer, OptionsPattern[]] :=
Module[
  {polys, monoExps, polyExps, vars, eps, kinSyms,
   selectedRays, mMatrix, detM, n,
   transformedPolys, clearedPolys, minExponents,
   rawAVals, effectiveAVals,
   flattenedPolys, prefactor,
   isDivergent, divVar, verbose, sectorData,
   parsedPolys},

  verbose = OptionValue["Verbose"];

  (* Extract fields from integrand spec *)
  polys    = integrandSpec["Polynomials"];
  monoExps = integrandSpec["MonomialExponents"];
  polyExps = integrandSpec["PolynomialExponents"];
  vars     = integrandSpec["Variables"];
  kinSyms  = integrandSpec["KinematicSymbols"];
  eps      = integrandSpec["RegulatorSymbol"];
  n        = Length[vars];

  (* --- Step 1: Monomial change of variables --- *)

  (* Extract ray vectors for this simplex (rows of dualVertices) *)
  selectedRays = dualVertices[[#]] & /@ simplex;

  (* Check: simplicial condition *)
  If[Length[selectedRays] != n,
    Message[TropicalEval::notsimplicial, coneIndex,
            Length[selectedRays], n];
    Return[$Failed]
  ];

  (* M_{ij} = -rho_j[i], i.e. M = -Transpose[selectedRays] *)
  mMatrix = -Transpose[selectedRays];

  (* Check: non-degenerate *)
  detM = Det[mMatrix];
  If[detM === 0 || TrueQ[detM == 0],
    Message[TropicalEval::degenerate, coneIndex];
    Return[$Failed]
  ];

  (* Raw exponents from monomial part + Jacobian:
     rawA_i = sum_k (A_k + 1) * M_{ki} = (monoExps + 1) . M *)
  rawAVals = (monoExps + 1) . mMatrix;

  (* --- Transform polynomials --- *)
  parsedPolys = ParsePolynomial[#, vars] & /@ polys;

  transformedPolys = Table[
    Table[
      Module[{coeff, origExp, newExp},
        coeff   = mono[[1]];
        origExp = mono[[2]];
        newExp  = TransformExponents[origExp, mMatrix];
        {coeff, newExp}
      ],
      {mono, parsedPolys[[j]]}
    ],
    {j, Length[polys]}
  ];

  (* --- Step 1b: Tropical factoring --- *)
  (* For each polynomial, find min exponents and factor them out *)

  minExponents = Table[
    Table[
      Min[#[[2, i]] & /@ transformedPolys[[j]]],
      {i, n}
    ],
    {j, Length[polys]}
  ];

  (* Cleared polynomials: shift exponents so minimum is 0 *)
  clearedPolys = Table[
    Table[
      {mono[[1]], mono[[2]] - minExponents[[j]]},
      {mono, transformedPolys[[j]]}
    ],
    {j, Length[polys]}
  ];

  (* Effective exponents: a_i^eff = rawA_i + sum_j B_j * d_{j,i} *)
  effectiveAVals = rawAVals + Total[
    Table[polyExps[[j]] * minExponents[[j]], {j, Length[polys]}]
  ];

  If[verbose,
    Print["Sector ", coneIndex, ": det(M) = ", detM,
          ", rawA = ", rawAVals,
          ", minExp = ", minExponents,
          ", effA = ", effectiveAVals]
  ];

  (* --- Step 2: Flattening using effective exponents --- *)
  (* Check convergence: Re(a_i^eff) > 0 for all i at eps = 0 *)
  isDivergent = False;
  divVar = 0;

  Module[{a0vals},
    a0vals = effectiveAVals /. (If[eps =!= None, eps -> 0, {}]);
    Do[
      If[TrueQ[Re[a0vals[[i]]] <= 0] ||
         (NumericQ[a0vals[[i]]] && Re[a0vals[[i]]] <= 0),
        isDivergent = True;
        divVar = i;
      ],
      {i, n}
    ];
  ];

  If[isDivergent,
    (* Divergent: return pre-flattened data for Module 2 *)
    sectorData = <|
      "ConeIndex"           -> coneIndex,
      "RayMatrix"           -> mMatrix,
      "DetM"                -> detM,
      "SelectedRays"        -> selectedRays,
      "RawExponents"        -> rawAVals,
      "NewExponents"        -> effectiveAVals,
      "MinExponents"        -> minExponents,
      "TransformedPolys"    -> transformedPolys,
      "ClearedPolys"        -> clearedPolys,
      "Prefactor"           -> Abs[detM],
      "IsDivergent"         -> True,
      "DivergentVariable"   -> divVar,
      "Dimension"           -> n,
      "PolynomialExponents" -> polyExps,
      "MonomialExponents"   -> monoExps
    |>;
    If[verbose,
      Print["  -> Divergent in variable y_", divVar]
    ];
    Return[sectorData]
  ];

  (* Convergent: flatten y_i -> (y_i')^{1/a_i^eff} *)
  (* Each monomial {coeff, {e1,...,en}} of Q_j becomes
     {coeff, {e1/a1^eff,...,en/an^eff}} *)
  flattenedPolys = Table[
    Table[
      {mono[[1]],
       MapThread[#1/#2 &, {mono[[2]], effectiveAVals}]},
      {mono, clearedPolys[[j]]}
    ],
    {j, Length[clearedPolys]}
  ];

  (* Prefactor = |det(M)| / Prod(a_i^eff) *)
  prefactor = Abs[detM] / (Times @@ effectiveAVals);

  sectorData = <|
    "ConeIndex"            -> coneIndex,
    "RayMatrix"            -> mMatrix,
    "DetM"                 -> detM,
    "SelectedRays"         -> selectedRays,
    "RawExponents"         -> rawAVals,
    "NewExponents"         -> effectiveAVals,
    "MinExponents"         -> minExponents,
    "TransformedPolys"     -> transformedPolys,
    "ClearedPolys"         -> clearedPolys,
    "FlattenedPolys"       -> flattenedPolys,
    "Prefactor"            -> prefactor,
    "IsDivergent"          -> False,
    "DivergentVariable"    -> 0,
    "Dimension"            -> n,
    "PolynomialExponents"  -> polyExps,
    "MonomialExponents"    -> monoExps
  |>;

  If[verbose,
    Print["  -> Convergent, prefactor = ", prefactor]
  ];

  sectorData
];

(* --------------------------------------------------------------------------
   CheckFlatteningMagnitude
   Spot-check that the flattened integrand is O(1) at random points.
   -------------------------------------------------------------------------- *)

CheckFlatteningMagnitude[sectorData_Association, nSamples_Integer: 20,
                         testKinematics_List: {}] :=
Module[
  {flatPolys, polyExps, prefactor, dim, mags, y, polyVals, integrandVal,
   kinRules},

  If[sectorData["IsDivergent"],
    Print["CheckFlatteningMagnitude: sector ", sectorData["ConeIndex"],
          " is divergent, skipping."];
    Return[Null]
  ];

  flatPolys = sectorData["FlattenedPolys"];
  polyExps  = sectorData["PolynomialExponents"];
  prefactor = sectorData["Prefactor"];
  dim       = sectorData["Dimension"];
  kinRules  = If[testKinematics === {}, {}, testKinematics];

  mags = Table[
    y = RandomReal[{0.01, 0.99}, dim];
    polyVals = Table[
      Total[
        Table[
          Module[{coeff, alphas, logY},
            coeff  = mono[[1]] /. kinRules;
            alphas = mono[[2]] /. kinRules;
            logY   = Log[y];
            coeff * Exp[Total[alphas * logY]]
          ],
          {mono, flatPolys[[j]]}
        ]
      ],
      {j, Length[flatPolys]}
    ];
    integrandVal = (prefactor /. kinRules) *
      Times @@ MapThread[
        Exp[#2 * Log[#1]] &,
        {polyVals, polyExps /. kinRules}
      ];
    Abs[integrandVal],
    {nSamples}
  ];

  Module[{meanMag, maxMag, minMag},
    meanMag = Mean[mags];
    maxMag  = Max[mags];
    minMag  = Min[mags];
    If[maxMag > 10^3 || minMag < 10^(-6),
      Print["WARNING: Sector ", sectorData["ConeIndex"],
            " flattening check: min=", minMag, " max=", maxMag,
            " mean=", meanMag]
    ];
    <|"Mean" -> meanMag, "Max" -> maxMag, "Min" -> minMag,
      "Samples" -> mags|>
  ]
];

(* --------------------------------------------------------------------------
   ValidateDecomposition
   Cross-check sector sum against direct NIntegrate.
   -------------------------------------------------------------------------- *)

ValidateDecomposition[integrandSpec_Association, fanData_List,
                      testKinematics_List, precisionGoal_Integer: 3] :=
Module[
  {polys, monoExps, polyExps, vars, n, dualVertices, simplexList,
   directIntegrand, directResult, sectorResults, sectorSum,
   relError, kinRules, allSectorData},

  polys    = integrandSpec["Polynomials"];
  monoExps = integrandSpec["MonomialExponents"];
  polyExps = integrandSpec["PolynomialExponents"];
  vars     = integrandSpec["Variables"];
  n        = Length[vars];
  kinRules = testKinematics;

  {dualVertices, simplexList} = fanData;

  (* Direct NIntegrate of the original integrand *)
  directIntegrand = (Times @@ MapThread[Power, {vars, monoExps}]) *
    (Times @@ MapThread[Power, {polys, polyExps}]) /. kinRules;

  directResult = NIntegrate[
    directIntegrand,
    Evaluate[Sequence @@ ({#, 0, Infinity} & /@ vars)],
    MaxRecursion -> 20,
    PrecisionGoal -> precisionGoal + 1,
    Method -> "GlobalAdaptive"
  ];

  (* Process all sectors *)
  allSectorData = Table[
    ProcessSector[integrandSpec /. kinRules, dualVertices,
                  simplexList[[s]], s],
    {s, Length[simplexList]}
  ];

  (* For each sector, evaluate via NIntegrate on [0,1]^n *)
  sectorResults = Table[
    Module[{sd, flatPolys, pExps, pf, dim, yVars, integrand,
            polyVals, result},
      sd = allSectorData[[s]];
      If[sd === $Failed, 0,
        If[sd["IsDivergent"],
          Print["ValidateDecomposition: sector ", s,
                " is divergent, skipping"];
          0,
          flatPolys = sd["FlattenedPolys"];
          pExps     = sd["PolynomialExponents"] /. kinRules;
          pf        = sd["Prefactor"] /. kinRules;
          dim       = sd["Dimension"];
          yVars     = Table[Unique["yv"], {dim}];

          (* Evaluate using log-exp form to avoid numerical issues *)
          polyVals = Table[
            Total[
              Table[
                Module[{coeff, alphas, logY},
                  coeff  = mono[[1]] /. kinRules;
                  alphas = mono[[2]] /. kinRules;
                  logY   = Log /@ yVars;
                  coeff * Exp[Total[alphas * logY]]
                ],
                {mono, flatPolys[[j]]}
              ]
            ],
            {j, Length[flatPolys]}
          ];

          integrand = pf *
            Times @@ MapThread[
              Function[{pv, be}, Exp[be * Log[pv]]],
              {polyVals, pExps}
            ];

          result = Quiet@NIntegrate[
            integrand,
            Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
            MaxRecursion -> 15,
            PrecisionGoal -> precisionGoal,
            Method -> "GlobalAdaptive"
          ];
          result
        ]
      ]
    ],
    {s, Length[simplexList]}
  ];

  sectorSum = Total[sectorResults];
  relError  = Abs[(sectorSum - directResult) / directResult];

  If[relError > 10^(-precisionGoal + 1),
    Message[TropicalEval::validate, "Decomposition", relError,
            10^(-precisionGoal + 1)]
  ];

  <|"DirectResult" -> directResult, "SectorSum" -> sectorSum,
    "RelativeError" -> relError, "SectorResults" -> sectorResults|>
];

(* --------------------------------------------------------------------------
   2D benchmark unit test (hard-coded)
   P = 1 + 2 x1^2 + x2^2 + x1 x2^2 + 3 x1^2 x2
   rho_1 = (0,1), rho_2 = (1,1), A_1 = A_2 = 0, polynomial exponent -A
   Expected effective: a_1 = 2A - 1, a_2 = 3A - 2
   -------------------------------------------------------------------------- *)

RunBenchmark2D[] := Module[
  {dualVerts, simplex, integrandSpec, sd, aExpected, aGot, A, pass},

  dualVerts = {{0, 1}, {1, 1}};
  simplex   = {1, 2};

  A = Symbol["Abench"];

  integrandSpec = <|
    "Polynomials"       -> {1 + 2 x[1]^2 + x[2]^2 + x[1] x[2]^2 + 3 x[1]^2 x[2]},
    "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-A},
    "Variables"         -> {x[1], x[2]},
    "KinematicSymbols"  -> {},
    "RegulatorSymbol"   -> None
  |>;

  sd = ProcessSector[integrandSpec, dualVerts, simplex, 1];

  (* Expected effective exponents:
     Raw a = (monoExps+1).M = {1,1}.{{0,-1},{-1,-1}} = {-1,-2}
     Transformed polynomial monomials in y have exponents:
       1         -> {0,0}.M = {0,0}
       2 x1^2    -> {2,0}.M = {0,-2}
       x2^2      -> {0,2}.M = {-2,-2}
       x1 x2^2   -> {1,2}.M = {-2,-3}     <- dominant (minimum)
       3 x1^2 x2 -> {2,1}.M = {-1,-3}
     minExp = {-2, -3}
     effA = rawA + B*minExp = {-1,-2} + (-A)*{-2,-3} = {2A-1, 3A-2} *)

  aExpected = {2 A - 1, 3 A - 2};
  aGot      = sd["NewExponents"];

  pass = TrueQ[Simplify[aGot - aExpected] === {0, 0}] ||
         (aGot === aExpected);
  If[!pass,
    Print["BENCHMARK 2D FAILED: expected a_eff = ", aExpected,
          " got ", aGot];
    ,
    Print["Benchmark 2D: PASSED (a_eff = ", aGot, ")"];
  ];
  pass
];


(* ============================================================================
   MODULE 2: DIVERGENCE REGULATION (Tropical Subtraction)
   ============================================================================ *)

(* --------------------------------------------------------------------------
   IdentifyDivergences
   Uses the EFFECTIVE exponents (post tropical factoring).
   -------------------------------------------------------------------------- *)

IdentifyDivergences[sectorData_Association, eps_] :=
Module[
  {aVals, n, a0, a1, divVars, ck, ak0, k},

  aVals = sectorData["NewExponents"];  (* effective exponents *)
  n     = sectorData["Dimension"];

  (* Expand a_i(epsilon) = a_i^(0) + a_i^(1) * epsilon + ... *)
  a0 = aVals /. eps -> 0;
  a1 = D[aVals, eps] /. eps -> 0;

  (* Find divergent variables: Re(a_k^(0)) <= 0 *)
  divVars = {};
  Do[
    If[TrueQ[Re[a0[[i]]] <= 0] ||
       (NumericQ[a0[[i]]] && Re[a0[[i]]] <= 0),
      AppendTo[divVars, i]
    ],
    {i, n}
  ];

  (* Check: at most one divergent variable *)
  If[Length[divVars] > 1,
    Message[TropicalEval::nested, sectorData["ConeIndex"], divVars];
    Return[$Failed]
  ];

  If[Length[divVars] == 0,
    Return[<|"IsDivergent" -> False|>]
  ];

  k   = divVars[[1]];
  ck  = a1[[k]];
  ak0 = a0[[k]];

  (* Check: c_k != 0 *)
  If[TrueQ[ck == 0] || (NumericQ[ck] && ck == 0),
    Message[TropicalEval::badck, sectorData["ConeIndex"], k];
    Return[$Failed]
  ];

  (* Check: all non-divergent variables have Re(a_i^(0)) > 0 *)
  Do[
    If[i != k && (TrueQ[Re[a0[[i]]] <= 0] ||
                  (NumericQ[a0[[i]]] && Re[a0[[i]]] <= 0)),
      Print["WARNING: Sector ", sectorData["ConeIndex"],
            ": non-divergent variable y_", i,
            " has Re(a_", i, "^(0)) = ", Re[a0[[i]]], " <= 0"]
    ],
    {i, n}
  ];

  Print["Sector ", sectorData["ConeIndex"],
        ": variable y_", k, " divergent, a_", k,
        "(eps) = ", ck, " * eps + ..., c_k = ", ck];

  <|"IsDivergent" -> True,
    "DivergentVariable" -> k,
    "ck" -> ck,
    "ak0" -> ak0,
    "a0" -> a0,
    "a1" -> a1|>
];

(* --------------------------------------------------------------------------
   ProcessDivergentSector
   Constructs G0, G1, remainder, and analytic pole.
   Works with the CLEARED polynomials Q_j (non-negative y-exponents)
   and the EFFECTIVE exponents.
   -------------------------------------------------------------------------- *)

ProcessDivergentSector[sectorData_Association, integrandSpec_Association] :=
Module[
  {eps, k, aVals, a0, a1, ck, n, polyExps,
   clearedPolys, detM,
   B0, B1, simpPolys, fullPolys,
   g0FlatPolys, g0Prefactor, g0Avals,
   g1LogInsertions,
   remainderData,
   analyticPole,
   divInfo},

  eps = integrandSpec["RegulatorSymbol"];
  n   = sectorData["Dimension"];

  (* Identify the divergent variable using effective exponents *)
  divInfo = IdentifyDivergences[sectorData, eps];
  If[divInfo === $Failed, Return[$Failed]];
  If[!divInfo["IsDivergent"],
    Print["ProcessDivergentSector called on non-divergent sector"];
    Return[$Failed]
  ];

  k      = divInfo["DivergentVariable"];
  ck     = divInfo["ck"];
  a0     = divInfo["a0"];
  a1     = divInfo["a1"];
  aVals  = sectorData["NewExponents"];  (* effective *)
  detM   = sectorData["DetM"];

  polyExps    = sectorData["PolynomialExponents"];
  clearedPolys = sectorData["ClearedPolys"];

  (* Polynomial exponents at eps=0 and first derivative *)
  B0 = polyExps /. eps -> 0;
  B1 = D[polyExps, eps] /. eps -> 0;

  (* --- Step 2: Construct F_simple by setting y_k = 0 --- *)
  (* In the cleared polynomials Q_j, keep only monomials where
     the exponent of y_k is 0 (these survive when y_k = 0). *)
  simpPolys = Table[
    Select[clearedPolys[[j]],
      (#[[2, k]] === 0 || TrueQ[#[[2, k]] == 0]) &
    ],
    {j, Length[clearedPolys]}
  ];

  fullPolys = clearedPolys;

  (* --- Step 3: G0 and G1 --- *)
  Module[{ndVars, g0aVals, g0Dim},
    ndVars  = DeleteCases[Range[n], k];
    g0aVals = a0[[ndVars]];
    g0Dim   = n - 1;

    (* For G0: flatten non-divergent variables *)
    g0FlatPolys = Table[
      Table[
        Module[{coeff, exps, ndExps, flatExps},
          coeff  = mono[[1]];
          exps   = mono[[2]];
          ndExps = exps[[ndVars]];
          flatExps = MapThread[#1/#2 &, {ndExps, g0aVals}];
          {coeff, flatExps}
        ],
        {mono, simpPolys[[j]]}
      ],
      {j, Length[simpPolys]}
    ];

    g0Prefactor = Abs[detM] / (Times @@ g0aVals);
    g0Avals = g0aVals;

    (* G1 log insertion factors *)
    g1LogInsertions = <|
      "VariableTerms" -> Table[
        {a1[[ndVars[[i]]]] / g0aVals[[i]], i},
        {i, g0Dim}
      ],
      "PolynomialTerms" -> Table[
        {B1[[j]], j},
        {j, Length[polyExps]}
      ]
    |>;
  ];

  (* --- Step 4: Finite remainder --- *)
  Module[{ndVars, remAvals, remFlatFullPolys, remFlatSimpPolys,
          remPrefactor},
    ndVars   = DeleteCases[Range[n], k];
    remAvals = a0[[ndVars]];

    (* Partially flatten: only non-div vars get flattened *)
    remFlatFullPolys = Table[
      Table[
        Module[{coeff, exps, flatExps},
          coeff = mono[[1]];
          exps  = mono[[2]];
          flatExps = Table[
            If[i == k,
              exps[[i]],
              exps[[i]] / a0[[i]]
            ],
            {i, n}
          ];
          {coeff, flatExps}
        ],
        {mono, fullPolys[[j]]}
      ],
      {j, Length[fullPolys]}
    ];

    remFlatSimpPolys = Table[
      Table[
        Module[{coeff, exps, flatExps},
          coeff = mono[[1]];
          exps  = mono[[2]];
          flatExps = Table[
            If[i == k,
              exps[[i]],
              exps[[i]] / a0[[i]]
            ],
            {i, n}
          ];
          {coeff, flatExps}
        ],
        {mono, simpPolys[[j]]}
      ],
      {j, Length[simpPolys]}
    ];

    remPrefactor = Abs[detM] / (Times @@ remAvals);

    remainderData = <|
      "FullPolys"     -> remFlatFullPolys,
      "SimplifiedPolys" -> remFlatSimpPolys,
      "Prefactor"     -> remPrefactor,
      "DivVarIndex"   -> k,
      "DivVarExp"     -> a0[[k]],
      "Dimension"     -> n,
      "NonDivVars"    -> DeleteCases[Range[n], k],
      "PolynomialExponents" -> B0
    |>;
  ];

  (* --- Step 5: Analytic contributions --- *)
  analyticPole = 1 / ck;

  <|
    "ConeIndex"           -> sectorData["ConeIndex"],
    "IsDivergent"         -> True,
    "DivergentVariable"   -> k,
    "ck"                  -> ck,
    "a0"                  -> a0,
    "a1"                  -> a1,
    "B0"                  -> B0,
    "B1"                  -> B1,
    "G0FlatPolys"         -> g0FlatPolys,
    "G0Prefactor"         -> g0Prefactor,
    "G0Avals"             -> g0Avals,
    "G0Dimension"         -> n - 1,
    "G0PolyExponents"     -> B0,
    "G1LogInsertions"     -> g1LogInsertions,
    "Remainder"           -> remainderData,
    "AnalyticPole"        -> analyticPole,
    "DetM"                -> detM,
    "Dimension"           -> n,
    "ClearedPolys"        -> clearedPolys,
    "SimplifiedPolys"     -> simpPolys,
    "TransformedPolys"    -> sectorData["TransformedPolys"],
    "MinExponents"        -> sectorData["MinExponents"],
    "PolynomialExponents" -> polyExps,
    "MonomialExponents"   -> sectorData["MonomialExponents"],
    "RawExponents"        -> sectorData["RawExponents"]
  |>
];

(* --------------------------------------------------------------------------
   ValidateSubtraction
   Self-consistency check at finite epsilon.
   Uses the ORIGINAL transformed polynomials (not cleared) to compute
   the true sector integral, then compares against the subtracted form.
   -------------------------------------------------------------------------- *)

ValidateSubtraction[divSectorData_Association, sectorData_Association,
                    integrandSpec_Association,
                    testKinematics_List, testEpsilon_: 0.01] :=
Module[
  {eps, n, k, ck, a0, a1, aVals, polyExps,
   kinRules, epsRules, fullRules,
   clearedPolys, simpPolys, minExps,
   detM, yVars,
   originalIntegral, g0Val, g1Val, remVal,
   divContrib, reconstructed, relError,
   rawAVals, B0, B1},

  eps      = integrandSpec["RegulatorSymbol"];
  n        = sectorData["Dimension"];
  k        = divSectorData["DivergentVariable"];
  ck       = divSectorData["ck"];
  a0       = divSectorData["a0"];
  a1       = divSectorData["a1"];
  aVals    = sectorData["NewExponents"];  (* effective, symbolic in eps *)
  rawAVals = sectorData["RawExponents"];
  polyExps = sectorData["PolynomialExponents"];
  detM     = sectorData["DetM"];
  minExps  = sectorData["MinExponents"];

  kinRules = testKinematics;
  epsRules = {eps -> testEpsilon};
  fullRules = Join[kinRules, epsRules];

  clearedPolys = divSectorData["ClearedPolys"];
  simpPolys    = divSectorData["SimplifiedPolys"];

  yVars = Table[Unique["yv"], {n}];

  (* Original sector integral at finite epsilon.
     Use the cleared polynomial form:
     integrand = |det(M)| * prod y_i^{a_i^eff - 1} * prod Q_j^{B_j}
     where Q_j are cleared polys with non-negative exponents. *)
  Module[{aNum, polyValsExpr, integrand},
    aNum = aVals /. fullRules;
    polyValsExpr = Table[
      Total[
        Table[
          Module[{coeff, exps, logY},
            coeff = mono[[1]] /. kinRules;
            exps  = mono[[2]];
            logY  = Log /@ yVars;
            coeff * Exp[Total[exps * logY]]
          ],
          {mono, clearedPolys[[j]]}
        ]
      ],
      {j, Length[clearedPolys]}
    ];
    integrand = Abs[detM] *
      Exp[Total[(aNum - 1) * Log /@ yVars]] *
      Times @@ MapThread[
        Function[{pv, be}, Exp[be * Log[pv]]],
        {polyValsExpr, polyExps /. fullRules}
      ];

    originalIntegral = Quiet@NIntegrate[
      integrand,
      Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
      MaxRecursion -> 20,
      PrecisionGoal -> 4,
      Method -> "GlobalAdaptive"
    ];
  ];

  B0 = polyExps /. eps -> 0 /. kinRules;
  B1 = (D[polyExps, eps] /. eps -> 0) /. kinRules;

  (* G0: (n-1)-dim integral at eps=0 *)
  Module[{ndVars, g0yVars, g0aNum, g0PolyVals, g0Integrand},
    ndVars    = DeleteCases[Range[n], k];
    g0yVars   = Table[Unique["gy"], {n - 1}];
    g0aNum    = (a0 /. kinRules)[[ndVars]];

    g0PolyVals = Table[
      Total[
        Table[
          Module[{coeff, ndExps, logY},
            coeff  = mono[[1]] /. kinRules;
            ndExps = mono[[2]][[ndVars]];
            logY   = Log /@ g0yVars;
            coeff * Exp[Total[ndExps * logY]]
          ],
          {mono, simpPolys[[j]]}
        ]
      ],
      {j, Length[simpPolys]}
    ];

    g0Integrand = Abs[detM] *
      Exp[Total[(g0aNum - 1) * Log /@ g0yVars]] *
      Times @@ MapThread[
        Function[{pv, be}, Exp[be * Log[pv]]],
        {g0PolyVals, B0}
      ];

    g0Val = Quiet@NIntegrate[
      g0Integrand,
      Evaluate[Sequence @@ ({#, 0, 1} & /@ g0yVars)],
      MaxRecursion -> 20,
      PrecisionGoal -> 4,
      Method -> "GlobalAdaptive"
    ];
  ];

  (* G1: G0 integrand * log insertion sum *)
  Module[{ndVars, g1yVars, g1aNum, g1PolyVals,
          g1BaseIntegrand, logInsertionSum, g1Integrand},
    ndVars  = DeleteCases[Range[n], k];
    g1yVars = Table[Unique["hy"], {n - 1}];
    g1aNum  = (a0 /. kinRules)[[ndVars]];

    g1PolyVals = Table[
      Total[
        Table[
          Module[{coeff, ndExps, logY},
            coeff  = mono[[1]] /. kinRules;
            ndExps = mono[[2]][[ndVars]];
            logY   = Log /@ g1yVars;
            coeff * Exp[Total[ndExps * logY]]
          ],
          {mono, simpPolys[[j]]}
        ]
      ],
      {j, Length[simpPolys]}
    ];

    g1BaseIntegrand = Abs[detM] *
      Exp[Total[(g1aNum - 1) * Log /@ g1yVars]] *
      Times @@ MapThread[
        Function[{pv, be}, Exp[be * Log[pv]]],
        {g1PolyVals, B0}
      ];

    logInsertionSum =
      Total[Table[
        (a1 /. kinRules)[[ndVars[[i]]]] / g1aNum[[i]] * Log[g1yVars[[i]]],
        {i, n - 1}
      ]] +
      Total[Table[
        B1[[j]] * Log[g1PolyVals[[j]]],
        {j, Length[polyExps]}
      ]];

    g1Integrand = g1BaseIntegrand * logInsertionSum;

    g1Val = Quiet@NIntegrate[
      g1Integrand,
      Evaluate[Sequence @@ ({#, 0, 1} & /@ g1yVars)],
      MaxRecursion -> 20,
      PrecisionGoal -> 3,
      Method -> "GlobalAdaptive"
    ];
  ];

  (* Remainder integral *)
  Module[{remYVars, remAnum, fullPolyVals, simpPolyVals,
          bracket, remIntegrand},
    remYVars   = Table[Unique["rv"], {n}];
    remAnum    = a0 /. kinRules;

    fullPolyVals = Table[
      Total[
        Table[
          Module[{coeff, exps, logY},
            coeff = mono[[1]] /. kinRules;
            exps  = mono[[2]];
            logY  = Log /@ remYVars;
            coeff * Exp[Total[exps * logY]]
          ],
          {mono, clearedPolys[[j]]}
        ]
      ],
      {j, Length[clearedPolys]}
    ];

    simpPolyVals = Table[
      Total[
        Table[
          Module[{coeff, exps, logY},
            coeff = mono[[1]] /. kinRules;
            exps  = mono[[2]];
            logY  = Log /@ remYVars;
            coeff * Exp[Total[exps * logY]]
          ],
          {mono, simpPolys[[j]]}
        ]
      ],
      {j, Length[simpPolys]}
    ];

    bracket = Times @@ MapThread[
        Function[{pv, be}, Exp[be * Log[pv]]],
        {fullPolyVals, B0}
      ] -
      Times @@ MapThread[
        Function[{pv, be}, Exp[be * Log[pv]]],
        {simpPolyVals, B0}
      ];

    remIntegrand = Abs[detM] *
      Exp[Total[(remAnum - 1) * Log /@ remYVars]] *
      bracket;

    remVal = Quiet@NIntegrate[
      remIntegrand,
      Evaluate[Sequence @@ ({#, 0, 1} & /@ remYVars)],
      MaxRecursion -> 20,
      PrecisionGoal -> 3,
      Method -> "GlobalAdaptive"
    ];
  ];

  (* Reconstruct: divContrib + remainder = original *)
  divContrib   = g0Val / (ck * testEpsilon) + g1Val / ck;
  reconstructed = divContrib + remVal;
  relError      = Abs[(reconstructed - originalIntegral) / originalIntegral];

  Print["ValidateSubtraction for sector ", sectorData["ConeIndex"], ":"];
  Print["  Original integral (eps=", testEpsilon, "): ", originalIntegral];
  Print["  G0/(ck*eps) + G1/ck = ", divContrib];
  Print["  Remainder = ", remVal];
  Print["  Reconstructed = ", reconstructed];
  Print["  Relative error = ", relError];

  If[relError > 0.01,
    Print["  WARNING: relative error exceeds 1%"]
  ];

  <|"OriginalIntegral" -> originalIntegral,
    "DivergentContribution" -> divContrib,
    "Remainder" -> remVal,
    "Reconstructed" -> reconstructed,
    "RelativeError" -> relError,
    "G0" -> g0Val, "G1" -> g1Val|>
];


(* ============================================================================
   MODULE 3: C++ CODE GENERATION
   ============================================================================ *)

(* --------------------------------------------------------------------------
   MmaToC: Convert Mathematica expression to C++ string
   Recursive pattern-matching converter (not CForm-based).
   -------------------------------------------------------------------------- *)

MmaToC[expr_, paramMap_Association: <||>] := mmaToCInternal[expr, paramMap];

(* Integers -> double literals *)
mmaToCInternal[n_Integer, _] := ToString[n] <> ".0";

(* Rationals -> (a.0/b.0) *)
mmaToCInternal[r_Rational, _] :=
  "(" <> ToString[Numerator[r]] <> ".0/" <>
  ToString[Denominator[r]] <> ".0)";

(* Reals -> string *)
mmaToCInternal[r_Real, _] := ToString[CForm[r]];

(* Complex numbers -> cx(re, im) *)
mmaToCInternal[Complex[re_, im_], pm_] :=
  "cx(" <> mmaToCInternal[re, pm] <> ", " <>
  mmaToCInternal[im, pm] <> ")";

(* Symbols -> parameter lookup or literal *)
mmaToCInternal[s_Symbol, pm_] :=
  If[KeyExistsQ[pm, s],
    pm[s],
    ToString[s]
  ];

(* Subscripted variables like x[i] *)
mmaToCInternal[s_Symbol[i_Integer], pm_] :=
  If[KeyExistsQ[pm, s[i]],
    pm[s[i]],
    ToString[s] <> "[" <> ToString[i] <> "]"
  ];

(* Power *)
mmaToCInternal[Power[base_, -1], pm_] :=
  "(1.0/" <> mmaToCInternal[base, pm] <> ")";

mmaToCInternal[Power[base_, 1/2], pm_] :=
  "std::sqrt(" <> mmaToCInternal[base, pm] <> ")";

mmaToCInternal[Power[base_, -1/2], pm_] :=
  "(1.0/std::sqrt(" <> mmaToCInternal[base, pm] <> "))";

mmaToCInternal[Power[E, exp_], pm_] :=
  "std::exp(" <> mmaToCInternal[exp, pm] <> ")";

mmaToCInternal[Power[base_, n_Integer], pm_] :=
  "std::pow(" <> mmaToCInternal[base, pm] <> ", " <>
  ToString[n] <> ".0)";

mmaToCInternal[Power[base_, exp_], pm_] :=
  "std::pow(" <> mmaToCInternal[base, pm] <> ", " <>
  mmaToCInternal[exp, pm] <> ")";

(* Log *)
mmaToCInternal[Log[arg_], pm_] :=
  "std::log(" <> mmaToCInternal[arg, pm] <> ")";

(* Exp *)
mmaToCInternal[Exp[arg_], pm_] :=
  "std::exp(" <> mmaToCInternal[arg, pm] <> ")";

(* Abs *)
mmaToCInternal[Abs[arg_], pm_] :=
  "std::abs(" <> mmaToCInternal[arg, pm] <> ")";

(* Re, Im *)
mmaToCInternal[Re[arg_], pm_] :=
  "(" <> mmaToCInternal[arg, pm] <> ").real()";

mmaToCInternal[Im[arg_], pm_] :=
  "(" <> mmaToCInternal[arg, pm] <> ").imag()";

(* Plus -> infix with parens.
   HoldPattern is needed because Mathematica evaluates Plus[args__] and
   Times[args__] before storing the DownValue inside BeginPackage. *)
mmaToCInternal[HoldPattern[Plus[args__]], pm_] :=
  "(" <> StringRiffle[mmaToCInternal[#, pm] & /@ {args}, " + "] <> ")";

(* Times: handle negation and general products *)
mmaToCInternal[HoldPattern[Times[-1, rest__]], pm_] :=
  "(-" <> mmaToCInternal[Times[rest], pm] <> ")";

mmaToCInternal[HoldPattern[Times[args__]], pm_] :=
  "(" <> StringRiffle[mmaToCInternal[#, pm] & /@ {args}, " * "] <> ")";

(* Fallback: use CForm and warn *)
mmaToCInternal[expr_, pm_] := Module[{str},
  str = ToString[CForm[expr]];
  str = StringReplace[str, {
    "Power(" ~~ a__ ~~ "," ~~ b__ ~~ ")" :>
      "std::pow(" <> a <> ", " <> b <> ")",
    "Sqrt(" ~~ a__ ~~ ")" :> "std::sqrt(" <> a <> ")"
  }];
  str
];

(* --------------------------------------------------------------------------
   Helper: Generate C++ for one monomial sum (polynomial evaluation)
   -------------------------------------------------------------------------- *)

GenerateMonomialSumCpp[flatPolys_List, polyIndex_Integer,
                       paramMap_Association, dim_Integer,
                       varPrefix_String: "log_y"] :=
Module[{lines, polyVar},
  polyVar = "P" <> ToString[polyIndex];
  lines = {"    cx " <> polyVar <> "(0.0, 0.0);"};

  Do[
    Module[{coeff, alphas, coeffStr, expTerms, expStr},
      coeff    = mono[[1]];
      alphas   = mono[[2]];
      coeffStr = mmaToCInternal[coeff, paramMap];

      expTerms = Table[
        If[TrueQ[alphas[[i]] == 0],
          Nothing,
          mmaToCInternal[alphas[[i]], paramMap] <> " * " <>
            varPrefix <> "[" <> ToString[i - 1] <> "]"
        ],
        {i, dim}
      ];

      expStr = If[Length[expTerms] == 0,
        "0.0",
        StringRiffle[expTerms, " + "]
      ];

      AppendTo[lines,
        "    " <> polyVar <> " += " <> coeffStr <>
        " * std::exp(" <> expStr <> ");"
      ];
    ],
    {mono, flatPolys}
  ];

  StringRiffle[lines, "\n"]
];

(* --------------------------------------------------------------------------
   GenerateCppMonteCarlo
   Main code generation function.
   -------------------------------------------------------------------------- *)

Options[GenerateCppMonteCarlo] = {
  "NSamples"   -> 1000000,
  "MaxDim"     -> 20,
  "SeedBase"   -> 42
};

GenerateCppMonteCarlo[convergentSectors_List, divergentSectors_List,
                      integrandSpec_Association, outputFile_String,
                      OptionsPattern[]] :=
Module[
  {kinSyms, paramMap, nParams,
   code, integrandFuncs, integrandDims, nIntegrands,
   nConvergent, nG0, nG1, nRemainder,
   maxDim, seedBase, nSamples},

  kinSyms  = integrandSpec["KinematicSymbols"];
  nParams  = Length[kinSyms];
  maxDim   = OptionValue["MaxDim"];
  seedBase = OptionValue["SeedBase"];
  nSamples = OptionValue["NSamples"];

  (* Build parameter map: kinematic symbol -> params[i] *)
  paramMap = Association @@ Table[
    kinSyms[[i]] -> ("params[" <> ToString[i - 1] <> "]"),
    {i, nParams}
  ];

  integrandFuncs = {};
  integrandDims  = {};
  nConvergent = 0; nG0 = 0; nG1 = 0; nRemainder = 0;

  (* --- Generate convergent sector integrands --- *)
  Do[
    Module[{sd, flatPolys, polyExps, prefactor, dim, funcName,
            funcCode, polyCode, prodCode},
      sd        = convergentSectors[[s]];
      flatPolys = sd["FlattenedPolys"];
      polyExps  = sd["PolynomialExponents"];
      prefactor = sd["Prefactor"];
      dim       = sd["Dimension"];
      funcName  = "integrand_conv_" <> ToString[s - 1];

      funcCode = "inline cx " <> funcName <>
        "(const double* y, const double* params) {\n";
      funcCode = funcCode <>
        "    // Convergent sector " <> ToString[sd["ConeIndex"]] <> "\n";
      funcCode = funcCode <>
        "    double log_y[" <> ToString[dim] <> "];\n";
      funcCode = funcCode <>
        "    for (int i = 0; i < " <> ToString[dim] <>
        "; i++)\n";
      funcCode = funcCode <>
        "        log_y[i] = (y[i] > 1e-300) ? std::log(y[i]) : -700.0;\n\n";

      Do[
        funcCode = funcCode <>
          GenerateMonomialSumCpp[flatPolys[[j]], j - 1, paramMap, dim] <>
          "\n\n";,
        {j, Length[flatPolys]}
      ];

      funcCode = funcCode <> "    cx result = " <>
        mmaToCInternal[prefactor, paramMap] <> ";\n";

      Do[
        funcCode = funcCode <>
          "    result *= std::exp(" <>
          mmaToCInternal[polyExps[[j]], paramMap] <>
          " * std::log(P" <> ToString[j - 1] <> "));\n";,
        {j, Length[polyExps]}
      ];

      funcCode = funcCode <> "    return result;\n}\n";

      AppendTo[integrandFuncs, funcCode];
      AppendTo[integrandDims, dim];
      nConvergent++;
    ],
    {s, Length[convergentSectors]}
  ];

  (* --- Generate G0 integrands for divergent sectors --- *)
  Do[
    Module[{dd, g0Polys, g0PolyExps, g0Pf, g0Dim, funcName, funcCode},
      dd       = divergentSectors[[s]];
      g0Polys  = dd["G0FlatPolys"];
      g0PolyExps = dd["G0PolyExponents"];
      g0Pf     = dd["G0Prefactor"];
      g0Dim    = dd["G0Dimension"];
      funcName = "integrand_g0_" <> ToString[s - 1];

      funcCode = "inline cx " <> funcName <>
        "(const double* y, const double* params) {\n";
      funcCode = funcCode <>
        "    // G0 for divergent sector " <>
        ToString[dd["ConeIndex"]] <> "\n";
      funcCode = funcCode <>
        "    double log_y[" <> ToString[g0Dim] <> "];\n";
      funcCode = funcCode <>
        "    for (int i = 0; i < " <> ToString[g0Dim] <>
        "; i++)\n";
      funcCode = funcCode <>
        "        log_y[i] = (y[i] > 1e-300) ? std::log(y[i]) : -700.0;\n\n";

      Do[
        funcCode = funcCode <>
          GenerateMonomialSumCpp[g0Polys[[j]], j - 1, paramMap, g0Dim] <>
          "\n\n";,
        {j, Length[g0Polys]}
      ];

      funcCode = funcCode <> "    cx result = " <>
        mmaToCInternal[g0Pf, paramMap] <> ";\n";

      Do[
        funcCode = funcCode <>
          "    result *= std::exp(" <>
          mmaToCInternal[g0PolyExps[[j]], paramMap] <>
          " * std::log(P" <> ToString[j - 1] <> "));\n";,
        {j, Length[g0PolyExps]}
      ];

      funcCode = funcCode <> "    return result;\n}\n";

      AppendTo[integrandFuncs, funcCode];
      AppendTo[integrandDims, g0Dim];
      nG0++;
    ],
    {s, Length[divergentSectors]}
  ];

  (* --- Generate G1 integrands for divergent sectors --- *)
  Do[
    Module[{dd, g0Polys, g0PolyExps, g0Pf, g0Dim, funcName, funcCode,
            logIns, varTerms, polyTerms},
      dd       = divergentSectors[[s]];
      g0Polys  = dd["G0FlatPolys"];
      g0PolyExps = dd["G0PolyExponents"];
      g0Pf     = dd["G0Prefactor"];
      g0Dim    = dd["G0Dimension"];
      logIns   = dd["G1LogInsertions"];
      varTerms  = logIns["VariableTerms"];
      polyTerms = logIns["PolynomialTerms"];
      funcName = "integrand_g1_" <> ToString[s - 1];

      funcCode = "inline cx " <> funcName <>
        "(const double* y, const double* params) {\n";
      funcCode = funcCode <>
        "    // G1 for divergent sector " <>
        ToString[dd["ConeIndex"]] <> "\n";
      funcCode = funcCode <>
        "    double log_y[" <> ToString[g0Dim] <> "];\n";
      funcCode = funcCode <>
        "    for (int i = 0; i < " <> ToString[g0Dim] <>
        "; i++)\n";
      funcCode = funcCode <>
        "        log_y[i] = (y[i] > 1e-300) ? std::log(y[i]) : -700.0;\n\n";

      Do[
        funcCode = funcCode <>
          GenerateMonomialSumCpp[g0Polys[[j]], j - 1, paramMap, g0Dim] <>
          "\n\n";,
        {j, Length[g0Polys]}
      ];

      funcCode = funcCode <> "    cx g0_val = " <>
        mmaToCInternal[g0Pf, paramMap] <> ";\n";

      Do[
        funcCode = funcCode <>
          "    g0_val *= std::exp(" <>
          mmaToCInternal[g0PolyExps[[j]], paramMap] <>
          " * std::log(P" <> ToString[j - 1] <> "));\n";,
        {j, Length[g0PolyExps]}
      ];

      funcCode = funcCode <> "\n    // Log insertion factors\n";
      funcCode = funcCode <> "    cx log_sum(0.0, 0.0);\n";

      Do[
        Module[{coeff, varIdx},
          coeff  = vt[[1]];
          varIdx = vt[[2]] - 1;
          funcCode = funcCode <>
            "    log_sum += " <>
            mmaToCInternal[coeff, paramMap] <>
            " * log_y[" <> ToString[varIdx] <> "];\n";
        ],
        {vt, varTerms}
      ];

      Do[
        Module[{coeff, polyIdx},
          coeff   = pt[[1]];
          polyIdx = pt[[2]] - 1;
          funcCode = funcCode <>
            "    log_sum += " <>
            mmaToCInternal[coeff, paramMap] <>
            " * std::log(P" <> ToString[polyIdx] <> ");\n";
        ],
        {pt, polyTerms}
      ];

      funcCode = funcCode <>
        "\n    return g0_val * log_sum;\n}\n";

      AppendTo[integrandFuncs, funcCode];
      AppendTo[integrandDims, g0Dim];
      nG1++;
    ],
    {s, Length[divergentSectors]}
  ];

  (* --- Generate remainder integrands for divergent sectors --- *)
  Do[
    Module[{dd, rem, fullPolys, simpPolys, remPf, remDim, k,
            divVarExp, polyExpsRem, funcName, funcCode},
      dd        = divergentSectors[[s]];
      rem       = dd["Remainder"];
      fullPolys = rem["FullPolys"];
      simpPolys = rem["SimplifiedPolys"];
      remPf     = rem["Prefactor"];
      remDim    = rem["Dimension"];
      k         = rem["DivVarIndex"] - 1; (* 0-based *)
      divVarExp = rem["DivVarExp"];
      polyExpsRem = rem["PolynomialExponents"];
      funcName = "integrand_rem_" <> ToString[s - 1];

      funcCode = "inline cx " <> funcName <>
        "(const double* y, const double* params) {\n";
      funcCode = funcCode <>
        "    // Remainder for divergent sector " <>
        ToString[dd["ConeIndex"]] <> "\n";
      funcCode = funcCode <>
        "    double log_y[" <> ToString[remDim] <> "];\n";
      funcCode = funcCode <>
        "    for (int i = 0; i < " <> ToString[remDim] <>
        "; i++)\n";
      funcCode = funcCode <>
        "        log_y[i] = (y[i] > 1e-300) ? std::log(y[i]) : -700.0;\n\n";

      funcCode = funcCode <> "    // Full polynomials\n";
      Do[
        funcCode = funcCode <>
          StringReplace[
            GenerateMonomialSumCpp[fullPolys[[j]], j - 1, paramMap,
                                   remDim, "log_y"],
            "P" <> ToString[j - 1] -> "Pfull" <> ToString[j - 1]
          ] <> "\n\n";,
        {j, Length[fullPolys]}
      ];

      funcCode = funcCode <> "    // Simplified polynomials (y_k = 0)\n";
      Do[
        funcCode = funcCode <>
          StringReplace[
            GenerateMonomialSumCpp[simpPolys[[j]], j - 1, paramMap,
                                   remDim, "log_y"],
            "P" <> ToString[j - 1] -> "Psimp" <> ToString[j - 1]
          ] <> "\n\n";,
        {j, Length[simpPolys]}
      ];

      funcCode = funcCode <> "    cx full_prod = ";
      If[Length[polyExpsRem] == 1,
        funcCode = funcCode <>
          "std::exp(" <> mmaToCInternal[polyExpsRem[[1]], paramMap] <>
          " * std::log(Pfull0));\n";,
        funcCode = funcCode <> "cx(1.0, 0.0);\n";
        Do[
          funcCode = funcCode <>
            "    full_prod *= std::exp(" <>
            mmaToCInternal[polyExpsRem[[j]], paramMap] <>
            " * std::log(Pfull" <> ToString[j - 1] <> "));\n";,
          {j, Length[polyExpsRem]}
        ];
      ];

      funcCode = funcCode <> "    cx simp_prod = ";
      If[Length[polyExpsRem] == 1,
        funcCode = funcCode <>
          "std::exp(" <> mmaToCInternal[polyExpsRem[[1]], paramMap] <>
          " * std::log(Psimp0));\n";,
        funcCode = funcCode <> "cx(1.0, 0.0);\n";
        Do[
          funcCode = funcCode <>
            "    simp_prod *= std::exp(" <>
            mmaToCInternal[polyExpsRem[[j]], paramMap] <>
            " * std::log(Psimp" <> ToString[j - 1] <> "));\n";,
          {j, Length[polyExpsRem]}
        ];
      ];

      funcCode = funcCode <>
        "\n    cx yk_factor = std::exp(" <>
        mmaToCInternal[divVarExp - 1, paramMap] <>
        " * log_y[" <> ToString[k] <> "]);\n";

      funcCode = funcCode <>
        "    return " <> mmaToCInternal[remPf, paramMap] <>
        " * yk_factor * (full_prod - simp_prod);\n}\n";

      AppendTo[integrandFuncs, funcCode];
      AppendTo[integrandDims, remDim];
      nRemainder++;
    ],
    {s, Length[divergentSectors]}
  ];

  nIntegrands = Length[integrandFuncs];

  (* --- Assemble the full C++ file --- *)
  code = "// Auto-generated by TropicalEval`GenerateCppMonteCarlo\n";
  code = code <> "// " <> ToString[nConvergent] <> " convergent, " <>
    ToString[nG0] <> " G0, " <> ToString[nG1] <> " G1, " <>
    ToString[nRemainder] <> " remainder integrands\n\n";

  code = code <> "#include <complex>\n";
  code = code <> "#include <cmath>\n";
  code = code <> "#include <random>\n";
  code = code <> "#include <fstream>\n";
  code = code <> "#include <vector>\n";
  code = code <> "#include <iostream>\n";
  code = code <> "#include <string>\n";
  code = code <> "#include <cassert>\n";
  code = code <> "#include <cstdlib>\n";
  code = code <> "#include <array>\n";
  code = code <> "#ifdef _OPENMP\n";
  code = code <> "#include <omp.h>\n";
  code = code <> "#endif\n\n";

  code = code <> "using cx = std::complex<double>;\n\n";

  Do[
    code = code <> integrandFuncs[[i]] <> "\n";,
    {i, nIntegrands}
  ];

  code = code <>
    "// Function pointer type\n" <>
    "using IntegrandFunc = cx(*)(const double*, const double*);\n\n";

  code = code <> "IntegrandFunc integrand_table[] = {\n";
  Module[{allNames},
    allNames = {};
    Do[AppendTo[allNames, "integrand_conv_" <> ToString[i - 1]],
       {i, nConvergent}];
    Do[AppendTo[allNames, "integrand_g0_" <> ToString[i - 1]],
       {i, nG0}];
    Do[AppendTo[allNames, "integrand_g1_" <> ToString[i - 1]],
       {i, nG1}];
    Do[AppendTo[allNames, "integrand_rem_" <> ToString[i - 1]],
       {i, nRemainder}];
    code = code <> "    " <>
      StringRiffle[allNames, ",\n    "] <> "\n";
  ];
  code = code <> "};\n\n";

  code = code <> "int integrand_dim[] = {" <>
    StringRiffle[ToString /@ integrandDims, ", "] <> "};\n";
  code = code <> "const int N_INTEGRANDS = " <>
    ToString[nIntegrands] <> ";\n";
  code = code <> "const int N_PARAMS = " <>
    ToString[nParams] <> ";\n";
  code = code <> "const int MAX_DIM = " <>
    ToString[maxDim] <> ";\n\n";

  (* Main function *)
  code = code <> "int main(int argc, char* argv[]) {\n";
  code = code <> "    if (argc < 3) {\n";
  code = code <> "        std::cerr << \"Usage: \" << argv[0] << \" <input_file> <output_file> [n_samples] [n_threads]\" << std::endl;\n";
  code = code <> "        return 1;\n";
  code = code <> "    }\n\n";

  code = code <> "    std::string input_file = argv[1];\n";
  code = code <> "    std::string output_file = argv[2];\n";
  code = code <> "    int n_samples = (argc > 3) ? std::atoi(argv[3]) : " <>
    ToString[nSamples] <> ";\n";
  code = code <> "    int n_threads = (argc > 4) ? std::atoi(argv[4]) : 1;\n";
  code = code <> "#ifdef _OPENMP\n";
  code = code <> "    if (n_threads == 1) n_threads = omp_get_max_threads();\n";
  code = code <> "    omp_set_num_threads(n_threads);\n";
  code = code <> "#endif\n\n";

  code = code <> "    // Read kinematic data\n";
  code = code <> "    std::ifstream fin(input_file);\n";
  code = code <> "    if (!fin) {\n";
  code = code <> "        std::cerr << \"Cannot open \" << input_file << std::endl;\n";
  code = code <> "        return 1;\n";
  code = code <> "    }\n\n";

  code = code <> "    std::vector<std::vector<double>> kinematic_data;\n";
  code = code <> "    if (N_PARAMS == 0) {\n";
  code = code <> "        // No kinematic parameters: read count from file, default 1\n";
  code = code <> "        int count = 1;\n";
  code = code <> "        fin >> count;\n";
  code = code <> "        if (count < 1) count = 1;\n";
  code = code <> "        for (int i = 0; i < count; i++)\n";
  code = code <> "            kinematic_data.push_back({});\n";
  code = code <> "    } else {\n";
  code = code <> "        double val;\n";
  code = code <> "        std::vector<double> row;\n";
  code = code <> "        while (fin >> val) {\n";
  code = code <> "            row.push_back(val);\n";
  code = code <> "            if ((int)row.size() == N_PARAMS) {\n";
  code = code <> "                kinematic_data.push_back(row);\n";
  code = code <> "                row.clear();\n";
  code = code <> "            }\n";
  code = code <> "        }\n";
  code = code <> "    }\n";
  code = code <> "    fin.close();\n";
  code = code <> "    int n_kp = (int)kinematic_data.size();\n";
  code = code <> "    std::cerr << \"Read \" << n_kp << \" kinematic points\" << std::endl;\n\n";

  code = code <> "    std::vector<std::array<double, 4>> results(n_kp);\n\n";

  code = code <> "    #pragma omp parallel for schedule(dynamic)\n";
  code = code <> "    for (int kp = 0; kp < n_kp; kp++) {\n";
  code = code <> "        const double* params = kinematic_data[kp].data();\n";
  code = code <> "        uint64_t seed = " <> ToString[seedBase] <> "ULL + (uint64_t)kp;\n";
  code = code <> "        std::mt19937_64 rng(seed);\n";
  code = code <> "        std::uniform_real_distribution<double> dist(0.0, 1.0);\n\n";

  code = code <> "        double total_re = 0.0, total_im = 0.0;\n";
  code = code <> "        double total_var_re = 0.0, total_var_im = 0.0;\n\n";

  code = code <> "        for (int s = 0; s < N_INTEGRANDS; s++) {\n";
  code = code <> "            int dim = integrand_dim[s];\n";
  code = code <> "            double mean_re = 0.0, mean_im = 0.0;\n";
  code = code <> "            double M2_re = 0.0, M2_im = 0.0;\n";

  code = code <> "#ifdef TROPICAL_MC_DEBUG\n";
  code = code <> "            int nan_count = 0;\n";
  code = code <> "            double max_mag = 0.0;\n";
  code = code <> "#endif\n\n";

  code = code <> "            for (int k = 0; k < n_samples; k++) {\n";
  code = code <> "                double y[MAX_DIM];\n";
  code = code <> "                for (int i = 0; i < dim; i++) y[i] = dist(rng);\n\n";

  code = code <> "                cx val = integrand_table[s](y, params);\n\n";

  code = code <> "#ifdef TROPICAL_MC_DEBUG\n";
  code = code <> "                if (!std::isfinite(val.real()) || !std::isfinite(val.imag())) {\n";
  code = code <> "                    nan_count++;\n";
  code = code <> "                    if (nan_count <= 5) {\n";
  code = code <> "                        std::cerr << \"NaN/Inf in sector \" << s << \" kp=\" << kp << \" y=[\";";
  code = code <> "\n                        for (int i = 0; i < dim; i++) std::cerr << y[i] << \" \";\n";
  code = code <> "                        std::cerr << \"]\" << std::endl;\n";
  code = code <> "                    }\n";
  code = code <> "                    continue;\n";
  code = code <> "                }\n";
  code = code <> "                double mag = std::abs(val);\n";
  code = code <> "                if (mag > max_mag) max_mag = mag;\n";
  code = code <> "#endif\n\n";

  code = code <> "                double d_re = val.real() - mean_re;\n";
  code = code <> "                mean_re += d_re / (k + 1);\n";
  code = code <> "                M2_re += d_re * (val.real() - mean_re);\n";
  code = code <> "                double d_im = val.imag() - mean_im;\n";
  code = code <> "                mean_im += d_im / (k + 1);\n";
  code = code <> "                M2_im += d_im * (val.imag() - mean_im);\n";
  code = code <> "            }\n\n";

  code = code <> "#ifdef TROPICAL_MC_DEBUG\n";
  code = code <> "            if (kp == 0) {\n";
  code = code <> "                std::cerr << \"Sector \" << s << \": mean=(\" << mean_re << \",\" << mean_im\n";
  code = code <> "                          << \") max_mag=\" << max_mag;\n";
  code = code <> "                if (nan_count > 0)\n";
  code = code <> "                    std::cerr << \" NaN_count=\" << nan_count;\n";
  code = code <> "                double mean_mag = std::sqrt(mean_re*mean_re + mean_im*mean_im);\n";
  code = code <> "                if (mean_mag > 0 && max_mag / mean_mag > 1000)\n";
  code = code <> "                    std::cerr << \" WARNING: large fluctuations\";\n";
  code = code <> "                std::cerr << std::endl;\n";
  code = code <> "            }\n";
  code = code <> "            if ((double)nan_count / n_samples > 0.001)\n";
  code = code <> "                std::cerr << \"WARNING: >0.1%% NaN in sector \" << s << \" kp=\" << kp << std::endl;\n";
  code = code <> "#endif\n\n";

  code = code <> "            total_re += mean_re;\n";
  code = code <> "            total_im += mean_im;\n";
  code = code <> "            total_var_re += M2_re / ((double)n_samples * (n_samples - 1));\n";
  code = code <> "            total_var_im += M2_im / ((double)n_samples * (n_samples - 1));\n";
  code = code <> "        }\n\n";

  code = code <> "        results[kp] = {total_re, total_im,\n";
  code = code <> "                       std::sqrt(total_var_re), std::sqrt(total_var_im)};\n";

  code = code <> "#ifdef TROPICAL_MC_DEBUG\n";
  code = code <> "        if (kp == 0) {\n";
  code = code <> "            std::cerr << \"KP 0 total: (\" << total_re << \", \" << total_im\n";
  code = code <> "                      << \") +/- (\" << std::sqrt(total_var_re) << \", \"\n";
  code = code <> "                      << std::sqrt(total_var_im) << \")\" << std::endl;\n";
  code = code <> "        }\n";
  code = code <> "#endif\n";

  code = code <> "    }\n\n";

  code = code <> "    // Write results\n";
  code = code <> "    std::ofstream fout(output_file);\n";
  code = code <> "    if (!fout) {\n";
  code = code <> "        std::cerr << \"Cannot open \" << output_file << std::endl;\n";
  code = code <> "        return 1;\n";
  code = code <> "    }\n";
  code = code <> "    fout.precision(17);\n";
  code = code <> "    for (int kp = 0; kp < n_kp; kp++) {\n";
  code = code <> "        fout << results[kp][0] << \" \" << results[kp][1] << \" \"\n";
  code = code <> "             << results[kp][2] << \" \" << results[kp][3] << \"\\n\";\n";
  code = code <> "    }\n";
  code = code <> "    fout.close();\n\n";

  code = code <> "    std::cerr << \"Done. Processed \" << n_kp << \" kinematic points.\" << std::endl;\n";
  code = code <> "    return 0;\n";
  code = code <> "}\n";

  Export[outputFile, code, "Text"];

  Print["Generated C++ Monte Carlo code: ", outputFile];
  Print["  ", nConvergent, " convergent sectors"];
  Print["  ", nG0, " G0 integrands"];
  Print["  ", nG1, " G1 integrands"];
  Print["  ", nRemainder, " remainder integrands"];
  Print["  Total: ", nIntegrands, " integrand functions"];

  Module[{codeStr, badPatterns, warnings},
    codeStr = code;
    badPatterns = {"Sin[", "Cos[", "Sqrt[", "Plus[", "Times[",
                   "Power[", "Rule[", "List["};
    warnings = Select[badPatterns, StringContainsQ[codeStr, #] &];
    If[Length[warnings] > 0,
      Print["WARNING: unresolved Mathematica symbols in C++ output: ",
            warnings]
    ];
  ];

  <|"Code" -> code, "OutputFile" -> outputFile,
    "NConvergent" -> nConvergent, "NG0" -> nG0,
    "NG1" -> nG1, "NRemainder" -> nRemainder,
    "NTotal" -> nIntegrands,
    "Dimensions" -> integrandDims|>
];

(* --------------------------------------------------------------------------
   CompileCpp: Compile generated C++ code
   -------------------------------------------------------------------------- *)

CompileCpp[sourceFile_String, outputBinary_String,
           debug_: False] :=
Module[{compiler, flags, cmd, result},
  compiler = "g++";

  flags = If[debug,
    {"-std=c++17", "-O2", "-fopenmp", "-DTROPICAL_MC_DEBUG",
     "-Wall", "-Wextra"},
    {"-std=c++17", "-O3", "-fopenmp", "-Wall", "-Wextra"}
  ];

  cmd = {compiler, Sequence @@ flags, "-o", outputBinary,
         sourceFile, "-lm"};

  Print["Compiling: ", StringRiffle[cmd, " "]];
  result = RunProcess[cmd];

  (* If compilation fails due to -fopenmp (e.g. macOS clang), retry without it *)
  If[result["ExitCode"] != 0 &&
     StringContainsQ[result["StandardError"], "fopenmp"],
    Print["  OpenMP not supported, retrying without -fopenmp..."];
    flags = DeleteCases[flags, "-fopenmp"];
    cmd = {compiler, Sequence @@ flags, "-o", outputBinary,
           sourceFile, "-lm"};
    Print["Compiling: ", StringRiffle[cmd, " "]];
    result = RunProcess[cmd];
  ];

  If[result["ExitCode"] != 0,
    Print["Compilation FAILED:"];
    Print[result["StandardError"]];
    Return[$Failed]
  ];

  If[StringLength[StringTrim[result["StandardError"]]] > 0,
    Print["Compiler warnings:"];
    Print[result["StandardError"]]
  ];

  Print["Compilation successful: ", outputBinary];
  outputBinary
];


(* ============================================================================
   MODULE 4: EvaluateTropicalMC (Driver)
   ============================================================================ *)

Options[EvaluateTropicalMC] = {
  "NSamples"       -> 1000000,
  "NThreads"       -> Automatic,
  "RunChecks"      -> True,
  "EpsilonValue"   -> None,
  "TestEpsilon"    -> 0.01,
  "PrecisionGoal"  -> 3,
  "WorkingDirectory" -> Automatic,
  "Verbose"        -> True
};

EvaluateTropicalMC[integrandSpec_Association, fanData_List,
                   kinematicPoints_List, OptionsPattern[]] :=
Module[
  {dualVertices, simplexList, n, nKP, nParams,
   allSectorData, convergentSectors, divergentSectors,
   processedDivergent, analyticContributions,
   cppFile, cppBinary, kinFile, resultFile,
   cppResult, mcResults, finalResults,
   runChecks, verbose, nSamples, nThreads,
   workDir, epsVal, testEps, precGoal, eps},

  runChecks  = OptionValue["RunChecks"];
  verbose    = OptionValue["Verbose"];
  nSamples   = OptionValue["NSamples"];
  nThreads   = OptionValue["NThreads"];
  workDir    = OptionValue["WorkingDirectory"];
  epsVal     = OptionValue["EpsilonValue"];
  testEps    = OptionValue["TestEpsilon"];
  precGoal   = OptionValue["PrecisionGoal"];
  eps        = integrandSpec["RegulatorSymbol"];

  If[workDir === Automatic,
    workDir = DirectoryName[$InputFileName];
    If[workDir === "" || !StringQ[workDir],
      workDir = Directory[]
    ];
    workDir = FileNameJoin[{workDir, "INTERFILES"}]
  ];
  If[!DirectoryQ[workDir], Quiet[CreateDirectory[workDir]]];

  {dualVertices, simplexList} = fanData;
  n      = Length[integrandSpec["Variables"]];
  nKP    = Length[kinematicPoints];
  nParams = Length[integrandSpec["KinematicSymbols"]];

  If[verbose,
    Print["EvaluateTropicalMC: ", Length[simplexList], " sectors, ",
          nKP, " kinematic points, ", n, " variables"]
  ];

  (* --- Step 1: Validate inputs --- *)
  If[nKP == 0,
    Print["ERROR: no kinematic points provided"];
    Return[$Failed]
  ];

  If[nParams > 0,
    If[!AllTrue[kinematicPoints,
                (ListQ[#] && Length[#] == nParams) &],
      Print["ERROR: kinematic points must have ", nParams, " parameters each"];
      Return[$Failed]
    ];
    If[!AllTrue[Flatten[kinematicPoints], (NumericQ[#] && Im[#] == 0 &&
                                          Abs[#] < 10^15) &],
      Print["ERROR: kinematic values must be real and finite"];
      Return[$Failed]
    ];
  ];

  (* --- Step 2: Process all sectors --- *)
  If[verbose, Print["Processing ", Length[simplexList], " sectors..."]];

  allSectorData = Table[
    Module[{sd, specToUse},
      specToUse = If[epsVal =!= None && eps =!= None,
        MapAt[# /. eps -> epsVal &, integrandSpec,
              {Key["MonomialExponents"]}] //
        MapAt[# /. eps -> epsVal &, #,
              {Key["PolynomialExponents"]}] &,
        integrandSpec
      ];
      sd = ProcessSector[specToUse, dualVertices,
                         simplexList[[s]], s, "Verbose" -> verbose];
      sd
    ],
    {s, Length[simplexList]}
  ];

  If[Count[allSectorData, _Association] != Length[simplexList],
    Print["WARNING: ", Count[allSectorData, $Failed],
          " sectors failed to process"];
  ];

  convergentSectors = Select[allSectorData,
    (AssociationQ[#] && !#["IsDivergent"]) &];
  divergentSectors  = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  If[verbose,
    Print["  ", Length[convergentSectors], " convergent, ",
          Length[divergentSectors], " divergent sectors"]
  ];

  (* --- Step 3: Validation checks --- *)
  If[runChecks && Length[convergentSectors] > 0 && nParams > 0,
    Module[{testKP, kinRules},
      testKP = Take[kinematicPoints, Min[3, nKP]];
      Do[
        kinRules = Thread[
          integrandSpec["KinematicSymbols"] -> testKP[[i]]
        ];
        If[verbose,
          Print["Validating decomposition at kinematic point ", i, "..."]
        ];
        Module[{vr},
          vr = Quiet@ValidateDecomposition[integrandSpec, fanData,
                                           kinRules, precGoal];
          If[AssociationQ[vr],
            Print["  Point ", i, ": rel error = ", vr["RelativeError"]]
          ];
        ];,
        {i, Length[testKP]}
      ];
    ]
  ];

  (* --- Step 4: Process divergent sectors --- *)
  processedDivergent = {};
  analyticContributions = {};

  If[Length[divergentSectors] > 0,
    If[verbose, Print["Processing ", Length[divergentSectors],
                      " divergent sectors..."]];
    processedDivergent = Table[
      ProcessDivergentSector[divergentSectors[[s]], integrandSpec],
      {s, Length[divergentSectors]}
    ];

    analyticContributions = Table[
      Module[{dd, ck},
        dd = processedDivergent[[s]];
        If[!AssociationQ[dd], 0,
          ck = dd["ck"];
          <|"PoleCoeff" -> 1/ck, "FiniteCoeff" -> 1/ck|>
        ]
      ],
      {s, Length[processedDivergent]}
    ];
  ];

  (* --- Step 5: Validate subtraction --- *)
  If[runChecks && Length[processedDivergent] > 0 && nParams > 0,
    Module[{testKP, kinRules},
      testKP = Take[kinematicPoints, Min[2, nKP]];
      Do[
        kinRules = Thread[
          integrandSpec["KinematicSymbols"] -> testKP[[i]]
        ];
        Do[
          If[AssociationQ[processedDivergent[[s]]],
            If[verbose,
              Print["Validating subtraction sector ", s,
                    " at point ", i, "..."]
            ];
            Quiet@ValidateSubtraction[
              processedDivergent[[s]], divergentSectors[[s]],
              integrandSpec, kinRules, testEps
            ];
          ],
          {s, Length[processedDivergent]}
        ],
        {i, Length[testKP]}
      ];
    ]
  ];

  (* --- Step 6: Generate C++ code --- *)
  cppFile    = FileNameJoin[{workDir, "tropical_mc_generated.cpp"}];
  cppBinary  = FileNameJoin[{workDir, "tropical_mc"}];
  kinFile    = FileNameJoin[{workDir, "kinematic_data.txt"}];
  resultFile = FileNameJoin[{workDir, "mc_results.txt"}];

  Module[{convForCpp, divForCpp, specForCpp},
    convForCpp = If[epsVal =!= None && eps =!= None,
      Map[
        MapAt[# /. eps -> epsVal &, #,
              {Key["PolynomialExponents"]}] &,
        convergentSectors
      ],
      convergentSectors
    ];

    divForCpp = If[epsVal =!= None && eps =!= None,
      Map[
        Function[dd,
          If[AssociationQ[dd],
            dd /. eps -> epsVal,
            dd
          ]
        ],
        processedDivergent
      ],
      processedDivergent
    ];

    specForCpp = If[epsVal =!= None && eps =!= None,
      integrandSpec /. eps -> epsVal,
      integrandSpec
    ];

    cppResult = GenerateCppMonteCarlo[
      convForCpp,
      Select[divForCpp, AssociationQ],
      specForCpp, cppFile,
      "NSamples" -> nSamples
    ];
  ];

  If[!AssociationQ[cppResult],
    Print["ERROR: C++ code generation failed"];
    Return[$Failed]
  ];

  (* --- Step 7: Debug compile and test --- *)
  If[runChecks,
    Module[{dbgBinary, dbgResult, dbgKinFile},
      dbgBinary  = FileNameJoin[{workDir, "tropical_mc_dbg"}];
      dbgKinFile = FileNameJoin[{workDir, "kinematic_data_dbg.txt"}];

      If[CompileCpp[cppFile, dbgBinary, True] =!= $Failed,
        Module[{testKP},
          testKP = Take[kinematicPoints, Min[5, nKP]];
          Export[dbgKinFile,
            StringRiffle[
              StringRiffle[ToString[CForm[#]] & /@ #, " "] & /@ testKP,
              "\n"
            ] <> "\n",
            "Text"
          ];

          dbgResult = RunProcess[{dbgBinary, dbgKinFile,
            FileNameJoin[{workDir, "mc_results_dbg.txt"}],
            "100000", "2"}];

          If[dbgResult["ExitCode"] == 0,
            Print["Debug run successful"];
            If[StringLength[StringTrim[dbgResult["StandardError"]]] > 0,
              Print["Debug output:\n", dbgResult["StandardError"]]
            ];,
            Print["Debug run FAILED:"];
            Print[dbgResult["StandardError"]];
          ];
        ];
      ];
    ]
  ];

  (* --- Step 8: Write kinematic data --- *)
  If[nParams > 0,
    Export[kinFile,
      StringRiffle[
        StringRiffle[
          ToString[CForm[#]] & /@ #, " "
        ] & /@ kinematicPoints,
        "\n"
      ] <> "\n",
      "Text"
    ];,
    Export[kinFile,
      StringRiffle[
        Table["0", {nKP}], "\n"
      ] <> "\n",
      "Text"
    ];
  ];

  (* --- Step 9: Release compile and run --- *)
  If[CompileCpp[cppFile, cppBinary, False] === $Failed,
    Print["ERROR: Release compilation failed"];
    Return[$Failed]
  ];

  Module[{nThreadsStr, result},
    nThreadsStr = If[nThreads === Automatic,
      ToString[$ProcessorCount],
      ToString[nThreads]
    ];

    If[verbose, Print["Running Monte Carlo (", nSamples,
                      " samples, ", nThreadsStr, " threads)..."]];

    result = RunProcess[{cppBinary, kinFile, resultFile,
      ToString[nSamples], nThreadsStr}];

    If[result["ExitCode"] != 0,
      Print["ERROR: Monte Carlo execution failed:"];
      Print[result["StandardError"]];
      Return[$Failed]
    ];

    If[verbose && StringLength[StringTrim[result["StandardError"]]] > 0,
      Print[result["StandardError"]]
    ];
  ];

  (* --- Step 10: Read results --- *)
  Module[{rawResults, lines, parsed},
    rawResults = Import[resultFile, "Text"];
    If[rawResults === $Failed,
      Print["ERROR: cannot read results file"];
      Return[$Failed]
    ];

    lines = Select[StringSplit[rawResults, "\n"],
                   StringLength[StringTrim[#]] > 0 &];

    (* Use Read[StringToStream[...], Number] instead of ToExpression
       because C++ outputs scientific notation like 2.05e-05 which
       ToExpression misparses (treats 'e' as a symbol). *)
    parsed = Table[
      Read[StringToStream[#], Number] & /@ StringSplit[line],
      {line, lines}
    ];

    If[Length[parsed] != nKP,
      Print["WARNING: expected ", nKP, " result rows, got ",
            Length[parsed]];
    ];

    Module[{badRows},
      badRows = Select[Range[Length[parsed]],
        !AllTrue[parsed[[#]], NumericQ[#] && Abs[#] < 10^30 &] &
      ];
      If[Length[badRows] > 0,
        Print["WARNING: ", Length[badRows],
              " rows contain non-finite values"]
      ];
    ];

    mcResults = parsed;
  ];

  (* --- Step 11: Assemble final results --- *)
  finalResults = Table[
    If[i <= Length[mcResults] && Length[mcResults[[i]]] >= 4,
      <|"KinematicPoint" -> If[nParams > 0, kinematicPoints[[i]], {}],
        "Re"     -> mcResults[[i, 1]],
        "Im"     -> mcResults[[i, 2]],
        "ReErr"  -> mcResults[[i, 3]],
        "ImErr"  -> mcResults[[i, 4]]|>,
      <|"KinematicPoint" -> If[nParams > 0, kinematicPoints[[i]], {}],
        "Re" -> 0., "Im" -> 0., "ReErr" -> 0., "ImErr" -> 0.|>
    ],
    {i, nKP}
  ];

  (* --- Step 12: Error summary --- *)
  If[verbose,
    Module[{reErrs, imErrs, reMags},
      reErrs = #["ReErr"] & /@ finalResults;
      imErrs = #["ImErr"] & /@ finalResults;
      reMags = Abs[#["Re"]] & /@ finalResults;

      Print["\n=== Monte Carlo Error Summary ==="];
      Print["  Re errors: mean=", Mean[reErrs],
            " median=", Median[reErrs],
            " max=", Max[reErrs]];
      Print["  Im errors: mean=", Mean[imErrs],
            " median=", Median[imErrs],
            " max=", Max[imErrs]];

      Module[{badPts},
        badPts = Select[Range[nKP],
          (reMags[[#]] > 0 &&
           reErrs[[#]] / reMags[[#]] > 0.1) &
        ];
        If[Length[badPts] > 0,
          Print["  WARNING: ", Length[badPts],
                " points have Re error > 10% of result magnitude"]
        ];
      ];
    ]
  ];

  If[runChecks && nParams > 0,
    Module[{testKP, kinRules},
      testKP = Take[kinematicPoints, Min[3, nKP]];
      Print["\n=== NIntegrate Cross-Check ==="];
      Do[
        kinRules = Thread[
          integrandSpec["KinematicSymbols"] -> testKP[[i]]
        ];
        Module[{directResult, mcRe, relErr},
          directResult = Quiet@ValidateDecomposition[
            integrandSpec, fanData, kinRules, precGoal
          ];
          If[AssociationQ[directResult],
            mcRe = finalResults[[i, "Re"]] + I * finalResults[[i, "Im"]];
            relErr = Abs[(mcRe - directResult["DirectResult"]) /
                        directResult["DirectResult"]];
            Print["  Point ", i, ": MC = ", mcRe,
                  ", NIntegrate = ", directResult["DirectResult"],
                  ", rel err = ", relErr];
          ];
        ];,
        {i, Length[testKP]}
      ];
    ]
  ];

  <|"Results"              -> finalResults,
    "ConvergentSectors"    -> Length[convergentSectors],
    "DivergentSectors"     -> Length[divergentSectors],
    "CppFile"              -> cppFile,
    "ResultFile"           -> resultFile,
    "AnalyticContributions" -> analyticContributions|>
];


(* --------------------------------------------------------------------------
   Package end
   -------------------------------------------------------------------------- *)

End[]

EndPackage[]
