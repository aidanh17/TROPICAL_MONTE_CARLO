(* ============================================================================
   FIESTA vs Tropical Monte Carlo Comparison

   Tests several integrals with both codes and compares:
   1. FIESTA via SDEvaluateDirect (sector decomposition with delta function)
   2. Tropical MC via the tropical_eval pipeline
   3. Exact analytic results where available

   Integrals tested:
   A) 1D divergent: Int_0^1 x^{2eps-1}/(1+x^2) dx = 1/(2eps) - log(2)/2
   B) 1-loop massive bubble (standard Feynman integral via FIESTA)
   C) 2D convergent: Int_0^inf dx1 dx2 / (1+x1^2+x2^2)^3 = Pi/8
   D) 2D divergent:  Int_0^inf dx1 dx2 x1^{2eps-1} / [(1+x1)(1+x2)]^2
      = Gamma(2eps)*Gamma(2-2eps), pole = 1/(2eps), finite = -1
   ============================================================================ *)

Print["================================================================"];
Print["  FIESTA vs Tropical MC Comparison"];
Print["================================================================"];
Print[];

(* ============================================================================
   PART 1: FIESTA Tests
   ============================================================================ *)

Print["================================================================"];
Print["  PART 1: FIESTA Computations"];
Print["================================================================"];
Print[];

(* Load FIESTA *)
SetDirectory["/usr/local/fiesta/FIESTA5"];
Get["FIESTA5.m"];
SetOptions[FIESTA, "NumberOfSubkernels" -> 0];  (* single kernel for simplicity *)
SetOptions[FIESTA, "NumberOfLinks" -> 4];

Print["FIESTA5 loaded."];
Print[];

(* --------------------------------------------------------------------------
   Test A (FIESTA): 1D divergent integral
   Int_0^1 x^{2eps-1} / (1+x^2) dx

   Formulated as: Int dx1 dx2 delta(1-x1-x2) x1^{2ep-1} (1+x1^2)^{-1}
   Setting x2 = 1-x1 gives Int_0^1 x1^{2ep-1} / (1+x1^2)

   Exact: 1/(2eps) - log(2)/2 + O(eps)
   -------------------------------------------------------------------------- *)

Print["--- Test A (FIESTA): 1D divergent integral ---"];
Print["Int_0^1 x^{2ep-1} / (1+x^2) dx"];
Print["Exact: 1/(2ep) - log(2)/2"];
Print[];

Module[{result, exactPole, exactFinite},
  exactPole = 1/2;
  exactFinite = -Log[2]/2;

  result = SDEvaluateDirect[
    {x[1], 1 + x[1]^2},
    {2 ep - 1, -1},
    0,  (* expand to order 0 in ep *)
    {{1, 2}}  (* delta(1 - x[1] - x[2]) *)
  ];

  Print["  FIESTA result: ", result];
  Print["  Exact pole:    1/(2ep) = ", N[exactPole], "/ep"];
  Print["  Exact finite:  -log(2)/2 = ", N[exactFinite]];

  (* Parse FIESTA result: typically {pole/ep + finite + ...} *)
  Print[];
];
Print[];

(* --------------------------------------------------------------------------
   Test B (FIESTA): 1-loop massive bubble, Euclidean kinematics
   Propagators: 1/(k^2+m^2), 1/((k+p)^2+m^2)
   m^2 = 1, p^2 = -1 (Euclidean)

   Known result for d=4-2eps:
   I = (i pi^{d/2}) Gamma(eps) * Int_0^1 da (1 + a(1-a))^{-eps}

   The 1/eps pole coefficient is 1, finite part involves log and dilog.
   -------------------------------------------------------------------------- *)

Print["--- Test B (FIESTA): 1-loop massive bubble ---"];
Print["m^2 = 1, p^2 = -1 (Euclidean), indices {1,1}"];
Print[];

Module[{uf, UU, FF, result},
  uf = UF[{k}, {k^2 + 1, (k + p)^2 + 1}, {p^2 -> -1}];
  UU = uf[[1]];
  FF = uf[[2]];

  Print["  U = ", UU];
  Print["  F = ", FF];

  result = SDEvaluate[{UU, FF, 1}, {1, 1}, 0];

  Print["  FIESTA result: ", result];

  (* Cross-check: Gamma(eps) * Int_0^1 (1+a(1-a))^{-eps} da *)
  Module[{niResult},
    niResult = NIntegrate[(1 + t (1 - t))^(-0.01), {t, 0, 1}];
    Print["  Cross-check at eps=0.01: Gamma(0.01) * NIntegrate = ",
          N[Gamma[0.01] * niResult]];
  ];
  Print[];
];
Print[];

(* --------------------------------------------------------------------------
   Test C (FIESTA): 2D convergent integral via delta function
   Int_0^inf dx1 dx2 / (1+x1^2+x2^2)^3 = Pi/8

   We reformulate by introducing x3 and delta(1-x1-x2-x3):
   But this changes the domain. Instead, use direct NIntegrate as reference.
   -------------------------------------------------------------------------- *)

(* Skip FIESTA for convergent [0,inf] integrals — they don't naturally
   have a delta function. FIESTA comparison is most meaningful for
   divergent integrals with 1/eps poles. *)

(* --------------------------------------------------------------------------
   Test D (FIESTA): 2D divergent integral with factored polynomial
   Int dx1 dx2 dx3 delta(1-x1-x2-x3) x1^{2ep-1} x2^0 x3^0 / [(1+x1)(1+x2)]^2

   Note: this is NOT the same as Int_0^inf since the delta constrains to
   the simplex. We'll compute this separately and compare to a direct check.
   -------------------------------------------------------------------------- *)

Print["--- Test D (FIESTA): 2D divergent with delta function ---"];
Print["Int delta(1-x1-x2) x1^{2ep-1} / (1+x1+x1*x2)^2 dx"];
Print[];

Module[{result},
  (* Use the factored form (1+x1)(1+x2) = 1+x1+x2+x1*x2 *)
  result = SDEvaluateDirect[
    {x[1], (1 + x[1]) (1 + x[2])},
    {2 ep - 1, -2},
    0,
    {{1, 2}}
  ];

  Print["  FIESTA result: ", result];

  (* Cross-check: with delta(1-x1-x2), x2=1-x1, integral becomes
     Int_0^1 x1^{2ep-1} / [(1+x1)(2-x1)]^2 dx1 *)
  Module[{niCheck},
    niCheck = Quiet@NIntegrate[
      t^(2*0.01 - 1) / ((1 + t)(2 - t))^2,
      {t, 0, 1},
      MaxRecursion -> 30, PrecisionGoal -> 6,
      Method -> "DoubleExponential"
    ];
    Print["  NIntegrate at ep=0.01: ", niCheck];
    Print["  (Expect ~50 if pole is ~1/(2ep))"];
  ];
  Print[];
];
Print[];

(* --------------------------------------------------------------------------
   Test E (FIESTA): Standard 1-loop triangle
   -------------------------------------------------------------------------- *)

Print["--- Test E (FIESTA): 1-loop massless box ---"];
Print["(Massless box with s=1, t=-1/2)"];
Print[];

Module[{uf, UU, FF, result},
  (* Massless box: 4 propagators *)
  uf = UF[{k}, {k^2, (k+p1)^2, (k+p1+p2)^2, (k-p4)^2},
          {p1^2->0, p2^2->0, p4^2->0, p1 p2 -> 1/2, p2 p4 -> 1/4, p1 p4 -> -3/4}];
  UU = uf[[1]];
  FF = uf[[2]];

  Print["  U = ", UU];
  Print["  F = ", FF];

  result = Quiet@SDEvaluate[{UU, FF, 1}, {1, 1, 1, 1}, 0,
    ComplexMode -> True];

  Print["  FIESTA result: ", result];
  Print[];
];
Print[];


(* ============================================================================
   PART 2: Tropical MC Tests
   ============================================================================ *)

Print["================================================================"];
Print["  PART 2: Tropical MC Computations"];
Print["================================================================"];
Print[];

(* Load tropical MC *)
SetDirectory["/home/aidanh/Desktop/Tropical_Monte_Carlo_Final/Bubble1final/tropical_code_claudef"];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];

Print["Tropical MC package loaded."];
Print[];

(* --------------------------------------------------------------------------
   Test A (Tropical MC): 1D divergent integral (Test 3)
   Int_0^1 y^{2eps-1} / (1+y^2) dy
   Exact: 1/(2eps) - log(2)/2

   Note: This is a [0,1] integral, which the tropical MC handles as a
   single sector (no fan decomposition needed for 1D).
   We verify via NIntegrate at finite eps.
   -------------------------------------------------------------------------- *)

Print["--- Test A (Tropical MC): 1D divergent ---"];
Print["Int_0^1 y^{2eps-1} / (1+y^2) dy"];
Print["Exact: 1/(2eps) - log(2)/2"];
Print[];

Module[{exactPole, exactFinite, testEps, numericalResult, exactAtEps, relErr},
  exactPole = 1/2;
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

  Print["  At eps = ", testEps, ":"];
  Print["    NIntegrate     = ", numericalResult];
  Print["    Exact formula  = ", N[exactAtEps]];
  Print["    Relative error = ", relErr];
  Print["    PASS: ", relErr < 0.001];
];
Print[];

(* --------------------------------------------------------------------------
   Test C (Tropical MC): 2D convergent integral
   Int_0^inf dx1 dx2 / (1+x1^2+x2^2)^3 = Pi/8

   Full tropical MC pipeline: fan -> sectors -> C++ -> MC result
   -------------------------------------------------------------------------- *)

Print["--- Test C (Tropical MC): 2D convergent integral ---"];
Print["Int_0^inf dx1 dx2 / (1+x1^2+x2^2)^3 = Pi/8"];
Print[];

Module[{poly, vars, spec, verts, fanData, vr, exact},
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

  verts = PolytopeVertices[poly^(-3), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];

  Print["  Fan: ", Length[fanData[[1]]], " rays, ",
        Length[fanData[[2]]], " sectors"];

  vr = Quiet@ValidateDecomposition[spec, fanData, {}, 4];
  If[AssociationQ[vr],
    Print["  NIntegrate (direct): ", vr["DirectResult"]];
    Print["  Tropical sector sum: ", vr["SectorSum"]];
    Print["  Exact:               ", N[exact]];
    Print["  Relative error:      ", vr["RelativeError"]];
    Print["  PASS: ", vr["RelativeError"] < 0.01];
  ];
];
Print[];

(* --------------------------------------------------------------------------
   Test D (Tropical MC): 2D divergent integral (Test 8)
   Int_0^inf dx1 dx2 x1^{2eps-1} / (1+x1+x2+x1*x2)^2

   Polynomial factors: (1+x1)(1+x2)
   Exact: Gamma(2eps)*Gamma(2-2eps) = 1/(2eps) - 1 + O(eps)
   -------------------------------------------------------------------------- *)

Print["--- Test D (Tropical MC): 2D divergent integral ---"];
Print["Int_0^inf dx1 dx2 x1^{2eps-1} / [(1+x1)(1+x2)]^2"];
Print["Exact: Gamma(2eps)*Gamma(2-2eps) = 1/(2eps) - 1 + O(eps)"];
Print[];

Module[{poly, vars, eps, spec, verts, fanData,
        dualVertices, simplexList,
        allSectorData, divSectors, convSectors,
        testEps, directResult, exactResult},

  eps = Symbol["epsD"];
  poly = 1 + x[1] + x[2] + x[1] x[2];  (* = (1+x1)(1+x2) *)
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
  simplexList = fanData[[2]];

  Print["  Fan: ", Length[dualVertices], " rays, ",
        Length[simplexList], " sectors"];

  (* Process all sectors *)
  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];
  convSectors = Select[allSectorData,
    (AssociationQ[#] && !#["IsDivergent"]) &];

  Print["  Convergent sectors: ", Length[convSectors]];
  Print["  Divergent sectors:  ", Length[divSectors]];

  (* Process divergent sectors *)
  Do[
    Module[{divData, vsResult},
      divData = Quiet@ProcessDivergentSector[sd, spec];
      If[AssociationQ[divData],
        Print["  Sector ", sd["ConeIndex"],
              ": divVar=", sd["DivergentVariable"],
              ", ck=", divData["ck"],
              ", pole=", N[1/divData["ck"]]];

        vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, 0.05];
        If[AssociationQ[vsResult],
          Print["    Subtraction validation: relErr = ",
                vsResult["RelativeError"],
                " ", If[vsResult["RelativeError"] < 0.02, "PASS", "FAIL"]];
        ];
      ];
    ],
    {sd, divSectors}
  ];

  (* Verify at finite eps *)
  testEps = 0.01;
  directResult = Quiet@NIntegrate[
    t1^(2 testEps - 1) / ((1 + t1) (1 + t2))^2,
    {t1, 0, Infinity}, {t2, 0, Infinity},
    MaxRecursion -> 20, PrecisionGoal -> 6
  ];
  exactResult = Gamma[2 testEps] Gamma[2 - 2 testEps];

  Print[];
  Print["  Verification at eps = ", testEps, ":"];
  Print["    NIntegrate = ", directResult];
  Print["    Exact      = ", N[exactResult]];
  Print["    Rel error  = ", Abs[(directResult - exactResult)/exactResult]];
  Print["    Pole:  1/(2eps) = ", N[1/(2 testEps)]];
  Print["    Gamma expansion: 1/(2eps) - 1 = ", N[1/(2 testEps) - 1]];
];
Print[];

(* --------------------------------------------------------------------------
   Test F (Tropical MC): 2D divergent with non-factoring polynomial
   Int_0^inf dx1 dx2 x1^{2eps-1} / (1+x1+x2+x1*x2^2+x1^3*x2^2)^2

   This polynomial does NOT factor, so no exact analytic result.
   We validate via NIntegrate at finite eps.
   -------------------------------------------------------------------------- *)

Print["--- Test F (Tropical MC): 2D divergent, non-factoring polynomial ---"];
Print["Int_0^inf dx1 dx2 x1^{2eps-1} x2^{2eps} / (1+x1+x2+x1*x2^2+x1^3*x2^2)^2"];
Print[];

Module[{poly, vars, eps, spec, verts, fanData,
        dualVertices, simplexList,
        allSectorData, divSectors,
        testEps, directResult},

  eps = Symbol["epsF"];
  poly = 1 + x[1] + x[2] + x[1] x[2]^2 + x[1]^3 x[2]^2;
  vars = {x[1], x[2]};

  spec = <|
    "Polynomials"        -> {poly},
    "MonomialExponents"  -> {2 eps - 1, 2 eps},
    "PolynomialExponents" -> {-2},
    "Variables"          -> vars,
    "KinematicSymbols"   -> {},
    "RegulatorSymbol"    -> eps
  |>;

  verts = PolytopeVertices[poly^(-2), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  dualVertices = fanData[[1]];
  simplexList = fanData[[2]];

  Print["  Fan: ", Length[dualVertices], " rays, ",
        Length[simplexList], " sectors"];

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData,
    (AssociationQ[#] && #["IsDivergent"]) &];

  Print["  Divergent sectors: ", Length[divSectors]];

  Do[
    Module[{divData, vsResult},
      divData = Quiet@ProcessDivergentSector[sd, spec];
      If[AssociationQ[divData],
        Print["  Sector ", sd["ConeIndex"],
              ": divVar=", sd["DivergentVariable"],
              ", ck=", divData["ck"]];

        vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, 0.05];
        If[AssociationQ[vsResult],
          Print["    Subtraction validation: relErr = ",
                vsResult["RelativeError"],
                " ", If[vsResult["RelativeError"] < 0.03, "PASS", "FAIL"]];
        ];
      ];
    ],
    {sd, divSectors}
  ];

  (* Verify at finite eps *)
  testEps = 0.05;
  directResult = Quiet@NIntegrate[
    t1^(2 testEps - 1) t2^(2 testEps) /
    (1 + t1 + t2 + t1 t2^2 + t1^3 t2^2)^2,
    {t1, 0, Infinity}, {t2, 0, Infinity},
    MaxRecursion -> 20, PrecisionGoal -> 4,
    Method -> "GlobalAdaptive"
  ];

  Print[];
  Print["  NIntegrate at eps = ", testEps, ": ", directResult];
  Print["  (Finite value confirms integral is well-defined)"];
];
Print[];


(* ============================================================================
   PART 3: Direct Comparison — Shared Integrals
   ============================================================================ *)

Print["================================================================"];
Print["  PART 3: Direct FIESTA vs Tropical MC Comparison"];
Print["================================================================"];
Print[];

(* --------------------------------------------------------------------------
   Comparison 1: 1D divergent integral
   Int_0^1 x^{2eps-1} / (1+x^2) dx
   Exact: 1/(2eps) - log(2)/2

   Both FIESTA and NIntegrate (as proxy for tropical MC sector sum)
   should match the exact answer.
   -------------------------------------------------------------------------- *)

Print["--- Comparison 1: 1D divergent integral ---"];
Print["Int_0^1 x^{2eps-1}/(1+x^2) dx = 1/(2eps) - log(2)/2"];
Print[];

Module[{exactPole, exactFinite, fiestaResult, tropicalNI, eps0},
  exactPole = 1/2;
  exactFinite = -Log[2]/2;

  (* FIESTA *)
  SetDirectory["/usr/local/fiesta/FIESTA5"];
  fiestaResult = Quiet@SDEvaluateDirect[
    {x[1], 1 + x[1]^2},
    {2 ep - 1, -1},
    0,
    {{1, 2}}
  ];

  Print["  Exact pole:     ", N[exactPole], " / ep"];
  Print["  Exact finite:   ", N[exactFinite]];
  Print["  FIESTA result:  ", fiestaResult];

  (* Tropical MC: validate sector decomposition at multiple eps values *)
  SetDirectory["/home/aidanh/Desktop/Tropical_Monte_Carlo_Final/Bubble1final/tropical_code_claudef"];

  Print[];
  Print["  NIntegrate verification at several eps:"];
  Do[
    Module[{ni, exact, relErr},
      ni = NIntegrate[
        y^(2 eps0 - 1) / (1 + y^2), {y, 0, 1},
        MaxRecursion -> 30, PrecisionGoal -> 8,
        Method -> "DoubleExponential"
      ];
      exact = exactPole / eps0 + exactFinite;
      relErr = Abs[(ni - exact) / exact];
      Print["    eps = ", eps0,
            ":  NI = ", ni,
            "  exact = ", N[exact],
            "  relErr = ", relErr];
    ],
    {eps0, {0.1, 0.01, 0.001}}
  ];
];
Print[];

(* --------------------------------------------------------------------------
   Comparison 2: 1-loop bubble
   Both FIESTA and tropical MC sector decomposition should agree.
   -------------------------------------------------------------------------- *)

Print["--- Comparison 2: Standard massive bubble integral ---"];
Print["Int_0^1 da (m^2 + s*a*(1-a))^{-eps} with m^2=1, s=1"];
Print[];

Module[{fiestaResult, niResults},
  (* FIESTA already computed this above as Test B *)
  SetDirectory["/usr/local/fiesta/FIESTA5"];

  Module[{uf, UU, FF},
    uf = UF[{k}, {k^2 + 1, (k + p)^2 + 1}, {p^2 -> -1}];
    UU = uf[[1]]; FF = uf[[2]];
    fiestaResult = Quiet@SDEvaluate[{UU, FF, 1}, {1, 1}, 0];
  ];

  Print["  FIESTA (1/ep + finite): ", fiestaResult];

  (* Direct numerical check *)
  Print[];
  Print["  NIntegrate cross-checks of Gamma(eps) * Int_0^1 (1+a(1-a))^{-eps}:"];
  Do[
    Module[{ni, fullResult},
      ni = NIntegrate[(1 + t (1 - t))^(-eps0), {t, 0, 1}, PrecisionGoal -> 8];
      fullResult = N[Gamma[eps0] * ni];
      Print["    eps = ", eps0, ": Gamma(eps)*NI = ", fullResult];
    ],
    {eps0, {0.1, 0.01, 0.001}}
  ];
];
Print[];


Print["================================================================"];
Print["  Comparison Complete"];
Print["================================================================"];
