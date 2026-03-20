(* ============================================================================
   Bubblepp Presector 1: FIESTA vs Tropical MC comparison

   Takes presector 1 of the bubble pp diagram (which has UV divergences),
   picks one kinematic point, and compares:
   - Tropical MC: tropical decomposition + C++ MC
   - FIESTA: SDEvaluateDirect on the GL(1)-converted simplex integral

   Both should agree on the 1/eps pole coefficient and finite part.
   ============================================================================ *)

Print["================================================================"];
Print["  Bubblepp Presector 1: FIESTA vs Tropical MC"];
Print["================================================================"];
Print[];

(* ============================================================ *)
(* Load packages                                                 *)
(* ============================================================ *)

baseDir = "/home/aidanh/Desktop/Tropical_Monte_Carlo_Final/Bubble1final";

SetDirectory[baseDir];
Get[FileNameJoin[{baseDir, "bispectrum_config.wl"}]];
Get[FileNameJoin[{baseDir, "bispectrum_utils.wl"}]];
Get[FileNameJoin[{baseDir, "tropical_code_claudef", "tropical_eval.wl"}]];

Print["mu = ", $BispectrumMu, "  delp = ", $BispectrumDelp, "  delm = ", $BispectrumDelm];
Print[];

(* Load presector integrands *)
exportData = Get[FileNameJoin[{baseDir, "bubblepp", "presector_integrands_bubblepp.m"}]];
sectorIntegrands = exportData["SectorIntegrands"];
Print["Loaded ", Length[sectorIntegrands], " presectors."];

(* Pick presector 1 (has divergent sectors) *)
sd1 = sectorIntegrands[[1]];
spec = sd1["Spec"];
prefactor = sd1["Prefactor"];
gamma = sd1["PrefactorLogDeriv"];

Print["Presector 1:"];
Print["  Variables: ", spec["Variables"]];
Print["  Kinematics: ", spec["KinematicSymbols"]];
Print["  Regulator: ", spec["RegulatorSymbol"]];
Print["  # polynomials: ", Length[spec["Polynomials"]]];
Print["  Prefactor: ", N[prefactor]];
Print["  PrefactorLogDeriv (gamma): ", gamma];
Print[];

(* Pick a kinematic point: equilateral k1=k2=k3=1 *)
k1Val = 1; k2Val = 1; k3Val = 1;
kinRules = {k1 -> k1Val, k2 -> k2Val, k3 -> k3Val};
prefNum = prefactor /. kinRules // N;

Print["Kinematic point: k1=", k1Val, ", k2=", k2Val, ", k3=", k3Val];
Print["Prefactor at this point: ", prefNum];
Print[];


(* ============================================================ *)
(* TROPICAL MC: Run tropical decomposition for presector 1       *)
(* ============================================================ *)

Print["================================================================"];
Print["  Tropical MC: Processing presector 1"];
Print["================================================================"];
Print[];

Module[{fanPoly, verts, fanData, dualVertices, simplexList,
        allSectorData, convergentSectors, divergentSectors,
        processedDivergent, eps, specNum,
        mainResult, g0Result},

  fanPoly = sd1["FanPoly"];
  eps = spec["RegulatorSymbol"];

  (* Compute tropical fan *)
  verts = PolytopeVertices[fanPoly^(-1), spec["Variables"]];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  {dualVertices, simplexList} = fanData;

  Print["Fan: ", Length[dualVertices], " rays, ", Length[simplexList], " cones"];

  (* Process all cones *)
  allSectorData = Table[
    ProcessSector[spec, dualVertices, simplexList[[c]], c],
    {c, Length[simplexList]}
  ];
  allSectorData = Select[allSectorData, AssociationQ];

  convergentSectors = Select[allSectorData, !#["IsDivergent"] &];
  divergentSectors = Select[allSectorData, #["IsDivergent"] &];

  Print["Convergent: ", Length[convergentSectors],
        ", Divergent: ", Length[divergentSectors]];

  (* Process divergent sectors *)
  processedDivergent = {};
  If[Length[divergentSectors] > 0,
    processedDivergent = Table[
      Module[{dd},
        dd = ProcessDivergentSector[divergentSectors[[j]], spec];
        If[AssociationQ[dd],
          Module[{ck = dd["ck"]},
            Print["  Div cone ", dd["ConeIndex"], ": ck=", ck];
            dd = <|dd, "G0Prefactor" -> dd["G0Prefactor"] / ck|>;
          ];
        ];
        dd
      ],
      {j, Length[divergentSectors]}
    ];
    processedDivergent = Select[processedDivergent, AssociationQ];
  ];
  Print[];

  (* Generate C++ for main (TOTAL) integral *)
  Module[{convForCpp, specForCpp, fakeG0Sectors,
          mainCpp, g0Cpp, cppResult, mainBin, g0Bin,
          sectorDir, kinFile, resultFile, g0ResultFile,
          runResult},

    sectorDir = FileNameJoin[{baseDir, "tropical_code_claudef", "FIESTA_Compare", "bubblepp_test"}];
    If[!DirectoryQ[sectorDir], CreateDirectory[sectorDir]];

    (* Convergent sectors at eps=0 *)
    convForCpp = Map[
      Function[sd2, <|sd2,
        "FlattenedPolys" -> (sd2["FlattenedPolys"] /. eps -> 0),
        "Prefactor" -> (sd2["Prefactor"] /. eps -> 0),
        "PolynomialExponents" -> (sd2["PolynomialExponents"] /. eps -> 0)|>],
      convergentSectors
    ];

    specForCpp = <|spec, "RegulatorSymbol" -> None|>;

    (* Main C++ *)
    mainCpp = FileNameJoin[{sectorDir, "main_mc.cpp"}];
    cppResult = GenerateCppMonteCarlo[
      convForCpp, processedDivergent, specForCpp, mainCpp,
      "NSamples" -> 2000000
    ];

    (* G0 C++ *)
    g0Cpp = Null;
    If[Length[processedDivergent] > 0,
      fakeG0Sectors = Table[
        Module[{dd = processedDivergent[[j]]},
          <|"FlattenedPolys" -> dd["G0FlatPolys"],
            "PolynomialExponents" -> dd["G0PolyExponents"],
            "Prefactor" -> dd["G0Prefactor"],
            "Dimension" -> dd["G0Dimension"],
            "ConeIndex" -> dd["ConeIndex"],
            "IsDivergent" -> False,
            "DivergentVariable" -> 0|>
        ],
        {j, Length[processedDivergent]}
      ];

      g0Cpp = FileNameJoin[{sectorDir, "g0_mc.cpp"}];
      GenerateCppMonteCarlo[fakeG0Sectors, {}, specForCpp, g0Cpp,
        "NSamples" -> 2000000];
    ];

    (* Compile both *)
    mainBin = FileNameJoin[{sectorDir, "main_mc"}];
    CompileCpp[mainCpp, mainBin, False];

    g0Bin = Null;
    If[g0Cpp =!= Null,
      g0Bin = FileNameJoin[{sectorDir, "g0_mc"}];
      CompileCpp[g0Cpp, g0Bin, False];
    ];

    (* Write kinematic data — single point *)
    kinFile = FileNameJoin[{sectorDir, "kin.txt"}];
    Export[kinFile,
      StringJoin[ToString[CForm[N[k1Val]]], " ",
                 ToString[CForm[N[k2Val]]], " ",
                 ToString[CForm[N[k3Val]]], "\n"],
      "Text"
    ];

    (* Run main *)
    resultFile = FileNameJoin[{sectorDir, "results_main.txt"}];
    Print["Running main MC (2M samples)..."];
    runResult = RunProcess[{mainBin, kinFile, resultFile, "2000000", "4"}];
    If[runResult["ExitCode"] == 0,
      Module[{output, vals},
        output = Import[resultFile, "Text"];
        vals = ToExpression /@ StringSplit[StringTrim[output]];
        mainResult = <|"Re" -> vals[[1]], "Im" -> vals[[2]],
                       "ReErr" -> vals[[3]], "ImErr" -> vals[[4]]|>;
        Print["  Main MC: ", mainResult["Re"], " + ", mainResult["Im"], " I"];
        Print["         ± ", mainResult["ReErr"], " + ", mainResult["ImErr"], " I"];
      ];,
      Print["  Main MC FAILED: ", runResult["StandardError"]];
      mainResult = <|"Re" -> 0., "Im" -> 0., "ReErr" -> 0., "ImErr" -> 0.|>;
    ];

    (* Run G0 *)
    g0Result = <|"Re" -> 0., "Im" -> 0., "ReErr" -> 0., "ImErr" -> 0.|>;
    If[g0Bin =!= Null,
      g0ResultFile = FileNameJoin[{sectorDir, "results_g0.txt"}];
      Print["Running G0 MC (2M samples)..."];
      runResult = RunProcess[{g0Bin, kinFile, g0ResultFile, "2000000", "4"}];
      If[runResult["ExitCode"] == 0,
        Module[{output, vals},
          output = Import[g0ResultFile, "Text"];
          vals = ToExpression /@ StringSplit[StringTrim[output]];
          g0Result = <|"Re" -> vals[[1]], "Im" -> vals[[2]],
                       "ReErr" -> vals[[3]], "ImErr" -> vals[[4]]|>;
          Print["  G0 MC: ", g0Result["Re"], " + ", g0Result["Im"], " I"];
          Print["       ± ", g0Result["ReErr"], " + ", g0Result["ImErr"], " I"];
        ];,
        Print["  G0 MC FAILED: ", runResult["StandardError"]];
      ];
    ];

    Print[];

    (* Combine: finite = pref*(main + (gamma-1)*g0), pole = pref*g0 *)
    Module[{mainVal, g0Val, finiteVal, poleVal,
            mainErrSq, g0ErrSq, finiteErr, poleErr},
      mainVal = mainResult["Re"] + I mainResult["Im"];
      g0Val = g0Result["Re"] + I g0Result["Im"];
      mainErrSq = mainResult["ReErr"]^2 + mainResult["ImErr"]^2;
      g0ErrSq = g0Result["ReErr"]^2 + g0Result["ImErr"]^2;

      finiteVal = prefNum * (mainVal + (gamma - 1) * g0Val);
      poleVal = prefNum * g0Val;
      finiteErr = Abs[prefNum] * Sqrt[mainErrSq + Abs[gamma - 1]^2 * g0ErrSq];
      poleErr = Abs[prefNum] * Sqrt[g0ErrSq];

      Print["--- Tropical MC Results (presector 1 only) ---"];
      Print["  Pole (1/eps coeff):  ", poleVal, " ± ", poleErr];
      Print["  Finite (eps^0):      ", finiteVal, " ± ", finiteErr];
      Print["  Raw main:            ", mainVal];
      Print["  Raw G0:              ", g0Val];
    ];
  ];
];
Print[];


(* ============================================================ *)
(* FIESTA: Same presector integral on the simplex                *)
(* ============================================================ *)

Print["================================================================"];
Print["  FIESTA: Same presector 1 integral (simplex form)"];
Print["================================================================"];
Print[];

(* The presector 1 integral is over [0,inf)^5 with variables x[1]..x[5]:
   f = prod x[i]^{A_i} * prod P_j^{B_j}

   To convert to simplex: introduce x[0], homogenize all polynomials,
   add x[0]^{D-E} where D = total homogeneity degree, E = 6 variables.
   Then use SDEvaluateDirect with delta(1 - sum x[i]).

   However, the polynomials have different degrees, so we need to
   homogenize each separately with appropriate powers of x[0].
*)

(* First, let's compute the structure *)
Module[{polys, monoExps, polyExps, vars, eps, nVars,
        polyDegs, totalDeg, homPolys, homMonoExps, homPolyExps,
        newVars, fiestaFunctions, fiestaDegrees},

  polys = spec["Polynomials"];
  monoExps = spec["MonomialExponents"];
  polyExps = spec["PolynomialExponents"];
  vars = spec["Variables"];
  eps = spec["RegulatorSymbol"];
  nVars = Length[vars];

  Print["Original integral over [0,inf)^", nVars, ":"];
  Print["  Polynomials: ", Length[polys]];
  Do[
    Module[{p, deg},
      p = polys[[j]];
      deg = Max[Table[Total[Exponent[term, #] & /@ vars],
        {term, If[Head[Expand[p]] === Plus, List @@ Expand[p], {Expand[p]}]}]];
      Print["  P", j, ": degree ", deg, " (first terms: ",
            Short[p, 2], ")"];
    ],
    {j, Length[polys]}
  ];

  Print["  Monomial exponents: ", monoExps];
  Print["  Polynomial exponents: ", polyExps];
  Print[];

  (* Compute degrees of each polynomial *)
  polyDegs = Table[
    Module[{p = polys[[j]], terms, degs},
      terms = If[Head[Expand[p]] === Plus, List @@ Expand[p], {Expand[p]}];
      degs = Table[Total[Exponent[term, #] & /@ vars], {term, terms}];
      Max[degs]
    ],
    {j, Length[polys]}
  ];
  Print["Polynomial degrees: ", polyDegs];

  (* Total homogeneity at eps=0:
     sum_i (A_i + 1) + sum_j B_j * deg_j
     should give -D where f(lambda*x) = lambda^{-D} f(x) *)
  Module[{monoExps0, polyExps0, D0},
    monoExps0 = monoExps /. eps -> 0;
    polyExps0 = polyExps /. eps -> 0;
    D0 = -(Total[monoExps0 + 1] + Total[polyExps0 * polyDegs]);
    Print["Total degree D (at eps=0): ", D0];
    Print["D - E = ", D0 - nVars, " (power of x[0] in simplex)"];
    Print[];

    (* For FIESTA: E+1 = nVars+1 variables (x[0]..x[5]) *)
    (* The simplex integral is:
       Int delta(1-x0-...-x5) x0^{D-E-1} prod x[i]^{A_i} prod P_hom_j^{B_j}
       Wait: the GL(1) conversion gives x0^{D-E}, and the measure absorbs one power
       since we have E+1 variables but E dimensional integration after delta.

       Actually: the affine-to-simplex conversion factor is x0^{D-(E+1)+1} = x0^{D-E}.
       But D counts the total scaling degree of the integrand including measure.

       Let me be more careful. The integrand is:
       f(x1,...,x5) = prod x_i^{A_i} * prod P_j(x)^{B_j}
       Under x_i -> lambda*x_i: f -> lambda^S * f where
       S = sum(A_i) + sum(B_j * deg_j)

       The measure d^5 x -> lambda^5 d^5 x.
       Total scaling of integrand*measure: lambda^{S+5}

       For convergence of projective integral: S+5 < 0, i.e. D = -(S+5) > 0.

       Affine chart (x0=1): Int_0^inf d^5 x f(x) [E=5 free variables]
       Simplex chart (sum=1): Int_simplex d^5 x (1-sum)^{D-1} f(x/(sum), ...)

       Hmm, let me use the formula I derived:
       I_affine = Int_simplex d^5 t * t0^{D-(E+1)} * f(t0, t1, ..., t5)|homogenized
       where E+1=6 total variables, D = total homogeneity degree of integrand*measure
       Wait no, D is the degree of the integrand alone (without measure).

       The correct formula: for f homogeneous of degree -D:
       Int_0^inf d^{n-1}x f(1,x1,...) = Int_simplex d^{n-1}t * t0^{D-n} * f(t0,...,t_{n-1})
       where n = total number of projective coordinates = nVars + 1 = 6.
    *)

    (* Homogenize each polynomial *)
    homPolys = Table[
      Module[{p = polys[[j]], deg = polyDegs[[j]], terms, homP, x0},
        x0 = Symbol["x0"];
        terms = If[Head[Expand[p]] === Plus, List @@ Expand[p], {Expand[p]}];
        homP = Total[Table[
          Module[{termDeg},
            termDeg = Total[Exponent[term, #] & /@ vars];
            term * x0^(deg - termDeg)
          ],
          {term, terms}
        ]];
        homP
      ],
      {j, Length[polys]}
    ];

    Print["Homogenized polynomials (first terms):"];
    Do[Print["  P", j, "_hom: ", Short[homPolys[[j]], 3]], {j, Length[homPolys]}];
    Print[];

    (* For FIESTA SDEvaluateDirect:
       Variables: x0 (=x[1] in FIESTA), x[1]..x[5] (=x[2]..x[6] in FIESTA)
       Functions: x0^{D0-nVars-1} * prod x_i^{A_i} and the homogenized polys

       Actually, SDEvaluateDirect takes {functions, degrees, order, deltas}
       where each function is a polynomial and its degree is the power.

       The integrand on the simplex is:
       x0^{D0-6} * prod x_i^{A_i(eps)} * prod P_hom_j^{B_j(eps)}

       Wait, D0 is at eps=0. For eps-dependent D, I need:
       D(eps) = -(sum(A_i(eps)) + sum(B_j(eps)*deg_j))
       and the x0 power is D(eps) - (nVars+1) = D(eps) - 6
    *)

    Module[{Deps, x0Power, fiestaVarMap, x0sym},
      Deps = -(Total[monoExps] + Total[polyExps * polyDegs]);
      x0Power = Simplify[Deps - (nVars + 1)];
      Print["D(eps) = ", Deps];
      Print["x0 power = D - (nVars+1) = ", x0Power];
      Print["x0 power at eps=0: ", x0Power /. eps -> 0];
      Print[];

      (* Now set up FIESTA call *)
      (* Map: x0 -> x[1], x[1]->x[2], ..., x[5]->x[6] *)
      x0sym = x[1];
      fiestaVarMap = Table[vars[[i]] -> x[i + 1], {i, nVars}];

      (* FIESTA functions and degrees:
         {x[1], x[2], x[3], x[4], x[5], x[6], P1_hom, P2_hom, ...}
         {x0Power, A1, A2, A3, A4, A5, B1, B2, ...} *)

      fiestaFunctions = Join[
        {x[1]},  (* x0 *)
        Table[x[i + 1], {i, nVars}],  (* x1..x5 *)
        (homPolys /. Append[fiestaVarMap, Symbol["x0"] -> x[1]])
      ];

      fiestaDegrees = Join[
        {x0Power},
        monoExps,
        polyExps
      ];

      Print["FIESTA setup:"];
      Print["  ", Length[fiestaFunctions], " functions, ",
            nVars + 1, " variables (x[1]..x[", nVars + 1, "])"];
      Print["  Delta: {{1,2,3,4,5,6}}"];
      Print["  Degrees: ", fiestaDegrees /. eps -> 0 // N];
      Print[];

      (* Substitute kinematic values *)
      Module[{fFuncs, fDegs, fResult},
        fFuncs = fiestaFunctions /. kinRules;
        fDegs = fiestaDegrees;

        Print["Running FIESTA SDEvaluateDirect..."];
        Print["  (This may take a while for 6D integral)"];

        SetDirectory["/usr/local/fiesta/FIESTA5"];
        If[!MemberQ[$Packages, "FIESTA`Private`"],
          Get["FIESTA5.m"];
        ];
        SetOptions[FIESTA, "NumberOfSubkernels" -> 4, "NumberOfLinks" -> 8];

        fResult = Quiet@SDEvaluateDirect[
          fFuncs, fDegs, 0,
          {{1, 2, 3, 4, 5, 6}},
          Integrator -> "vegasCuba",
          IntegratorOptions -> {{"maxeval", "5000000"}}
        ];

        Print[];
        Print["  FIESTA result: ", fResult];
        Print[];

        (* Extract pole and finite from FIESTA *)
        Print["--- COMPARISON ---"];
        Print["  FIESTA gives the per-presector integral (no prefactor)."];
        Print["  Tropical MC gives main and G0 (also no prefactor)."];
        Print["  Both should match after accounting for the GL(1) conversion."];
      ];
    ];
  ];
];

Print[];
Print["================================================================"];
Print["  Done"];
Print["================================================================"];
