(* ============================================================================
   FIESTA vs Tropical MC: Pole AND Finite Part Comparison

   Uses the GL(1) gauge transformation (Cheng-Wu theorem) to convert
   between the [0,inf) representation (tropical MC) and the simplex
   representation (FIESTA).

   For integrand f(x0,x1,...) homogeneous of degree -D in E variables:
     I_affine = Int_0^inf d^{E-1}x f(1,x1,...,x_{E-1})
              = Int_simplex d^{E-1}t  (1-Sum t_i)^{D-E}  f(t0,t1,...,t_{E-1})
   where t0 = 1 - t1 - ... - t_{E-1}.

   This means BOTH codes compute the SAME number.
   ============================================================================ *)

Print["================================================================"];
Print["  FIESTA vs Tropical MC: Pole + Finite Part"];
Print["  via GL(1) Gauge Transformation"];
Print["================================================================"];
Print[];

(* ============================================================================
   TEST 1: Convergent integral  Pi/8

   Tropical MC:  Int_0^inf dx1 dx2 / (1+x1^2+x2^2)^3
   Homogenize:   P_hom = x0^2+x1^2+x2^2,  f = P_hom^{-3}
   Degree:       D = 6, E = 3

   FIESTA (simplex, 3 variables with delta):
     Int delta(1-x0-x1-x2) x0^{D-E} P_hom^{-3}
   = Int delta(1-x0-x1-x2) x0^3 (x0^2+x1^2+x2^2)^{-3}
   ============================================================================ *)

Print["=== Test 1: Convergent 2D integral = Pi/8 ==="];
Print["Tropical MC: Int_0^inf dx1 dx2 / (1+x1^2+x2^2)^3"];
Print["FIESTA:      Int delta(1-x0-x1-x2) x0^3 (x0^2+x1^2+x2^2)^{-3}"];
Print[];

(* --- FIESTA --- *)
SetDirectory["/usr/local/fiesta/FIESTA5"];
Get["FIESTA5.m"];
SetOptions[FIESTA, "NumberOfSubkernels" -> 0, "NumberOfLinks" -> 4];

Module[{result},
  (* x[1]=x0, x[2]=x1, x[3]=x2; delta(1-x[1]-x[2]-x[3]) *)
  result = Quiet@SDEvaluateDirect[
    {x[1], x[1]^2 + x[2]^2 + x[3]^2},
    {3, -3},
    0,
    {{1, 2, 3}}
  ];
  Print["  FIESTA result:      ", result];
  Print["  Expected (Pi/8):    ", N[Pi/8]];
];
Print[];

(* --- Tropical MC --- *)
SetDirectory["/home/aidanh/Desktop/Tropical_Monte_Carlo_Final/Bubble1final/tropical_code_claudef"];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];

Module[{poly, vars, spec, verts, fanData, vr},
  poly = 1 + x[1]^2 + x[2]^2;
  vars = {x[1], x[2]};
  spec = <|"Polynomials" -> {poly}, "MonomialExponents" -> {0, 0},
    "PolynomialExponents" -> {-3}, "Variables" -> vars,
    "KinematicSymbols" -> {}, "RegulatorSymbol" -> None|>;
  verts = PolytopeVertices[poly^(-3), vars];
  fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
  vr = Quiet@ValidateDecomposition[spec, fanData, {}, 5];
  Print["  Tropical MC result: ", vr["SectorSum"]];
  Print["  NIntegrate direct:  ", vr["DirectResult"]];
];
Print[];


(* ============================================================================
   TEST 2: Divergent integral — pole AND finite part

   Tropical MC:  Int_0^inf dx1 dx2 x1^{2ep-1} / [(1+x1)(1+x2)]^2
   Homogenize:   P_hom = (x0+x1)(x0+x2), monomial x1^{2ep-1}
   Full integrand: f(x0,x1,x2) = x1^{2ep-1} * P_hom^{-2}
   Degree: D = 5-2ep, E = 3, D-E = 2-2ep

   FIESTA (simplex):
     Int delta(1-x0-x1-x2) x0^{2-2ep} x1^{2ep-1} [(x0+x1)(x0+x2)]^{-2}

   Exact: Gamma(2ep)*Gamma(2-2ep) = 1/(2ep) - 1 + O(ep)
   So: pole = 1/2, finite = -1
   ============================================================================ *)

Print["=== Test 2: Divergent 2D — pole AND finite ==="];
Print["Tropical MC: Int_0^inf dx1 dx2 x1^{2ep-1}/[(1+x1)(1+x2)]^2"];
Print["FIESTA:      Int delta x0^{2-2ep} x1^{2ep-1} [(x0+x1)(x0+x2)]^{-2}"];
Print["Exact:       Gamma(2ep)*Gamma(2-2ep) = 1/(2ep) - 1 + O(ep)"];
Print[];

(* --- FIESTA --- *)
SetDirectory["/usr/local/fiesta/FIESTA5"];

Module[{result},
  (* x[1]=x0, x[2]=x1, x[3]=x2 *)
  (* Integrand: x[1]^{2-2ep} * x[2]^{2ep-1} * [(x[1]+x[2])*(x[1]+x[3])]^{-2} *)
  result = Quiet@SDEvaluateDirect[
    {x[1], x[2], (x[1] + x[2]) (x[1] + x[3])},
    {2 - 2 ep, 2 ep - 1, -2},
    0,   (* expand to order 0 in ep *)
    {{1, 2, 3}}
  ];
  Print["  FIESTA result: ", result];
];
Print[];

(* --- Tropical MC: full sector decomposition + subtraction --- *)
SetDirectory["/home/aidanh/Desktop/Tropical_Monte_Carlo_Final/Bubble1final/tropical_code_claudef"];

Module[{poly, vars, eps, spec, verts, fanData,
        dualVertices, simplexList,
        allSectorData, divSectors, convSectors,
        totalPole, totalFinite, totalG0, totalG1, totalRemainder,
        totalConvContrib},

  eps = Symbol["epsCmp"];
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
  {dualVertices, simplexList} = fanData;

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData, (AssociationQ[#] && #["IsDivergent"]) &];
  convSectors = Select[allSectorData, (AssociationQ[#] && !#["IsDivergent"]) &];

  Print["  Tropical MC: ", Length[simplexList], " sectors (",
        Length[divSectors], " div, ", Length[convSectors], " conv)"];

  (* Compute convergent sector contributions *)
  totalConvContrib = 0;
  Do[
    Module[{sd, flatPolys, pExps, pf, dim, yVars, polyVals, integrand, result},
      sd = cs;
      flatPolys = sd["FlattenedPolys"];
      pExps = sd["PolynomialExponents"] /. eps -> 0;
      pf = sd["Prefactor"] /. eps -> 0;
      dim = sd["Dimension"];
      yVars = Table[Unique["yv"], {dim}];

      polyVals = Table[
        Total[Table[
          Module[{coeff, alphas, logY},
            coeff = mono[[1]] /. eps -> 0;
            alphas = mono[[2]] /. eps -> 0;
            logY = Log /@ yVars;
            coeff * Exp[Total[alphas * logY]]
          ], {mono, flatPolys[[j]]}]],
        {j, Length[flatPolys]}
      ];

      integrand = pf * Times @@ MapThread[
        Function[{pv, be}, Exp[be * Log[pv]]], {polyVals, pExps}];

      result = Quiet@NIntegrate[integrand,
        Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
        MaxRecursion -> 15, PrecisionGoal -> 5];
      totalConvContrib += result;
    ],
    {cs, convSectors}
  ];

  Print["  Convergent sector sum (eps=0): ", totalConvContrib];

  (* Compute divergent sector contributions: G0, G1, remainder *)
  totalG0 = 0;
  totalG1 = 0;
  totalRemainder = 0;
  totalPole = 0;

  Do[
    Module[{divData, vsResult},
      divData = Quiet@ProcessDivergentSector[sd, spec];
      If[AssociationQ[divData],
        Module[{ck = divData["ck"]},
          vsResult = Quiet@ValidateSubtraction[divData, sd, spec, {}, 0.001];
          If[AssociationQ[vsResult],
            totalG0 += vsResult["G0"] / ck;
            totalG1 += vsResult["G1"] / ck;
            (* Remainder from ValidateSubtraction uses eps=testEpsilon,
               but at eps->0 the remainder converges. Use small eps. *)
          ];

          (* Also compute G0 and G1 analytically via NIntegrate *)
          Module[{g0Val, g1Val, remVal, ndVars, k, g0aVals, g0Polys,
                  g0yVars, g0PolyVals, g0Integrand, g0Pf,
                  g1Base, g1LogSum, g1Integrand,
                  B0, B1, a0, a1, detM, n,
                  clearedPolys, simpPolys,
                  remYVars, fullPolyVals, simpPolyVals, bracket, remIntegrand},

            k = divData["DivergentVariable"];
            n = divData["Dimension"];
            ndVars = DeleteCases[Range[n], k];
            a0 = divData["a0"];
            a1 = divData["a1"];
            B0 = divData["B0"];
            B1 = divData["B1"];
            detM = divData["DetM"];
            clearedPolys = divData["ClearedPolys"];
            simpPolys = divData["SimplifiedPolys"];

            (* G0: (n-1)-dim integral at eps=0 *)
            g0yVars = Table[Unique["gy"], {n - 1}];
            g0aVals = a0[[ndVars]];
            g0Pf = Abs[detM];

            g0PolyVals = Table[
              Total[Table[
                Module[{coeff, ndExps, logY},
                  coeff = mono[[1]];
                  ndExps = mono[[2]][[ndVars]];
                  logY = Log /@ g0yVars;
                  coeff * Exp[Total[ndExps * logY]]
                ], {mono, simpPolys[[j]]}]],
              {j, Length[simpPolys]}
            ];

            g0Integrand = g0Pf *
              Exp[Total[(g0aVals - 1) * Log /@ g0yVars]] *
              Times @@ MapThread[
                Function[{pv, be}, Exp[be * Log[pv]]], {g0PolyVals, B0}];

            g0Val = Quiet@NIntegrate[g0Integrand,
              Evaluate[Sequence @@ ({#, 0, 1} & /@ g0yVars)],
              MaxRecursion -> 20, PrecisionGoal -> 6];

            (* G1: G0 integrand * log insertion *)
            Module[{g1yVars2, g1PolyVals2, g1BaseInt, logIns},
              g1yVars2 = Table[Unique["hy"], {n - 1}];

              g1PolyVals2 = Table[
                Total[Table[
                  Module[{coeff, ndExps, logY},
                    coeff = mono[[1]];
                    ndExps = mono[[2]][[ndVars]];
                    logY = Log /@ g1yVars2;
                    coeff * Exp[Total[ndExps * logY]]
                  ], {mono, simpPolys[[j]]}]],
                {j, Length[simpPolys]}
              ];

              g1BaseInt = g0Pf *
                Exp[Total[(g0aVals - 1) * Log /@ g1yVars2]] *
                Times @@ MapThread[
                  Function[{pv, be}, Exp[be * Log[pv]]], {g1PolyVals2, B0}];

              logIns = Total[Table[
                a1[[ndVars[[i]]]] / g0aVals[[i]] * Log[g1yVars2[[i]]],
                {i, n - 1}
              ]] + Total[Table[B1[[j]] * Log[g1PolyVals2[[j]]],
                {j, Length[B0]}]];

              g1Val = Quiet@NIntegrate[g1BaseInt * logIns,
                Evaluate[Sequence @@ ({#, 0, 1} & /@ g1yVars2)],
                MaxRecursion -> 20, PrecisionGoal -> 5];
            ];

            (* Remainder at eps=0 *)
            remYVars = Table[Unique["rv"], {n}];

            fullPolyVals = Table[
              Total[Table[
                Module[{coeff, exps, logY},
                  coeff = mono[[1]]; exps = mono[[2]];
                  logY = Log /@ remYVars;
                  coeff * Exp[Total[exps * logY]]
                ], {mono, clearedPolys[[j]]}]],
              {j, Length[clearedPolys]}
            ];

            simpPolyVals = Table[
              Total[Table[
                Module[{coeff, exps, logY},
                  coeff = mono[[1]]; exps = mono[[2]];
                  logY = Log /@ remYVars;
                  coeff * Exp[Total[exps * logY]]
                ], {mono, simpPolys[[j]]}]],
              {j, Length[simpPolys]}
            ];

            bracket = Times @@ MapThread[
              Function[{pv, be}, Exp[be * Log[pv]]], {fullPolyVals, B0}] -
              Times @@ MapThread[
              Function[{pv, be}, Exp[be * Log[pv]]], {simpPolyVals, B0}];

            remIntegrand = Abs[detM] *
              Exp[Total[(a0 - 1) * Log /@ remYVars]] * bracket;

            remVal = Quiet@NIntegrate[remIntegrand,
              Evaluate[Sequence @@ ({#, 0, 1} & /@ remYVars)],
              MaxRecursion -> 20, PrecisionGoal -> 4];

            Print["  Sector ", sd["ConeIndex"],
                  ": G0/ck=", N[g0Val/divData["ck"]],
                  ", G1/ck=", N[g1Val/divData["ck"]],
                  ", rem=", N[remVal]];

            totalPole += g0Val / divData["ck"];
            totalG1 += g1Val / divData["ck"];
            totalRemainder += remVal;
          ];
        ];
      ];
    ],
    {sd, divSectors}
  ];

  Print[];
  Print["  --- Tropical MC Laurent expansion ---"];
  Print["  1/eps coefficient (pole):  ", N[totalPole]];
  Print["  eps^0 coefficient (finite): ",
        N[totalG1 + totalRemainder + totalConvContrib]];
  Print["    (G1 contrib:  ", N[totalG1], ")"];
  Print["    (Remainder:   ", N[totalRemainder], ")"];
  Print["    (Conv sectors:", N[totalConvContrib], ")"];
  Print[];
  Print["  --- Exact ---"];
  Print["  Pole:   1/2 = ", N[1/2]];
  Print["  Finite: -1  = ", N[-1]];
];
Print[];


(* ============================================================================
   TEST 3: Non-trivial divergent integral (non-factoring polynomial)

   Tropical MC: Int_0^inf dx1 dx2 x1^{2ep-1} / (1+x1+x2+x1*x2^2)^2
   Homogenize:  P_hom(x0,x1,x2) = x0^3 + x0^2*x1 + x0^2*x2 + x1*x2^2
   (degree 3, D = 7-2ep, E = 3, D-E = 4-2ep)

   FIESTA: Int delta(1-x0-x1-x2) x0^{4-2ep} x1^{2ep-1} P_hom^{-2}
   ============================================================================ *)

Print["=== Test 3: Non-factoring divergent integral ==="];
Print["Tropical MC: Int_0^inf dx1 dx2 x1^{2ep-1}/(1+x1+x2+x1*x2^2)^2"];
Print[];

(* First verify the GL(1) conversion numerically at finite eps *)
Module[{eps0, niAffine, niSimplex, relErr},
  eps0 = 0.05;

  (* Affine (tropical MC domain) *)
  niAffine = Quiet@NIntegrate[
    t1^(2 eps0 - 1) / (1 + t1 + t2 + t1 t2^2)^2,
    {t1, 0, Infinity}, {t2, 0, Infinity},
    MaxRecursion -> 20, PrecisionGoal -> 5
  ];

  (* Simplex (FIESTA domain): homogenize P = x0^3+x0^2*x1+x0^2*x2+x1*x2^2 *)
  (* P(1,x1,x2) = 1+x1+x2+x1*x2^2, check degree: *)
  (* Monomials: 1 (deg 0), x1 (deg 1), x2 (deg 1), x1*x2^2 (deg 3) *)
  (* NOT homogeneous! Max degree = 3, so homogenize to degree 3: *)
  (* P_hom = x0^3 + x0^2*x1 + x0^2*x2 + x1*x2^2 *)
  (* Check: P_hom(1,x1,x2) = 1+x1+x2+x1*x2^2 ✓ *)
  (* D = 2*3 + (2*eps0-1) + 0 + 1 = 6+2*eps0. Wait let me redo. *)
  (* f(x0,x1,x2) = x1^{2ep-1} * P_hom^{-2}, P_hom has degree 3 *)
  (* Under scaling: x1^{2ep-1} -> lam^{2ep-1}, P_hom^{-2} -> lam^{-6} *)
  (* Total: lam^{2ep-7}, so D = 7-2ep, D-E = 4-2ep *)

  niSimplex = Quiet@NIntegrate[
    (1 - s1 - s2)^(4 - 2 eps0) * s1^(2 eps0 - 1) *
    ((1-s1-s2)^3 + (1-s1-s2)^2 s1 + (1-s1-s2)^2 s2 + s1 s2^2)^(-2),
    {s1, 0, 1}, {s2, 0, 1 - s1},
    MaxRecursion -> 20, PrecisionGoal -> 5,
    Method -> "GlobalAdaptive"
  ];

  relErr = Abs[(niAffine - niSimplex) / niAffine];
  Print["  GL(1) conversion check at eps=", eps0, ":"];
  Print["    Affine [0,inf):  ", niAffine];
  Print["    Simplex:         ", niSimplex];
  Print["    Relative error:  ", relErr];
  Print["    MATCH: ", relErr < 0.001];
];
Print[];

(* --- FIESTA --- *)
SetDirectory["/usr/local/fiesta/FIESTA5"];

Module[{result},
  (* P_hom = x[1]^3 + x[1]^2*x[2] + x[1]^2*x[3] + x[2]*x[3]^2 *)
  (* where x[1]=x0, x[2]=x1, x[3]=x2 *)
  result = Quiet@SDEvaluateDirect[
    {x[1], x[2], x[1]^3 + x[1]^2 x[2] + x[1]^2 x[3] + x[2] x[3]^2},
    {4 - 2 ep, 2 ep - 1, -2},
    0,
    {{1, 2, 3}}
  ];
  Print["  FIESTA result: ", result];
];
Print[];

(* --- Tropical MC --- *)
SetDirectory["/home/aidanh/Desktop/Tropical_Monte_Carlo_Final/Bubble1final/tropical_code_claudef"];

Module[{poly, vars, eps, spec, verts, fanData,
        dualVertices, simplexList,
        allSectorData, divSectors, convSectors},

  eps = Symbol["eps3b"];
  poly = 1 + x[1] + x[2] + x[1] x[2]^2;
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
  {dualVertices, simplexList} = fanData;

  allSectorData = Table[
    Quiet@ProcessSector[spec, dualVertices, simplexList[[i]], i],
    {i, Length[simplexList]}
  ];

  divSectors = Select[allSectorData, (AssociationQ[#] && #["IsDivergent"]) &];
  convSectors = Select[allSectorData, (AssociationQ[#] && !#["IsDivergent"]) &];

  Print["  Tropical MC: ", Length[simplexList], " sectors (",
        Length[divSectors], " div, ", Length[convSectors], " conv)"];

  (* For each divergent sector, extract G0, G1, remainder *)
  Module[{totalPole = 0, totalFinite = 0, convContrib = 0},
    (* Convergent sectors at eps=0 *)
    Do[
      Module[{sd, flatPolys, pExps, pf, dim, yVars, polyVals, integrand, result},
        sd = cs;
        flatPolys = sd["FlattenedPolys"];
        pExps = sd["PolynomialExponents"] /. eps -> 0;
        pf = sd["Prefactor"] /. eps -> 0;
        dim = sd["Dimension"];
        yVars = Table[Unique["yv"], {dim}];

        polyVals = Table[
          Total[Table[Module[{c, a, lY},
            c = mono[[1]] /. eps -> 0; a = mono[[2]] /. eps -> 0;
            lY = Log /@ yVars; c * Exp[Total[a * lY]]
          ], {mono, flatPolys[[j]]}]],
          {j, Length[flatPolys]}];

        integrand = pf * Times @@ MapThread[
          Function[{pv, be}, Exp[be * Log[pv]]], {polyVals, pExps}];

        result = Quiet@NIntegrate[integrand,
          Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
          MaxRecursion -> 15, PrecisionGoal -> 4];
        convContrib += result;
      ],
      {cs, convSectors}
    ];

    (* Divergent sectors *)
    Do[
      Module[{divData, ck, k, n, ndVars, a0, a1, B0, B1, detM,
              clearedPolys, simpPolys,
              g0Val, g1Val, remVal},
        divData = Quiet@ProcessDivergentSector[sd, spec];
        If[AssociationQ[divData],
          ck = divData["ck"];
          k = divData["DivergentVariable"];
          n = divData["Dimension"];
          ndVars = DeleteCases[Range[n], k];
          a0 = divData["a0"]; a1 = divData["a1"];
          B0 = divData["B0"]; B1 = divData["B1"];
          detM = divData["DetM"];
          clearedPolys = divData["ClearedPolys"];
          simpPolys = divData["SimplifiedPolys"];

          (* G0 *)
          Module[{yV, aV, pV, intg},
            yV = Table[Unique["g0"], {n-1}];
            aV = a0[[ndVars]];
            pV = Table[Total[Table[Module[{c, e},
              c = mono[[1]]; e = mono[[2]][[ndVars]];
              c * Exp[Total[e * Log /@ yV]]
            ], {mono, simpPolys[[j]]}]], {j, Length[simpPolys]}];
            intg = Abs[detM] * Exp[Total[(aV-1)*Log/@yV]] *
              Times@@MapThread[Function[{p,b},Exp[b*Log[p]]],{pV,B0}];
            g0Val = Quiet@NIntegrate[intg,
              Evaluate[Sequence@@({#,0,1}&/@yV)],
              MaxRecursion->20,PrecisionGoal->5];
          ];

          (* G1 *)
          Module[{yV, aV, pV, base, logI, intg},
            yV = Table[Unique["g1"], {n-1}];
            aV = a0[[ndVars]];
            pV = Table[Total[Table[Module[{c, e},
              c = mono[[1]]; e = mono[[2]][[ndVars]];
              c * Exp[Total[e * Log /@ yV]]
            ], {mono, simpPolys[[j]]}]], {j, Length[simpPolys]}];
            base = Abs[detM] * Exp[Total[(aV-1)*Log/@yV]] *
              Times@@MapThread[Function[{p,b},Exp[b*Log[p]]],{pV,B0}];
            logI = Total[Table[a1[[ndVars[[i]]]]/aV[[i]]*Log[yV[[i]]],
              {i,n-1}]] + Total[Table[B1[[j]]*Log[pV[[j]]],{j,Length[B0]}]];
            intg = base * logI;
            g1Val = Quiet@NIntegrate[intg,
              Evaluate[Sequence@@({#,0,1}&/@yV)],
              MaxRecursion->20,PrecisionGoal->4];
          ];

          (* Remainder *)
          Module[{yV, fPV, sPV, brk, intg},
            yV = Table[Unique["rm"], {n}];
            fPV = Table[Total[Table[Module[{c,e},
              c=mono[[1]]; e=mono[[2]];
              c*Exp[Total[e*Log/@yV]]
            ],{mono,clearedPolys[[j]]}]],{j,Length[clearedPolys]}];
            sPV = Table[Total[Table[Module[{c,e},
              c=mono[[1]]; e=mono[[2]];
              c*Exp[Total[e*Log/@yV]]
            ],{mono,simpPolys[[j]]}]],{j,Length[simpPolys]}];
            brk = Times@@MapThread[Function[{p,b},Exp[b*Log[p]]],{fPV,B0}] -
                  Times@@MapThread[Function[{p,b},Exp[b*Log[p]]],{sPV,B0}];
            intg = Abs[detM]*Exp[Total[(a0-1)*Log/@yV]]*brk;
            remVal = Quiet@NIntegrate[intg,
              Evaluate[Sequence@@({#,0,1}&/@yV)],
              MaxRecursion->20,PrecisionGoal->3];
          ];

          totalPole += g0Val/ck;
          totalFinite += g1Val/ck + remVal;

          Print["  Sector ", sd["ConeIndex"],
                ": pole=", N[g0Val/ck],
                ", G1/ck=", N[g1Val/ck],
                ", rem=", N[remVal]];
        ];
      ],
      {sd, divSectors}
    ];

    totalFinite += convContrib;

    Print[];
    Print["  --- Tropical MC Laurent expansion ---"];
    Print["  Pole (1/eps coeff):  ", N[totalPole]];
    Print["  Finite (eps^0):      ", N[totalFinite]];
    Print["    (conv sectors: ", N[convContrib], ")"];
  ];
];
Print[];


(* ============================================================================
   SUMMARY TABLE
   ============================================================================ *)

Print["================================================================"];
Print["  FINAL COMPARISON TABLE"];
Print["================================================================"];
Print[];
Print["The GL(1) gauge transformation converts between representations:"];
Print["  [0,inf) affine chart  <-->  simplex with delta function"];
Print["Both codes compute the SAME number for each integral."];
Print[];
Print["  Integral        | FIESTA pole | FIESTA finite | TropMC pole | TropMC finite | Exact"];
Print["  ----------------|-------------|---------------|-------------|---------------|------"];
Print["  Pi/8 (conv)     |    n/a      |   (see above) |    n/a      |  (see above)  | Pi/8"];
Print["  (1+x1)(1+x2)   | (see above) |  (see above)  | (see above) | (see above)   | 1/2, -1"];
Print["  Non-factoring   | (see above) |  (see above)  | (see above) | (see above)   | (numerical)"];
Print[];
Print["================================================================"];
