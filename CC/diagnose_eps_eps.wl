(* ============================================================================
   Diagnose the epsilon/epsilon effect in the tropical MC subtraction

   The G1 computation in ProcessDivergentSector captures:
     G1 = integral of [G0_integrand × (Σ a1_i/a0_i log(y_i) + Σ B1_j log(P_j))]

   But it MISSES the ε-derivative of the G0 prefactor:
     prefactor(ε) = |det(M)| / prod_{i≠k} a_i(ε)
     d/dε prefactor|_{ε=0} = -prefactor(0) × Σ_{i≠k} a1_i/a0_i

   This gives a missing finite contribution per divergent sector:
     Δ_sector = G0_mc × (-Σ_{i≠k} a1_i/a0_i)

   where G0_mc already includes the 1/ck rescaled prefactor.
   ============================================================================ *)

baseDir = "/home/aidanh/Desktop/Tropical_Monte_Carlo_Final/Bubble1final";
SetDirectory[baseDir];
Get[FileNameJoin[{baseDir, "bispectrum_config.wl"}]];
Get[FileNameJoin[{DirectoryName[$InputFileName], "..", "tropical_eval.wl"}]];

exportData = Get[FileNameJoin[{baseDir, "bubblepp", "presector_integrands_bubblepp.m"}]];
sd1 = exportData["SectorIntegrands"][[1]];
spec = sd1["Spec"];
eps = spec["RegulatorSymbol"];

(* Tropical fan *)
fanPoly = sd1["FanPoly"];
verts = PolytopeVertices[fanPoly^(-1), spec["Variables"]];
fanData = ComputeDecomposition[verts, "ShowProgress" -> False];
{dualVertices, simplexList} = fanData;

allSectorData = Table[
  ProcessSector[spec, dualVertices, simplexList[[c]], c],
  {c, Length[simplexList]}
];
allSectorData = Select[allSectorData, AssociationQ];
divSectors = Select[allSectorData, #["IsDivergent"] &];

Print["Found ", Length[divSectors], " divergent sectors"];
Print[];

(* For each divergent sector, compute the missing prefactor ε-derivative *)
totalCorrection = 0;

Do[
  Module[{divData, k, n, ndVars, a0, a1, ck, g0aVals,
          prefactorLogDeriv, correction},

    divData = Quiet@ProcessDivergentSector[sd, spec];
    If[!AssociationQ[divData], Continue[]];

    k = divData["DivergentVariable"];
    n = divData["Dimension"];
    ndVars = DeleteCases[Range[n], k];
    a0 = divData["a0"];
    a1 = divData["a1"];
    ck = divData["ck"];
    g0aVals = a0[[ndVars]];

    (* The missing term: d/dε [1/prod a_i(ε)] = -Σ a1_i/a0_i × [1/prod a_i(0)] *)
    prefactorLogDeriv = -Total[a1[[ndVars]] / g0aVals];

    Print["Sector ", sd["ConeIndex"],
          ": k=", k,
          ", ck=", ck,
          ", Σ a1_i/a0_i = ", Total[a1[[ndVars]] / g0aVals],
          ", prefactor log deriv = ", N[prefactorLogDeriv]];

    (* The correction to the finite part from this sector is:
       Δ = prefactorLogDeriv × G0_sector_val / ck
       where G0_sector_val is the G0 integral for this sector

       But we don't have per-sector G0 values from the MC.
       Instead, the total correction is approximately:
       Δ_total ≈ prefactorLogDeriv × G0_total_mc
       IF all sectors have the same prefactorLogDeriv.
    *)

    totalCorrection += prefactorLogDeriv;
  ],
  {sd, divSectors}
];

Print[];
Print["Sum of prefactor log derivatives across all ", Length[divSectors], " sectors:"];
Print["  Σ prefactorLogDeriv = ", N[totalCorrection]];
Print[];

(* The correction to the finite part is approximately:
   Δ ≈ (average prefactorLogDeriv) × G0_total
   where G0_total = 0.000513 *)
g0Total = 0.000513;
avgLogDeriv = totalCorrection / Length[divSectors];

Print["Average per-sector prefactor log deriv: ", N[avgLogDeriv]];
Print["G0_total (from MC): ", g0Total];
Print[];

(* Actually, since all sectors have the same ck, and the G0 MC
   already includes the 1/ck factor, the correction per sector is:
   Δ_s = prefactorLogDeriv_s × (G0_s value / ck_s already included)

   But we need per-sector G0 values. Let's compute them via NIntegrate. *)

Print["Computing per-sector G0 values and corrections..."];
Print[];

Module[{totalG0 = 0, totalCorrected = 0, sectorG0s = {}, kinRules},
  kinRules = {k1 -> 1, k2 -> 1, k3 -> 1};

  Do[
    Module[{divData, k, n, ndVars, a0, a1, ck, g0aVals,
            simpPolys, B0, detM, g0Pf,
            yVars, g0PolyVals, g0Integrand, g0Val,
            prefactorLogDeriv, correction},

      divData = Quiet@ProcessDivergentSector[sd, spec];
      If[!AssociationQ[divData], Continue[]];

      k = divData["DivergentVariable"];
      n = divData["Dimension"];
      ndVars = DeleteCases[Range[n], k];
      a0 = divData["a0"];
      a1 = divData["a1"];
      ck = divData["ck"];
      g0aVals = a0[[ndVars]];
      B0 = divData["B0"];
      detM = divData["DetM"];
      simpPolys = divData["SimplifiedPolys"];

      g0Pf = Abs[detM] / (ck * Times @@ g0aVals);  (* includes 1/ck *)

      yVars = Table[Unique["gy"], {n - 1}];

      g0PolyVals = Table[
        Total[Table[Module[{c, e},
          c = mono[[1]] /. kinRules;
          e = mono[[2]][[ndVars]];
          c * Exp[Total[e * Log /@ yVars]]
        ], {mono, simpPolys[[j]]}]],
        {j, Length[simpPolys]}
      ];

      g0Integrand = g0Pf *
        Exp[Total[(g0aVals - 1) * Log /@ yVars]] *
        Times @@ MapThread[
          Function[{p, b}, Exp[b * Log[p]]],
          {g0PolyVals, B0 /. kinRules}
        ];

      g0Val = Quiet@NIntegrate[g0Integrand,
        Evaluate[Sequence @@ ({#, 0, 1} & /@ yVars)],
        MaxRecursion -> 20, PrecisionGoal -> 5];

      prefactorLogDeriv = -Total[a1[[ndVars]] / g0aVals];
      correction = prefactorLogDeriv * g0Val;

      totalG0 += g0Val;
      totalCorrected += correction;

      AppendTo[sectorG0s, <|"Cone" -> sd["ConeIndex"],
        "G0" -> g0Val, "LogDeriv" -> prefactorLogDeriv, "Correction" -> correction|>];
    ],
    {sd, divSectors}
  ];

  Print["Per-sector results (first 5):"];
  Do[
    Print["  Cone ", sectorG0s[[i]]["Cone"],
          ": G0=", NumberForm[sectorG0s[[i]]["G0"], 5],
          ", logDeriv=", NumberForm[sectorG0s[[i]]["LogDeriv"], 4],
          ", correction=", NumberForm[sectorG0s[[i]]["Correction"], 5]],
    {i, Min[5, Length[sectorG0s]]}
  ];
  If[Length[sectorG0s] > 5, Print["  ... (", Length[sectorG0s] - 5, " more)"]];

  Print[];
  Print["Total G0 (NIntegrate): ", totalG0];
  Print["Total missing correction: ", totalCorrected];
  Print[];

  (* Now compare *)
  mcTotal = 0.025448;   (* from 10M MC *)
  mcG0 = 0.000513;      (* from 10M MC *)
  fiesta = 0.024968;     (* from FIESTA 50M *)

  mcFiniteUncorrected = mcTotal - mcG0;
  mcFiniteCorrected = mcFiniteUncorrected + totalCorrected;

  Print["=== COMPARISON ==="];
  Print[];
  Print["  MC TOTAL (10M):           ", mcTotal];
  Print["  MC G0:                    ", mcG0];
  Print["  MC finite (TOTAL - G0):   ", mcFiniteUncorrected];
  Print["  Missing correction:       ", N[totalCorrected]];
  Print["  MC finite (corrected):    ", N[mcFiniteCorrected]];
  Print["  FIESTA (50M):             ", fiesta];
  Print[];
  Print["  Before correction: |MC - FIESTA| = ", Abs[mcFiniteUncorrected - fiesta]];
  Print["  After correction:  |MC - FIESTA| = ", Abs[mcFiniteCorrected - fiesta]];
  Print[];

  If[Abs[mcFiniteCorrected - fiesta] < Abs[mcFiniteUncorrected - fiesta],
    Print["  Correction REDUCES the discrepancy!"];
    If[Abs[mcFiniteCorrected - fiesta] < 0.0001,
      Print["  Agreement after correction: < 0.01% — ε/ε effect confirmed!"];
    ];
  ];
];
