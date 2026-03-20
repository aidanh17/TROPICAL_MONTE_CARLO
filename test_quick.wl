(* Quick test: just Example 1 from tropical_eval_examples.wl *)
(* Verifies: package loading, Polymake fan computation, ProcessSector, ValidateDecomposition *)

SetDirectory[DirectoryName[$InputFileName]];
Get[FileNameJoin[{Directory[], "tropical_eval.wl"}]];

Print["tropical_eval.wl loaded successfully"];
Print[];

Print["=== Example 1: Basic 2D convergent integral ==="];
Print["Integral[0,Inf] dx1 dx2 / (1 + x1^2 + x2^2)^3"];
Print[];

Module[{poly, vars, spec, verts, fanData, dualVerts, simplices},

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

  Print["Fan: ", Length[dualVerts], " rays, ",
        Length[simplices], " sectors"];

  Do[
    Module[{sd},
      sd = ProcessSector[spec, dualVerts, simplices[[s]], s,
                         "Verbose" -> True];
      If[AssociationQ[sd],
        Print["  Effective exponents: ", sd["NewExponents"]];
        Print["  Prefactor: ", sd["Prefactor"]];
        Print[];
      ];
    ],
    {s, Length[simplices]}
  ];

  Print["Validating against NIntegrate..."];
  Module[{vr},
    vr = Quiet@ValidateDecomposition[spec, fanData, {}, 4];
    Print["  Direct NIntegrate:  ", vr["DirectResult"]];
    Print["  Sector sum:         ", vr["SectorSum"]];
    Print["  Relative error:     ", vr["RelativeError"]];
  ];
];
Print[];
Print["=== Quick test complete ==="];
