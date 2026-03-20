(* ::Package:: *)

(* ============================================================================
   tropical_fan.wl

   Computes the tropical fan (Newton polytope + simplicial decomposition)
   of Euler-type integrands using Polymake.

   Given an integrand of the form  (poly1)^a1 * (poly2)^a2 * ... ,
   this package:
     1. Extracts Newton polytope vertices
     2. Computes the convex hull
     3. Shifts the polytope so the origin is interior
     4. Triangulates the dual fan into simplicial sectors via Polymake

   Dependencies: Polymake  (set $TropicalFanPolymakePath if not auto-detected)

   Usage:
     Get["tropical_fan.wl"]
     verts = PolytopeVertices[integrand, {x[1], x[2]}];
     decomp = ComputeDecomposition[verts];
   ============================================================================ *)

BeginPackage["TropicalFan`"]

(* --- Configuration --- *)

$TropicalFanPolymakePath::usage =
  "$TropicalFanPolymakePath is the path to the polymake executable. \
Auto-detected from common locations or falls back to \"polymake\" on PATH.";

(* --- Public functions --- *)

convexHullVertices::usage =
  "convexHullVertices[points] returns the vertices of the convex hull \
using Polymake. Points is a list of coordinate lists.";

PolytopeVertices::usage =
  "PolytopeVertices[integrand, vars] extracts the Newton polytope vertices \
from an Euler-type integrand (a product of polynomial^exponent terms).";

translateToOriginInteger::usage =
  "translateToOriginInteger[vertices] shifts a polytope so that the origin \
is an interior point.";

ComputeDecomposition::usage =
  "ComputeDecomposition[vertexList] computes the simplicial fan decomposition \
from Newton polytope vertices. Returns {dualVertices, simplexList} where each \
simplex is a list of vertex indices.";

ComputeDecompositiony::usage =
  "ComputeDecompositiony[vertexList] computes the simplicial fan decomposition \
intersected with the positive orthant (y >= 0). Returns {dualVertices, simplexList}.";

(* --- Error messages --- *)

TropicalFan::polymake =
  "Polymake error: `1`. Check that $TropicalFanPolymakePath is set correctly.";

(* ============================================================================
   PRIVATE IMPLEMENTATION
   ============================================================================ *)

Begin["`Private`"]

(* --------------------------------------------------------------------------
   Configuration
   -------------------------------------------------------------------------- *)

If[!ValueQ[$TropicalFanPolymakePath],
  $TropicalFanPolymakePath = Which[
    FileExistsQ["/opt/homebrew/bin/polymake"], "/opt/homebrew/bin/polymake",
    FileExistsQ["/usr/local/bin/polymake"], "/usr/local/bin/polymake",
    FileExistsQ["/usr/bin/polymake"], "/usr/bin/polymake",
    True, "polymake"
  ]
];

(* Unique session ID for temp files *)
sessionID = ToString[$ProcessID] <> "_" <>
  IntegerString[RandomInteger[{100000, 999999}]];

(* --------------------------------------------------------------------------
   Helpers
   -------------------------------------------------------------------------- *)

polymakeFormat[x_] := ToString[x, InputForm];

(* Extract {base, exponent} from a factor *)
extractFactors[Power[base_, exp_]] := {base, exp};
extractFactors[x_] := {x, 1};

(* Parse Euler integrand into list of {base, exponent} pairs *)
parseEulerIntegrand[expr_] := Module[{factors},
  factors = Switch[Head[expr],
    Times, List @@ expr,
    Power, {expr},
    _, {expr}
  ];
  extractFactors /@ factors
];

(* Extract exponent vectors from a polynomial *)
newtonExponents[poly_, vars_List] := Module[{expanded, terms},
  expanded = Expand[poly];
  terms = If[Head[expanded] === Plus, List @@ expanded, {expanded}];
  DeleteDuplicates[
    Table[Exponent[term, v], {term, terms}, {v, vars}]
  ]
];

(* Raise each monomial term to a power separately *)
monomialPower[poly_, n_, vars_List] := Module[
  {expanded, terms, raiseTerm},
  expanded = Expand[poly];
  terms = If[Head[expanded] === Plus, List @@ expanded, {expanded}];
  raiseTerm[term_] := Module[{coeff, exponents},
    coeff = term /. Thread[vars -> 1];
    exponents = Exponent[term, #] & /@ vars;
    coeff^n * Times @@ MapThread[Power, {vars, n * exponents}]
  ];
  Total[raiseTerm /@ terms]
];

monomialPowerAuto[poly_, n_] := monomialPower[poly, n, Variables[poly]];

smallestMultiplier[list_] := LCM @@ (Denominator /@ Flatten[list]);

(* --------------------------------------------------------------------------
   Polymake interface: convex hull
   -------------------------------------------------------------------------- *)

convexHullVertices[points_List] := Module[
  {inFile, outFile, script, result, output, hullVerts},

  inFile = FileNameJoin[{$TemporaryDirectory, "points_" <> sessionID <> ".poly"}];
  outFile = FileNameJoin[{$TemporaryDirectory, "hull_verts_" <> sessionID <> ".txt"}];

  Export[inFile,
    "_type Polytope<Rational>\nPOINTS\n" <>
    StringRiffle[
      StringRiffle[polymakeFormat /@ Prepend[#, 1], " "] & /@ points,
      "\n"
    ] <> "\n",
    "Text"
  ];

  script = "use application 'polytope'; " <>
    "my $p = load('" <> inFile <> "'); " <>
    "open(my $fh, '>', '" <> outFile <> "'); " <>
    "print $fh $p->VERTICES; " <>
    "close($fh);";

  result = RunProcess[{$TropicalFanPolymakePath, script}];
  If[result["ExitCode"] != 0,
    Message[TropicalFan::polymake, result["StandardError"]];
    Return[$Failed]
  ];

  output = Import[outFile, "Text"];
  Quiet[DeleteFile /@ {inFile, outFile}];

  hullVerts = (Rest[#]/First[#]) & /@ (
    ToExpression[StringSplit[#]] & /@
    Select[StringSplit[output, "\n"], StringLength[StringTrim[#]] > 0 &]
  );

  hullVerts
];

(* --------------------------------------------------------------------------
   Polymake interface: interior lattice point
   -------------------------------------------------------------------------- *)

findIntegerPointPolymake[vertices_] := Module[
  {vStr, result, stdout, pt, script},

  vStr = StringRiffle[
    StringRiffle[Prepend[#, 1], ","] & /@ vertices, "],["
  ];

  script = "my $p = new Polytope(POINTS => [[" <> vStr <>
    "]]); print $p->INTERIOR_LATTICE_POINTS->row(0);";

  result = RunProcess[{$TropicalFanPolymakePath, script}];
  stdout = StringTrim[result["StandardOutput"]];

  If[result["ExitCode"] != 0 || stdout == "", $Failed,
    pt = ToExpression /@ StringSplit[stdout];
    Rest[pt]
  ]
];

translateToOriginInteger[vertices_] := Module[{interiorPt},
  interiorPt = findIntegerPointPolymake[vertices];
  If[interiorPt === $Failed,
    interiorPt = findIntegerPointPolymake[
      4 * smallestMultiplier[vertices] * vertices
    ];
    (# - interiorPt) & /@ vertices
    ,
    (# - interiorPt) & /@ vertices
  ]
];

(* --------------------------------------------------------------------------
   Polymake interface: triangulation / fan decomposition
   -------------------------------------------------------------------------- *)

runPolymakeDecomposition[facetVectors_, showProgress_] := Module[
  {dim, polyFile, scriptFile, script, result,
   outputLines, vertexStart, vertexEnd, triangStart, triangEnd,
   nVertices, vertices, nSimplices, simplices,
   process, output, errOutput, startTime,
   progressCell, currentStatus},

  dim = Length[First[facetVectors]] - 1;
  polyFile = FileNameJoin[{$TemporaryDirectory,
    "polytope_input_" <> sessionID <> ".poly"}];
  scriptFile = FileNameJoin[{$TemporaryDirectory,
    "polymake_script_" <> sessionID <> ".pl"}];

  Export[polyFile,
    StringJoin[
      "_type Polytope<Rational>\n\n",
      "INEQUALITIES\n",
      StringRiffle[
        Map[StringRiffle[ToString /@ #, " "] &, facetVectors],
        "\n"
      ],
      "\n"
    ],
    "Text"
  ];

  script = "use application 'polytope';
use strict;
use warnings;

$| = 1;
select(STDERR); $| = 1; select(STDOUT);

print STDERR \"PROGRESS:Loading polytope...\\n\";
my $p = load('" <> polyFile <> "');

print STDERR \"PROGRESS:Computing vertices...\\n\";
my $V = $p->VERTICES;
my $n = $V->rows;
my $d = $V->cols - 1;
print STDERR \"PROGRESS:Found $n vertices in dimension $d\\n\";

print \"VERTICES_START\\n\";
print $n, \"\\n\";
for my $i (0..$n-1) {
    my @coords;
    for my $j (1..$V->cols-1) {
        push @coords, $V->elem($i,$j);
    }
    print join(' ', @coords), \"\\n\";
}
print \"VERTICES_END\\n\";

my @simplices;

print STDERR \"PROGRESS:Checking if polytope is simplicial...\\n\";
if ($p->SIMPLICIAL) {
    print STDERR \"PROGRESS:Polytope is simplicial - extracting facets...\\n\";
    my $vif = $p->VERTICES_IN_FACETS;
    my $nf = $vif->rows;
    for my $i (0..$nf-1) {
        if ($i % 100 == 0 || $i == $nf-1) {
            my $pct = int(100 * ($i+1) / $nf);
            print STDERR \"PROGRESS:Processing facet \".($i+1).\" of $nf ($pct%)\\n\";
        }
        my @facet = @{$vif->row($i)};
        push @simplices, \\@facet;
    }
} else {
    print STDERR \"PROGRESS:Polytope is non-simplicial - triangulating facets...\\n\";
    my $facets = $p->FACETS;
    my $vif = $p->VERTICES_IN_FACETS;
    my $nf = $vif->rows;

    for my $i (0..$nf-1) {
        my $pct = int(100 * ($i+1) / $nf);
        print STDERR \"PROGRESS:Triangulating facet \".($i+1).\" of $nf ($pct%)\\n\";

        my @facet_verts = @{$vif->row($i)};

        if (scalar(@facet_verts) == $d) {
            push @simplices, \\@facet_verts;
        } else {
            my @pts;
            for my $vi (@facet_verts) {
                push @pts, $V->row($vi);
            }
            my $facet_poly = new Polytope<Rational>(POINTS => \\@pts);
            my $triang = $facet_poly->TRIANGULATION;
            my $tf = $triang->FACETS;
            my $num_simplices = scalar(@{$tf});

            for my $j (0..$num_simplices-1) {
                my @local_simplex = @{$tf->[$j]};
                my @global_simplex = map { $facet_verts[$_] } @local_simplex;
                push @simplices, \\@global_simplex;
            }
        }
    }
}

print STDERR \"PROGRESS:Writing output...\\n\";
print \"TRIANGULATION_START\\n\";
print scalar(@simplices), \"\\n\";
for my $s (@simplices) {
    my @full_simplex = (0, map { $_ + 1 } @{$s});
    print join(' ', @full_simplex), \"\\n\";
}
print \"TRIANGULATION_END\\n\";
print STDERR \"PROGRESS:Done!\\n\";
";

  Export[scriptFile, script, "Text"];
  startTime = AbsoluteTime[];

  If[showProgress,
    (* Interactive mode with progress indicator *)
    process = StartProcess[{$TropicalFanPolymakePath, "--script", scriptFile}];
    currentStatus = "Starting polymake...";

    progressCell = PrintTemporary[
      Dynamic[
        Row[{
          ProgressIndicator[Appearance -> "Indeterminate"],
          "  ", currentStatus,
          "  [", ToString[Round[AbsoluteTime[] - startTime, 0.1]], "s]"
        }],
        UpdateInterval -> 0.5
      ]
    ];

    While[ProcessStatus[process] === "Running",
      errOutput = ReadString[ProcessConnection[process, "StandardError"], EndOfBuffer];
      If[StringQ[errOutput] && StringLength[errOutput] > 0,
        With[{lines = StringSplit[errOutput, "\n"]},
          Do[
            If[StringStartsQ[line, "PROGRESS:"],
              currentStatus = StringDrop[line, 9]
            ],
            {line, lines}
          ]
        ]
      ];
      Pause[0.1];
    ];

    output = ReadString[ProcessConnection[process, "StandardOutput"], EndOfFile];
    errOutput = ReadString[ProcessConnection[process, "StandardError"], EndOfFile];
    NotebookDelete[progressCell];

    If[ProcessInformation[process]["ExitCode"] =!= 0,
      Print["Polymake error after ", Round[AbsoluteTime[] - startTime, 0.1], "s:"];
      Print[errOutput];
      KillProcess[process];
      Return[$Failed]
    ];
    KillProcess[process];
    Print["Completed in ", Round[AbsoluteTime[] - startTime, 0.1], " seconds"];
    ,
    (* Non-interactive mode *)
    result = RunProcess[{$TropicalFanPolymakePath, "--script", scriptFile}];
    If[result["ExitCode"] != 0,
      Print["Polymake error: ", result["StandardError"]];
      Return[$Failed]
    ];
    output = result["StandardOutput"];
  ];

  (* Parse output *)
  outputLines = StringSplit[output, "\n"];

  vertexStart = FirstPosition[outputLines, "VERTICES_START"][[1]] + 1;
  vertexEnd = FirstPosition[outputLines, "VERTICES_END"][[1]] - 1;
  nVertices = ToExpression[outputLines[[vertexStart]]];
  vertices = Table[
    ToExpression /@ StringSplit[outputLines[[vertexStart + i]]],
    {i, 1, nVertices}
  ];

  triangStart = FirstPosition[outputLines, "TRIANGULATION_START"][[1]] + 1;
  triangEnd = FirstPosition[outputLines, "TRIANGULATION_END"][[1]] - 1;
  nSimplices = ToExpression[outputLines[[triangStart]]];
  simplices = Table[
    ToExpression /@ StringSplit[outputLines[[triangStart + i]]],
    {i, 1, nSimplices}
  ];

  simplices = DeleteCases[#, 0] & /@ simplices;
  {dim, vertices, simplices}
];

(* --------------------------------------------------------------------------
   High-level decomposition functions
   -------------------------------------------------------------------------- *)

Options[ComputeDecomposition] = {"ShowProgress" -> True};

ComputeDecomposition[vertex_List, OptionsPattern[]] := Module[
  {hullVerts, shiftedVerts, facetVectors, result, vertices, simplices},

  hullVerts = convexHullVertices[vertex];
  If[hullVerts === $Failed, Return[$Failed]];

  shiftedVerts = translateToOriginInteger[hullVerts];
  If[shiftedVerts === $Failed, Return[$Failed]];

  facetVectors = Prepend[-#, 1] & /@ shiftedVerts;
  result = runPolymakeDecomposition[facetVectors, OptionValue["ShowProgress"]];
  If[result === $Failed, Return[$Failed]];

  {vertices, simplices} = result[[{2, 3}]];
  {vertices, simplices}
];

Options[ComputeDecompositiony] = {"ShowProgress" -> True};

ComputeDecompositiony[vertex_List, OptionsPattern[]] := Module[
  {hullVerts, shiftedVerts, facetVectorsPre, facetVectorsOrth,
   facetVectors, result, dim, vertices, simplices, badvertex},

  hullVerts = convexHullVertices[vertex];
  If[hullVerts === $Failed, Return[$Failed]];

  shiftedVerts = translateToOriginInteger[hullVerts];
  If[shiftedVerts === $Failed, Return[$Failed]];

  facetVectorsPre = Prepend[-#, 1] & /@ shiftedVerts;
  facetVectorsOrth = Table[
    ReplacePart[ConstantArray[0, Length[shiftedVerts[[1]]] + 1], i + 1 -> -1],
    {i, 1, Length[shiftedVerts[[1]]]}
  ];
  facetVectors = Join[facetVectorsPre, facetVectorsOrth];

  result = runPolymakeDecomposition[facetVectors, OptionValue["ShowProgress"]];
  If[result === $Failed, Return[$Failed]];

  {dim, vertices, simplices} = result;
  badvertex = Position[vertices, ConstantArray[0, dim]][[1, 1]];
  simplices = Select[simplices, !MemberQ[#, badvertex] &];

  {vertices, simplices}
];

(* --------------------------------------------------------------------------
   Newton polytope extraction from integrand
   -------------------------------------------------------------------------- *)

PolytopeVertices[integrand_, vars_] := Module[{polypre},
  polypre = {#[[1]], If[MemberQ[vars, #[[1]]], 0, 1]} & /@
    parseEulerIntegrand[integrand];
  newtonExponents[
    Expand[Times @@ (monomialPowerAuto[#[[1]], Re[#[[2]]]] & /@ polypre)],
    vars
  ]
];

End[]

EndPackage[]
