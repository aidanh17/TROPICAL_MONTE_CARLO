(* ============================================================================
   tropical_fan_examples.wl

   Examples demonstrating the tropical_fan package.
   Requires: Polymake installed and accessible.

   Load the package first:
     SetDirectory[NotebookDirectory[]];   (* or wherever tropical_fan.wl lives *)
     << tropical_fan`
   ============================================================================ *)

(* --- Load package --- *)
(* If running from a notebook, use: SetDirectory[NotebookDirectory[]] first *)
(* If running as a script, use the full path below *)
Get[FileNameJoin[{DirectoryName[$InputFileName], "tropical_fan.wl"}]];

(* Optional: set polymake path if not auto-detected *)
(* $TropicalFanPolymakePath = "/opt/homebrew/bin/polymake"; *)


(* ============================================================================
   Example 1: One-loop bubble integral (2 variables)

   Integrand:  (x[1] + x[2])^(-1) * (x[1] + 2 x[2])^(-1)
   ============================================================================ *)

Print["=== Example 1: One-loop bubble ==="];

integrand1 = (x[1] + x[2])^(-1) * (x[1] + 2 x[2])^(-1);
vars1 = {x[1], x[2]};

(* Step 1: Extract Newton polytope vertices *)
verts1 = PolytopeVertices[integrand1, vars1];
Print["Newton polytope vertices: ", verts1];

(* Step 2: Compute the simplicial decomposition (tropical fan) *)
decomp1 = ComputeDecomposition[verts1];
Print["Dual fan vertices: ", decomp1[[1]]];
Print["Simplices (sectors): ", decomp1[[2]]];
Print["Number of sectors: ", Length[decomp1[[2]]]];
Print[];


(* ============================================================================
   Example 2: Triangle integral (3 variables)

   Integrand:  (x[1] + x[2] + x[3])^(-2) * (x[1] x[2] + x[2] x[3] + x[1] x[3])^(-1)
   ============================================================================ *)

Print["=== Example 2: Triangle integral ==="];

integrand2 = (x[1] + x[2] + x[3])^(-2) *
             (x[1] x[2] + x[2] x[3] + x[1] x[3])^(-1);
vars2 = {x[1], x[2], x[3]};

verts2 = PolytopeVertices[integrand2, vars2];
Print["Newton polytope vertices: ", verts2];

decomp2 = ComputeDecomposition[verts2];
Print["Dual fan vertices: ", decomp2[[1]]];
Print["Number of sectors: ", Length[decomp2[[2]]]];
Print[];


(* ============================================================================
   Example 3: Simple 2D polytope from explicit points

   Directly provide points and compute the fan, bypassing
   integrand parsing. Useful when you already know the polytope.
   ============================================================================ *)

Print["=== Example 3: Direct polytope input ==="];

points3 = {{0, 0}, {2, 0}, {0, 2}, {1, 1}};

(* Convex hull *)
hull3 = convexHullVertices[points3];
Print["Convex hull vertices: ", hull3];

(* Shift so origin is interior *)
shifted3 = translateToOriginInteger[hull3];
Print["Shifted vertices: ", shifted3];

(* Full decomposition from original points *)
decomp3 = ComputeDecomposition[points3];
Print["Fan vertices: ", decomp3[[1]]];
Print["Simplices: ", decomp3[[2]]];
Print["Number of sectors: ", Length[decomp3[[2]]]];
Print[];


(* ============================================================================
   Example 4: Feynman-type integrand with exponents

   Integrand:  (x[1] + x[2])^(-1 + d\[Epsilon]) * (x[1]^2 + x[2]^2 + x[1] x[2])^(-1)
   The exponent d\[Epsilon] is symbolic; PolytopeVertices uses Re[exponent].
   ============================================================================ *)

Print["=== Example 4: Symbolic exponents (dimensional regularization) ==="];

integrand4 = (x[1] + x[2])^(-1 + d\[Epsilon]) *
             (x[1]^2 + x[2]^2 + x[1] x[2])^(-1);
vars4 = {x[1], x[2]};

verts4 = PolytopeVertices[integrand4, vars4];
Print["Newton polytope vertices: ", verts4];

decomp4 = ComputeDecomposition[verts4];
Print["Number of sectors: ", Length[decomp4[[2]]]];
Print[];


(* ============================================================================
   Example 5: Positive orthant decomposition

   ComputeDecompositiony intersects the fan with y >= 0, giving
   only the sectors in the positive orthant.
   ============================================================================ *)

Print["=== Example 5: Positive orthant decomposition ==="];

integrand5 = (x[1] + x[2] + x[3])^(-1) * (x[1] + x[3])^(-1) * (x[2] + x[3])^(-1);
vars5 = {x[1], x[2], x[3]};

verts5 = PolytopeVertices[integrand5, vars5];
Print["Newton polytope vertices: ", verts5];

decomp5 = ComputeDecompositiony[verts5];
Print["Positive-orthant fan vertices: ", decomp5[[1]]];
Print["Number of positive-orthant sectors: ", Length[decomp5[[2]]]];

(* Compare with full decomposition *)
decomp5full = ComputeDecomposition[verts5];
Print["Number of sectors (full fan): ", Length[decomp5full[[2]]]];
Print[];


(* ============================================================================
   Example 6: Higher-dimensional example (4 variables)

   Integrand:  (x[1]+x[2]+x[3]+x[4])^(-2) * (x[1] x[3]+x[2] x[4])^(-1)
   ============================================================================ *)

Print["=== Example 6: 4-dimensional example ==="];

integrand6 = (x[1] + x[2] + x[3] + x[4])^(-2) *
             (x[1] x[3] + x[2] x[4])^(-1);
vars6 = {x[1], x[2], x[3], x[4]};

verts6 = PolytopeVertices[integrand6, vars6];
Print["Newton polytope vertices: ", verts6];

decomp6 = ComputeDecomposition[verts6];
Print["Number of dual vertices: ", Length[decomp6[[1]]]];
Print["Number of sectors: ", Length[decomp6[[2]]]];
Print[];


Print["=== All examples complete ==="];
