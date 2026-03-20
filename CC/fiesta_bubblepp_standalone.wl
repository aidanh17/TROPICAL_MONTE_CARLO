(* ============================================================================
   FIESTA-only evaluation of bubblepp presector 1 (simplex form)
   Run in a clean session to avoid symbol conflicts with tropical MC.

   The [0,inf) integral has been GL(1)-converted to simplex form:
   6 original variables -> 7 simplex variables with delta(1-sum)

   Kinematic point: k1=k2=k3=1 (equilateral)
   ============================================================================ *)

SetDirectory["/usr/local/fiesta/FIESTA5"];
Get["FIESTA5.m"];
SetOptions[FIESTA, "NumberOfSubkernels" -> 4, "NumberOfLinks" -> 8];

Print["FIESTA5 loaded."];
Print[];

(* The homogenized integral on the simplex:
   Variables: x[1]=x0 (new), x[2]..x[6] = original x[1]..x[5]
   Delta: delta(1 - x[1] - x[2] - ... - x[6])

   Functions and their degrees (using ep = FIESTA's regulator):
   f1  = x[1]              degree = 6 + 4*ep           (GL(1) Jacobian)
   f2  = x[2]              degree = -1/4 - 6*ep        (monomial x1)
   f3  = x[3]              degree = -1/4 - 6*ep        (monomial x2)
   f4  = x[4]              degree = -1/4 - 6*ep        (monomial x3)
   f5  = x[5]              degree = -1/4 - 6*ep        (monomial x4)
   f6  = x[6]              degree = 3/2 - 12*ep        (monomial x5)
   f7  = x[1]+x[6]         degree = 7/2 - 12*ep        (poly: 1+x5)
   f8  = x[2]+x[3]         degree = 1/4 + 5*ep         (poly: x1+x2)
   f9  = x[4]+x[5]         degree = 1/4 + 5*ep         (poly: x3+x4)
   f10 = x[2]+x[3]+x[4]+x[5] degree = -3/2 + ep       (poly: x1+x2+x3+x4)
   f11 = P5_hom             degree = -5 + 11*ep         (big kinematic poly)
   f12 = 4(x[2]+x[3]+x[4]+x[5]) degree = 5 - 11*ep    (poly: 4*sum)
*)

(* Build the homogenized P5 at k1=k2=k3=1 *)
P5hom = x[1]^4 +
  8 x[1]^3 x[2] + 4 x[1]^2 x[2]^2 +
  8 x[1]^3 x[3] + 8 x[1]^2 x[2] x[3] + 4 x[1]^2 x[3]^2 +
  8 x[1]^3 x[4] + 8 x[1]^2 x[2] x[4] + 8 x[1]^2 x[3] x[4] + 4 x[1]^2 x[4]^2 +
  8 x[1]^3 x[5] + 8 x[1]^2 x[2] x[5] + 8 x[1]^2 x[3] x[5] + 8 x[1]^2 x[4] x[5] + 4 x[1]^2 x[5]^2 +
  12 x[1]^2 x[2] x[6] + 12 x[1]^2 x[3] x[6] + 16 x[1] x[2] x[3] x[6] + 16 x[1] x[3]^2 x[6] +
  12 x[1]^2 x[4] x[6] + 16 x[1] x[3] x[4] x[6] +
  12 x[1]^2 x[5] x[6] + 16 x[1] x[2] x[5] x[6] + 32 x[1] x[3] x[5] x[6] + 16 x[1] x[4] x[5] x[6] + 16 x[1] x[5]^2 x[6] +
  16 x[2] x[3] x[6]^2 + 16 x[3]^2 x[6]^2 + 16 x[3] x[4] x[6]^2 +
  16 x[2] x[5] x[6]^2 + 32 x[3] x[5] x[6]^2 + 16 x[4] x[5] x[6]^2 + 16 x[5]^2 x[6]^2;

Print["P5_hom has ", Length[If[Head[Expand[P5hom]] === Plus, List @@ Expand[P5hom], {P5hom}]], " terms"];
Print[];

(* Set up FIESTA call *)
functions = {
  x[1],                         (* x0 *)
  x[2], x[3], x[4], x[5],     (* x1..x4 *)
  x[6],                         (* x5 *)
  x[1] + x[6],                  (* 1+x5 homogenized *)
  x[2] + x[3],                  (* x1+x2 *)
  x[4] + x[5],                  (* x3+x4 *)
  x[2] + x[3] + x[4] + x[5],  (* x1+x2+x3+x4 *)
  P5hom,                         (* big kinematic poly *)
  4 (x[2] + x[3] + x[4] + x[5])  (* 4(x1+x2+x3+x4) *)
};

degrees = {
  6 + 4 ep,          (* x0 power *)
  -1/4 - 6 ep,       (* x1 mono *)
  -1/4 - 6 ep,       (* x2 mono *)
  -1/4 - 6 ep,       (* x3 mono *)
  -1/4 - 6 ep,       (* x4 mono *)
  3/2 - 12 ep,       (* x5 mono *)
  7/2 - 12 ep,       (* (1+x5) poly *)
  1/4 + 5 ep,        (* (x1+x2) poly *)
  1/4 + 5 ep,        (* (x3+x4) poly *)
  -3/2 + ep,         (* (x1+x2+x3+x4) poly *)
  -5 + 11 ep,        (* P5 poly *)
  5 - 11 ep          (* 4*sum poly *)
};

Print["Degrees at ep=0: ", degrees /. ep -> 0];
Print[];

Print["Running FIESTA SDEvaluateDirect on 6D simplex integral..."];
Print["  12 functions, delta(1-x[1]-...-x[6])"];
Print["  maxeval = 5000000"];
Print[];

result = SDEvaluateDirect[
  functions, degrees, 0,
  {{1, 2, 3, 4, 5, 6}},
  IntegratorOptions -> {{"maxeval", "5000000"}}
];

Print[];
Print["================================================================"];
Print["  FIESTA Result (presector 1, k1=k2=k3=1)"];
Print["================================================================"];
Print["  ", result];
Print[];

(* Compare with tropical MC *)
Print["================================================================"];
Print["  Comparison with Tropical MC"];
Print["================================================================"];
Print[];
Print["  Tropical MC raw values (2M samples, presector 1):"];
Print["    main (TOTAL) = 0.025485"];
Print["    G0           = 0.000513"];
Print[];
Print["  Tropical MC gives the integral WITHOUT the prefactor."];
Print["  FIESTA also gives the integral WITHOUT the prefactor."];
Print["  They should match directly (both are the bare presector integral)."];
Print[];
Print["  Expected FIESTA:"];
Print["    1/ep coeff ~ G0 = 0.000513  (from tropical MC)"];
Print["    ep^0 coeff ~ TOTAL integral = 0.0255  (from tropical MC)"];
Print[];
Print["  (The actual pole/finite decomposition may differ in detail"];
Print["   but the Laurent expansion coefficients should agree.)"];
