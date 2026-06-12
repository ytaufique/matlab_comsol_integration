(* convert_1d_to_2d.wl
   Mathematica script: 1D axisymmetric -> 2D Cartesian conversion
   for thin-film spin-coating equations.

   Usage: Get["path/to/convert_1d_to_2d.wl"]
   Or paste each section into a Mathematica notebook and evaluate. *)

sep  = StringRepeat["=", 68];
sep2 = StringRepeat["-", 68];

Print[""];  Print[sep];
Print[" 1D AXISYMMETRIC -> 2D CARTESIAN  |  THIN-FILM SPIN COATING"];
Print[sep];  Print[""];

(* =========================================================================
   PART 1: CORE DIVERGENCE IDENTITY
   ========================================================================= *)
Print["PART 1: CORE IDENTITY"];
Print[sep2];
Print[""];
Print[" Claim:  (1/r) d/dr[ r Ga(r) ]  =  d/dx[ Ga x/r ] + d/dy[ Ga y/r ]"];
Print["         where r = Sqrt[x^2 + y^2]"];
Print[""];
Print[" Proof:"];

(* LHS: cylindrical divergence, substituting r = sqrt(x^2+y^2) *)
lhs = Simplify[
  (1/r * D[r * Ga[r], r]) /. r -> Sqrt[x^2 + y^2],
  {x > 0, y > 0}
];

(* RHS: Cartesian divergence of the radial vector field (Ga*x/r, Ga*y/r) *)
rhs = Simplify[
  D[Ga[Sqrt[x^2 + y^2]] * x / Sqrt[x^2 + y^2], x] +
  D[Ga[Sqrt[x^2 + y^2]] * y / Sqrt[x^2 + y^2], y],
  {x > 0, y > 0}
];

Print["  LHS = ", lhs];
Print["  RHS = ", rhs];
Print["  LHS - RHS = ", FullSimplify[lhs - rhs, {x > 0, y > 0}]];
Print[""];
Print[" => Confirmed: difference = 0. r-factors dissolve into x/r, y/r components."];
Print[""];

(* Corollary: cylindrical Laplacian = Cartesian Laplacian for f(r) *)
Print[" Corollary — Laplacian identity:"];
Print["   (1/r) d/dr[ r df/dr ]  =  d2f/dx2 + d2f/dy2"];
Print[""];

lhsLap = Simplify[
  (1/r * D[r * D[f[r], r], r]) /. r -> Sqrt[x^2 + y^2],
  {x > 0, y > 0}
];
rhsLap = Simplify[
  D[f[Sqrt[x^2 + y^2]], x, x] + D[f[Sqrt[x^2 + y^2]], y, y],
  {x > 0, y > 0}
];
Print["  Cylindrical Laplacian = ", lhsLap];
Print["  Cartesian Laplacian   = ", rhsLap];
Print["  Difference = ", FullSimplify[lhsLap - rhsLap, {x > 0, y > 0}]];
Print[""];
Print[" => Confirmed: Laplacian identity holds."];
Print[""];

(* Chain rule: df/dr -> df/dx, df/dy *)
Print[" Chain rule:  d/dx[ f(r) ] = (df/dr)(x/r)"];
chainCheck = FullSimplify[
  D[f[Sqrt[x^2 + y^2]], x] -
  Derivative[1][f][Sqrt[x^2 + y^2]] * x / Sqrt[x^2 + y^2],
  {x > 0, y > 0}
];
Print["  d/dx[f(r)] - (df/dr)(x/r) = ", chainCheck];
Print[""];
Print[" => Confirmed: df/dx = (df/dr)(x/r),  df/dy = (df/dr)(y/r)"];
Print[""];

(* =========================================================================
   PART 2: SUBSTITUTION RULES TABLE
   ========================================================================= *)
Print["PART 2: SUBSTITUTION RULES"];
Print[sep2];
Print[""];
Print[" ", PaddedForm["1D axisymmetric", 38], "2D Cartesian"];
Print[" ", StringRepeat["-", 38], "  ", StringRepeat["-", 28]];
rules = {
  {"r",                       "Sqrt[x^2 + y^2]"},
  {"df/dr",                   "(x df/dx + y df/dy) / r"},
  {"(1/r) d/dr[r Ga(r)]",    "d/dx[Ga x/r] + d/dy[Ga y/r]"},
  {"Laplacian_cyl(f)",        "d2f/dx2 + d2f/dy2"},
  {"ds/dr",                   "ds/dx=(ds/dr)(x/r),  ds/dy=(ds/dr)(y/r)"},
  {"da[2] = 3r*eta(c)",       "divide eq by eta(c)  => da = 1"},
  {"da[9] = hh",              "divide eq by hh_safe => da = 1"}
};
Scan[Function[row,
  Print[" ", PaddedForm[row[[1]], 38], row[[2]]]
], rules];
Print[""];

(* =========================================================================
   PART 3: EQ 0 — THIN FILM EVOLUTION
   ========================================================================= *)
Print["PART 3: EQUATION CONVERSIONS"];
Print[sep2];
Print[""];
Print[" EQ 0 — Thin film evolution"];
Print[""];
Print[" 1D axisymmetric:"];
Print["   eta(c) dh/dt + (1/r) d/dr[ hh^3/3 (dP/dr - rho w^2 r) ] = -F(c) eta(c)"];
Print[""];
Print[" Radial flux:  Ga0(r) = hh^3/3 * (dP/dr - rho w^2 r)"];
Print[""];
Print[" Apply divergence identity:"];
Print["   Ga0 x/r = hh^3/3 * (dP/dx - rho w^2 x)   [since dP/dr * x/r = dP/dx for radial P]"];
Print["   Ga0 y/r = hh^3/3 * (dP/dy - rho w^2 y)"];
Print[""];
Print[" Divide by eta(c):  da goes from eta(c) -> 1 (constant)"];
Print[""];
Print[" Result — 2D Cartesian:"];
Print["   dh/dt + d/dx[ hh^3/(3 eta) (dP/dx - rho w^2 x) ]"];
Print["         + d/dy[ hh^3/(3 eta) (dP/dy - rho w^2 y) ]  =  -F(c)"];
Print[""];

(* Symbolic check: centrifugal part of Ga0 *)
Print[" Symbolic check on centrifugal term  Ga0_cent(r) = -rho w^2 r^2 * H^3/3:"];
lhsCent = Simplify[
  (1/r * D[r * (-\[Rho]*\[Omega]^2 * r * H^3/3), r]) /. r -> Sqrt[x^2 + y^2],
  {x > 0, y > 0}
];
rhsCent = Simplify[
  D[-\[Rho]*\[Omega]^2*x * H^3/3, x] + D[-\[Rho]*\[Omega]^2*y * H^3/3, y],
  {x > 0, y > 0}
];
Print["  Cylindrical: ", lhsCent];
Print["  Cartesian:   ", rhsCent];
Print["  Equal: ", FullSimplify[lhsCent - rhsCent] === 0];
Print[""];

(* =========================================================================
   PART 3: EQ 1 — CAPILLARY PRESSURE
   ========================================================================= *)
Print[" EQ 1 — Capillary pressure (quasi-static, da = 0)"];
Print[""];
Print[" 1D:   P = gamma * (1/r) d/dr[ r d(h+s)/dr ]  (cylindrical Laplacian)"];
Print[""];
Print[" The Laplacian identity (Part 1 corollary) gives directly:"];
Print[""];
Print[" Result — 2D Cartesian:"];
Print["   P = gamma * ( d2(h+s)/dx2 + d2(h+s)/dy2 )"];
Print[""];
Print[" Substrate gradient (chain rule):"];
Print["   1D: ds/dr = A * flc2hs'(r - r0, eps)"];
Print["   2D: ds/dx = ds/dr * x/r"];
Print["       ds/dy = ds/dr * y/r"];
Print[""];

(* Chain rule verification for a general radial function s(r) *)
chainS = FullSimplify[
  D[s[Sqrt[x^2+y^2]], x] - Derivative[1][s][Sqrt[x^2+y^2]] * x/Sqrt[x^2+y^2],
  {x > 0, y > 0}
];
Print[" Chain rule check: d/dx[s(r)] - (ds/dr)(x/r) = ", chainS];
Print[""];

(* =========================================================================
   PART 3: EQ 2 — SOLVENT CONCENTRATION
   ========================================================================= *)
Print[" EQ 2 — Solvent concentration"];
Print[""];
Print[" 1D:"];
Print["   hh dc/dt + (1/r) d/dr[-D hh r dc/dr]"];
Print["     = c F(c) - (hh^3/(3 eta)) (dP/dr - rho w^2 r) dc/dr"];
Print[""];
Print[" Step 1 — Diffusion flux Ga_c(r) = -D hh dc/dr:"];
Print["   Cartesian: d/dx[-D hh dc/dx] + d/dy[-D hh dc/dy]   (no cross term)"];
Print[""];
Print[" Step 2 — Advection source:"];
Print["   1D: (dP/dr - rho w^2 r) dc/dr          <- scalar product"];
Print["   2D: (dP/dx - rho w^2 x) dc/dx"];
Print["     + (dP/dy - rho w^2 y) dc/dy           <- vector dot product (no d2/dxdy)"];
Print[""];
Print[" Step 3 — da = hh  ->  divide by hh_safe  ->  da = 1"];
Print["   Correction term D*(dhh/dx dc/dx + dhh/dy dc/dy)/hh_safe moved to source"];
Print[""];
Print[" Result — 2D Cartesian (physical form):"];
Print["   hh dc/dt - D [ d/dx(hh dc/dx) + d/dy(hh dc/dy) ]"];
Print["     = c F(c) - (hh^3/(3 eta))[(dP/dx - rho w^2 x) dc/dx"];
Print["                               + (dP/dy - rho w^2 y) dc/dy]"];
Print[""];

(* =========================================================================
   PART 4: VERIFY NO CROSS DERIVATIVES
   ========================================================================= *)
Print["PART 4: NO CROSS DERIVATIVES — SYMBOLIC VERIFICATION"];
Print[sep2];
Print[""];

(* Thin film: div of (h^3/3*(dP/dx - rho w^2 x), h^3/3*(dP/dy - rho w^2 y)) *)
divHH = Expand[
  D[h[x,y]^3/3 * (D[P[x,y], x] - \[Rho]*\[Omega]^2*x), x] +
  D[h[x,y]^3/3 * (D[P[x,y], y] - \[Rho]*\[Omega]^2*y), y]
];

hasMixedP = ! FreeQ[divHH, Derivative[1,1][P]];
hasMixedH = ! FreeQ[divHH, Derivative[1,1][h]];
Print[" Thin film divergence — mixed d2P/dxdy present? ", hasMixedP];
Print[" Thin film divergence — mixed d2h/dxdy present? ", hasMixedH];
Print[""];

(* Concentration diffusion: div of (D*h*dc/dx, D*h*dc/dy) *)
divC = Expand[
  D[D*h[x,y]*D[c[x,y], x], x] + D[D*h[x,y]*D[c[x,y], y], y]
];
hasMixedC  = ! FreeQ[divC, Derivative[1,1][c]];
hasMixedHC = ! FreeQ[divC, Derivative[1,1][h]];
Print[" Concentration diffusion — mixed d2c/dxdy present? ", hasMixedC];
Print[" Concentration diffusion — mixed d2h/dxdy present? ", hasMixedHC];
Print[""];

(* Advection dot product *)
advTerm = (D[P[x,y], x] - \[Rho]*\[Omega]^2*x)*D[c[x,y], x] +
          (D[P[x,y], y] - \[Rho]*\[Omega]^2*y)*D[c[x,y], y];
hasMixedAdv = ! FreeQ[advTerm, Derivative[1,1][c]] ||
              ! FreeQ[advTerm, Derivative[1,1][P]];
Print[" Advection source — any mixed derivative present? ", hasMixedAdv];
Print[""];
Print[" => All False: no cross derivatives d2/dxdy appear anywhere."];
Print[" => The 1D->2D extension adds only x- and y-direction terms independently."];
Print[""];
Print[sep];
