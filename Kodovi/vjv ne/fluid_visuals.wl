(* ========================================================= *)
(*  TEMPLATE – Konsistentan monokromatski stil                *)
(* ========================================================= *)

ClearAll["Global`*"];

(* Osnovni stil *)
$plotStyle = {
   Frame -> True, Axes -> False, Background -> White,
   FrameStyle -> Directive[Black, 12],
   ImageSize -> 800,
   PlotRange -> All
};

$monoContours = Directive[Black, Thin];
$monoField = Directive[Black, Opacity[0.7]];

(* Kompleksni potencijal: F = W'(z) = u - i v *)
complexVelocity[W_, z_] := D[W, z];

(* Iz kompleksne brzine izvući realnu brzinu (u,v) *)
uvFromF[F_, x_, y_] := Module[{z = x + I y, f = F /. z -> x + I y},
  {Re[f], -Im[f]}
];

(* Strujnice (psi) i potencijalne linije (phi) *)
phiFromW[W_, x_, y_] := Re[W /. (z -> x + I y)];
psiFromW[W_, x_, y_] := Im[W /. (z -> x + I y)];

(* Streamlines plot *)
streamPlot[W_, {x_, xmin_, xmax_}, {y_, ymin_, ymax_}, opts___] := Module[
  {F = complexVelocity[W, z]},
  StreamPlot[
    Evaluate[uvFromF[F, x, y]],
    {x, xmin, xmax}, {y, ymin, ymax},
    StreamStyle -> $monoField,
    StreamPoints -> Fine,
    StreamScale -> None,
    Evaluate[$plotStyle],
    opts
  ]
];

(* Contour plot of streamfunction and potential *)
contourPlot[W_, {x_, xmin_, xmax_}, {y_, ymin_, ymax_}, n_:20, opts___] :=
  ContourPlot[
    {phiFromW[W, x, y], psiFromW[W, x, y]},
    {x, xmin, xmax}, {y, ymin, ymax},
    Contours -> n,
    ContourStyle -> {$monoContours, Directive[GrayLevel[0.5], Thin]},
    Evaluate[$plotStyle],
    opts
  ];

(* Konformno preslikavanje: mreža *)
conformalMapGrid[g_, {x_, xmin_, xmax_}, {y_, ymin_, ymax_},
   nx_:10, ny_:10, opts___] := Module[
  {xlines, ylines},
  xlines = Table[
    ParametricPlot[
      Evaluate[g[x0 + I t]],
      {t, ymin, ymax},
      PlotStyle -> $monoContours],
    {x0, Subdivide[xmin, xmax, nx]}
  ];
  ylines = Table[
    ParametricPlot[
      Evaluate[g[t + I y0]],
      {t, xmin, xmax},
      PlotStyle -> $monoContours],
    {y0, Subdivide[ymin, ymax, ny]}
  ];
  Show[Join[xlines, ylines], Evaluate[$plotStyle], opts]
];

(* ========================================================= *)
(*  POTENCIJALNI TOKOVI – primjeri iz potencijalni_tokovi.txt *)
(* ========================================================= *)

(* 2.6.1 Uniformni tok *)
W261[U_, V_] := (U - I V) z;

ex261[] := streamPlot[W261[1, 0], {x, -2, 2}, {y, -2, 2}];

(* 2.6.2 Optjecanje cilindra bez cirkulacije *)
W262[U_, r_] := U (z + r^2/z);

ex262[] := streamPlot[W262[1, 1], {x, -3, 3}, {y, -3, 3},
   RegionFunction -> Function[{x, y}, x^2 + y^2 >= 1.01]
];

(* 2.6.3 Točkasti vrtlog *)
W263[Γ_] := (Γ/(2 Pi I)) Log[z];

ex263[] := streamPlot[W263[1], {x, -2, 2}, {y, -2, 2},
   RegionFunction -> Function[{x, y}, x^2 + y^2 >= 0.1]
];

(* 2.6.5 Cilindar s cirkulacijom *)
W265[U_, r_, Γ_] := U (z + r^2/z) - (I Γ/(2 Pi)) Log[z];

ex265[] := streamPlot[W265[1, 1, 2], {x, -3, 3}, {y, -3, 3},
   RegionFunction -> Function[{x, y}, x^2 + y^2 >= 1.01]
];

(* 2.6.7 Dipol *)
W267[K_] := -K/(2 Pi z);

ex267[] := streamPlot[W267[1], {x, -2, 2}, {y, -2, 2},
   RegionFunction -> Function[{x, y}, x^2 + y^2 >= 0.1]
];

(* 2.6.8 Izvor/Ponor *)
W268[D_] := (D/(2 Pi)) Log[z];

ex268[] := streamPlot[W268[1], {x, -2, 2}, {y, -2, 2},
   RegionFunction -> Function[{x, y}, x^2 + y^2 >= 0.1]
];

(* 2.6.9 Tok u kutu *)
W269[a_] := a z^2;

ex269[] := streamPlot[W269[1], {x, 0, 2}, {y, 0, 2}];

(* 2.6.10 Aeroprofil (Kutta + Joukowski ideja) *)
(* Approximacija: Joukowski map + uniform flow + circulation *)
joukowski[z_] := z + 1/z;

W2610[U_, Γ_] := U z - (I Γ/(2 Pi)) Log[z];

ex2610[] := Module[{W = W2610[1, 1]},
  streamPlot[W /. z -> inverseJoukowski[w], {x, -3, 3}, {y, -3, 3},
   RegionFunction -> Function[{x, y}, x^2 + y^2 >= 1.01]]
];

(* ========================================================= *)
(*  KONFORMNA PRESLIKAVANJA – primjeri iz konformno.txt       *)
(* ========================================================= *)

(* 5.1 Translations *)
g51[β_] := (# + β) &;
ex51[] := conformalMapGrid[g51[0.5 + 0.3 I], {x, -1, 1}, {y, -1, 1}];

(* 5.2 Scaling + rotation *)
g52[ρ_, φ_] := (ρ Exp[I φ] #) &;
ex52[] := conformalMapGrid[g52[1.5, Pi/6], {x, -1, 1}, {y, -1, 1}];

(* 5.3 Affine *)
g53[α_, β_] := (α # + β) &;
ex53[] := conformalMapGrid[g53[1 + I, -0.2 + 0.1 I], {x, -1, 1}, {y, -1, 1}];

(* 5.4 Inversion *)
g54 := (1/#) &;
ex54[] := conformalMapGrid[g54, {x, -2, 2}, {y, -2, 2},
   RegionFunction -> Function[{x, y}, x^2 + y^2 >= 0.2]
];

(* 5.5 Exponential strip->wedge *)
g55 := Exp[#] &;
ex55[] := conformalMapGrid[g55, {x, -2, 2}, {y, -Pi/2, Pi/2}];

(* 5.6 Squaring map *)
g56 := (#^2) &;
ex56[] := conformalMapGrid[g56, {x, 0, 2}, {y, 0, 2}];

(* 5.7 Linear fractional map to disk *)
g57[z_] := (z - 1)/(z + 1);
ex57[] := conformalMapGrid[g57, {x, 0.1, 3}, {y, -2, 2}];

(* 5.8 Disk automorphism *)
g58[α_] := Function[z, (z - α)/(α z - 1)];
ex58[] := conformalMapGrid[g58[0.3], {x, -1, 1}, {y, -1, 1}];

(* 5.13 Joukowski map *)
g513 := (1/2) (# + 1/#) &;
ex513[] := conformalMapGrid[g513, {x, -2, 2}, {y, -2, 2},
   RegionFunction -> Function[{x, y}, x^2 + y^2 >= 0.2]
];

(* 5.18 Half-disk to disk *)
g518 := ( #^2 + 2 I # + 1 )/( #^2 - 2 I # + 1 ) &;
ex518[] := conformalMapGrid[g518, {x, -1, 1}, {y, 0, 1}];

(* 5.19 Nonconcentric annulus to concentric *)
g519[c_] := Function[z, (2 c z - 1)/(z - 2 c)];
ex519[] := conformalMapGrid[g519[0.2], {x, -1, 1}, {y, -1, 1}];

(* ========================================================= *)
(*  PRIMJENE – primjeri iz 6.x                                *)
(* ========================================================= *)

(* 6.3 Right half-plane -> disk; prikaz mapiranja *)
ex63[] := conformalMapGrid[g57, {x, 0.1, 3}, {y, -2, 2}];

(* 6.4 Non-coaxial cable: annulus -> concentric annulus *)
ex64[] := ex519[];

(* 6.5 Tok oko vanjskog kuta: W = (i z)^(2/3) *)
W65 := (I z)^(2/3);
ex65[] := streamPlot[W65, {x, -2, 2}, {y, -2, 2}];

(* 6.6 Horizontal plate: W = z *)
W66 := z;
ex66[] := streamPlot[W66, {x, -2, 2}, {y, -2, 2}];

(* 6.7 Disk: W = 1/2 (z + 1/z) *)
W67 := (1/2) (z + 1/z);
ex67[] := streamPlot[W67, {x, -3, 3}, {y, -3, 3},
   RegionFunction -> Function[{x, y}, x^2 + y^2 >= 1.01]
];

(* 6.8 Tilted plate: W = e^{iφ}(z cosφ - i sinφ sqrt(z^2 - e^{-2iφ})) *)
W68[φ_] := Exp[I φ] (z Cos[φ] - I Sin[φ] Sqrt[z^2 - Exp[-2 I φ]]);
ex68[] := streamPlot[W68[Pi/8], {x, -3, 3}, {y, -3, 3}];

(* 6.9 Airfoil via Joukowski *)
W69[α_, β_] := Module[{w = α z + β},
  (1/2) (w + 1/w)
];
ex69[] := conformalMapGrid[W69[1, -0.1 + 0.2 I], {x, -2, 2}, {y, -2, 2}];

(* 6.12 Green function for half-plane: use g57 and contour *)
G612[z_, ζ_] := (1/(2 Pi)) Log[( (ζ - 1)/(ζ + 1) (z - 1)/(z + 1) - 1 )/
   ( (z - 1)/(z + 1) - (ζ - 1)/(ζ + 1) )];

(* ========================================================= *)
(*  POZIVI (primjeri koje možeš renderirati)                  *)
(* ========================================================= *)

(* ex261[]; ex262[]; ex263[]; ex265[]; ex267[]; ex268[]; ex269[]; ex2610[]; *)
(* ex51[]; ex52[]; ex53[]; ex54[]; ex55[]; ex56[]; ex57[]; ex58[]; ex513[]; ex518[]; ex519[]; *)
(* ex63[]; ex64[]; ex65[]; ex66[]; ex67[]; ex68[]; ex69[]; *)