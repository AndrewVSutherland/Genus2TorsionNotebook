// claude_ov_erratum_ss2.m -- genus2red ERRATUM propagation (2026-07-25).
//
// WHY.  After the erratum, the paper's [2,2,14] sentence ("the displayed curve
// --- of conductor 2^2*9550095752925 --- the SMALLEST of ten generic
// witnesses") has exactly one soft spot:
//
//   curve 1 (displayed, (-10,-10/7,-1/2)):  N = 2^2 * 9550095752925 = 38200383011700
//                                           v_2(disc_min) = 10 < 12, Ogg unwarned
//   curve 2 ((-5,-15/8,-15/22)):            odd part 17934836205825,
//                                           v_2(disc_min) = 14, Magma's Ogg gives 2^2
//
// 2^2 * 17934836205825 = 71739344823300 > 38200383011700, fine -- but if curve
// 2's true 2-part were 2^1 we would get 35869672411650 < 38200383011700 and the
// "smallest" claim would flip.  Every other witness is orders of magnitude
// larger already in its ODD part, so this is the only comparison at risk.
//
// THE TEST.  f_2 = 2u + t + delta (u = unipotent rank, t = toric rank of the
// special fibre of the Neron model, delta = Swan).  f_2 = 1 forces
// (u,t,delta) = (0,1,0), i.e. J SEMISTABLE at 2, i.e. the special fibre of the
// minimal regular model of C over Z_2 is REDUCED.  So compute the minimal
// regular model at 2 and read the multiplicities: any multiplicity > 1 proves
// f_2 >= 2 with no appeal to Ogg's formula at all.  If instead the fibre is
// reduced then f_2 = t = b_1(dual graph) = e - v + 1, also computable here.
//
// Magma API note: the `multiplicity` field of M`abstract_fibres is NOT
// assigned until something forces it; calling IntersectionMatrix(M) first does.
// (IntersectionMatrix returns ONE value, not two.)
//
// Usage:   magma -b IDX:=2 MemGB:=24 code/claude_ov_erratum_ss2.m
//          IDX in [1..11] = the contact-7 three-root points; IDX := 0 = the
//          BLP C4-corrected [2,22] curve; IDX := -1 = control 12300.e.2.
// Markers: SS2 / ERRATUM_SS2_DONE
SetColumns(0);
if not assigned IDX   then IDX := 2;    elif Type(IDX)   eq MonStgElt then IDX   := StringToInteger(IDX);   end if;
if not assigned MemGB then MemGB := 24; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q); Z := Integers();

Gfun := func< v | -(v^5 - v^3 - v^2/2)/(v+1)^2 >;
trips := [
  [-10, -10/7, -1/2], [-5, -15/8, -15/22], [-3, -3/4, -3/5],
  [-15/8, -15/19, -1/2], [-5/18, -10/49, 4/17], [-4/9, -4/25, 4/17],
  [-511/61, -511/625, -1/2], [-165/41, -33/16, -165/289],
  [-164/297, -1/2, 164/361], [-17/50, -34/189, 34/121], [-1/2, -13/49, 13/50]
];

if IDX ge 1 then
    T := trips[IDX]; s := T[1]; t := T[2]; u := T[3];
    c4 := (Gfun(s)-Gfun(t))/(s^2-t^2); c0 := Gfun(s) - c4*s^2;
    b := c4 - 2; a := 9/2 - c0 - c4;
    h := 1 - (7/2)*x + a*x^2 + b*x^3;
    f := (h^2 + (x-1)^7) div x^2;
    for w in [s,t,u] do assert Evaluate(f, 1-w^2) eq 0; end for;
    den := LCM([Denominator(co) : co in Coefficients(f)]);
    C := HyperellipticCurve(P![ co*den^2 : co in Coefficients(f) ]);
    tag := Sprintf("contact7 idx=%o (s,t,u)=(%o, %o, %o)", IDX, s, t, u);
elif IDX eq 0 then
    C := HyperellipticCurve(x^6-18*x^5-4001*x^4-22524*x^3+859039*x^2-1926258*x-9043839);
    tag := "blpC4 [2,22]";
else
    C := HyperellipticCurve(x^6-7*x^5+8*x^4+16*x^3-5*x^2-5*x, x^2+x);
    tag := "control 12300.e.2 (published N = 12300, v_2(N) = 2)";
end if;

Cm := ReducedMinimalWeierstrassModel(C);
Dm := Z!Discriminant(Cm);
v2D := Valuation(AbsoluteValue(Dm), 2);
printf "SS2 %o\n", tag;
printf "SS2 minimal model = %o\n", Cm;
printf "SS2 v_2(disc_min) = %o  (Magma warns when >= 12: %o)\n", v2D, v2D ge 12;

t0 := Cputime();
M := RegularModel(Cm, 2);
IM := IntersectionMatrix(M);            // forces the multiplicities to be assigned
printf "SS2 RegularModel + IntersectionMatrix in %o s\n", Cputime(t0);
af := M`abstract_fibres;
mults := [ (assigned r`multiplicity) select r`multiplicity else -1 : r in af ];
m := #mults;
printf "SS2 components of the special fibre at 2 : m = %o\n", m;
printf "SS2 multiplicities = %o\n", mults;
red := &and[ mu eq 1 : mu in mults ];
printf "SS2 special fibre REDUCED = %o\n", red;
printf "SS2 => J semistable at 2 = %o ; f_2 = 1 possible = %o ; hence f_2 >= 2 proved = %o\n",
    red, red, not red;
printf "SS2 Ogg value f_2 = v_2(disc_min) + 1 - m = %o\n", v2D + 1 - m;
// Toric rank of a SEMISTABLE (= reduced) special fibre = b_1 of the dual graph
// = e - v + 1, computed straight off the intersection matrix.  This is
// independent of Ogg's formula: semistable => u = 0 and Swan delta = 0, so
// f_2 = t = b_1.
if red then
    nv := m; nedge := 0;
    for i in [1..m] do for j in [i+1..m] do nedge +:= Z!IM[i][j]; end for; end for;
    b1 := nedge - nv + 1;
    printf "SS2 dual graph: v = %o, e = %o, b_1 = e - v + 1 = %o\n", nv, nedge, b1;
    printf "SS2 f_2 (toric rank, NO Ogg) = %o ; agrees with Ogg = %o\n",
        b1, b1 eq (v2D+1-m);
end if;
printf "SS2 intersection matrix = %o\n", IM;
printf "SS2 component group = %o\n", Invariants(ComponentGroup(M));
print "ERRATUM_SS2_DONE";
quit;
