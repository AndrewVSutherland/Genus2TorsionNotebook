// Leprevost-style second engine for [2,22], step 0: the two-contact family.
//
// Framework (rebuilt from first principles; the 1995 scan's display math is
// image-only): C: y^2 = f, f monic sextic, with TWO contact representations
//   f = h1^2 + c1 * x^4 (x-1)      (exponents (4,1), m1 = 5)
//   f = h2^2 + c2 * x   (x-1)^3    (exponents (1,3), m2 = 4)
// h1, h2 monic cubics.  Relations alpha_i u' + beta_i v' = (m_i-3) w in J
// with minor 4*3 - 1*1 = 11 => suitable solution components carry
// 11-torsion of [P]-[oo] type (distinct from the BLP engine's D_inf class).
//
// System: h1^2 + c1 x^4(x-1) = h2^2 + c2 x(x-1)^3 -- 6 coefficient
// equations in 8 unknowns: a 2-parameter family.  This script SAMPLES the
// family: freeze (a, b) = random rationals, solve the 0-dimensional system
// in (e,p,q,r,c1,c2) by Variety over Q, and for every rational solution
// build the curve and compute exact TorsionSubgroup.  Components whose
// samples show 11 | #tors identify the live branch and its decoration.
//
// Run: magma -b code/lep22_engine_setup.m > results/lep22_engine_setup.log

SetColumns(0);
SetSeed(1);
SetMemoryLimit(4*10^9);

Q := Rationals();
PQ<X> := PolynomialRing(Q);

samples := [
    <1/2, 1/3>, <2, -1>, <1/5, 3>, <-1/2, 2/3>, <3, 1/2>, <-2/3, -1/4>
];

for smp in samples do
    a0 := Q!smp[1]; b0 := Q!smp[2];
    // freeze (a,b); unknowns (e,p,q,r,c1,c2) + Rabinowitsch w for c1*c2 != 0
    R6<e,p,q,r,c1,c2,w> := PolynomialRing(Q, 7);
    Px<x> := PolynomialRing(R6);
    h1 := x^3 + a0*x^2 + b0*x + e;
    h2 := x^3 + p*x^2 + q*x + r;
    dif := (h1^2 + c1*x^4*(x-1)) - (h2^2 + c2*x*(x-1)^3);
    eqs := [ Coefficient(dif, i) : i in [0..5] ] cat [ w*c1*c2 - 1 ];
    I := ideal< R6 | eqs >;
    if Dimension(I) ne 0 then
        printf "SAMPLE a=%o b=%o dim=%o (skip)\n", a0, b0, Dimension(I);
        continue;
    end if;
    V := Variety(I);
    printf "SAMPLE a=%o b=%o rational_solutions=%o\n", a0, b0, #V;
    for pt in V do
        e0 := pt[1]; p0 := pt[2]; q0 := pt[3]; r0 := pt[4]; cc1 := pt[5]; cc2 := pt[6];
        if cc1 eq 0 or cc2 eq 0 then
            printf "  SOL degenerate (c=0) e=%o p=%o\n", e0, p0;
            continue;
        end if;
        H1 := X^3 + a0*X^2 + b0*X + e0;
        f := H1^2 + cc1*X^4*(X-1);
        if Degree(f) ne 6 or Discriminant(f) eq 0 then
            printf "  SOL singular e=%o c1=%o\n", e0, cc1;
            continue;
        end if;
        C := HyperellipticCurve(f);
        try
            Cm := ReducedMinimalWeierstrassModel(C);
            C := SimplifiedModel(Cm);
        catch err ; end try;
        T := TorsionSubgroup(Jacobian(C));
        invs := Invariants(T);
        printf "  SOL e=%o p=%o q=%o r=%o c1=%o c2=%o TORSION=%o div11=%o\n",
               e0, p0, q0, r0, cc1, cc2, invs, #T mod 11 eq 0;
    end for;
end for;
printf "LEP22_SETUP_DONE\n";
quit;
