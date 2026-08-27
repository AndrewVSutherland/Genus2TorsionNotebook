// deep dive: the X1(11) quadratic point s0 = 4/5 over K = Q(sqrt 11) with
// E[2] ~ E^sigma[2]: is it isogeny-induced (degenerate glue) or genuine?
// If genuine: glue over K, check Igusa invariants rational, reconstruct /Q.
SetColumns(0);
SetMemoryLimit(12*10^9);
AttachSpec("/home/claude/Magma/magma.spec");
Attach("/home/claude/torsion_jac/genus2.m");
Q := Rationals();

K<w> := QuadraticField(11);
RK<xk> := PolynomialRing(K);
s0 := K!(4/5);
rr := Roots(xk^2 - (s0^3-3*s0^2+4*s0)*xk + s0);
printf "r-fiber roots: %o\n", #rr;
r0 := rr[1][1];
c0 := s0*(r0-1); b0 := r0*c0;
E := EllipticCurve([1-c0, -b0, -b0, 0, 0]);
printf "order of (0,0): %o\n", Order(E![0,0]);
printf "torsion of E(K): %o\n", Invariants(TorsionSubgroup(E));
sig := hom< K -> K | -w >;
Es := EllipticCurve([sig(a) : a in aInvariants(E)]);
jv := jInvariant(E); js := jInvariant(Es);
printf "j = %o\n", jv;
printf "j sigma-conjugate distinct: %o\n", jv ne js;

// isogeny between E and E^sigma? degrees 2,3,5,7,11,13
for d in [2,3,5,7,11,13] do
    ok := true; val := K!1;
    try val := Evaluate(ClassicalModularPolynomial(d), [jv, js]); catch e; ok := false; end try;
    if ok then printf "Phi_%o(j, js) = 0: %o\n", d, val eq 0; end if;
end for;

// 2-division cubic over K: Galois structure
a1 := 1-c0; b2 := a1^2-4*b0; b4 := -a1*b0; b6 := b0^2;
cub := xk^3 + b2*xk^2 + 8*b4*xk + 16*b6;
printf "cubic irreducible /K: %o, disc square in K: %o\n",
    IsIrreducible(cub), IsSquare(Discriminant(cub));

// attempt the gluing over K
glues := [];
try
    glues := Genus2Elliptic2(E, Es);
catch e;
    printf "Genus2Elliptic2 error: %o\n", e`Object;
end try;
printf "Genus2Elliptic2 returned %o curve(s) over K\n", #glues;
for i in [1..#glues] do
    C := glues[i];
    J := IgusaInvariants(C);
    inQ := &and[ ji in Q : ji in J ];
    printf "glue %o: Igusa invariants rational: %o\n", i, inQ;
    if inQ then
        printf "   J2..J10 = %o\n", [Q!ji : ji in J];
    end if;
end for;
quit;
