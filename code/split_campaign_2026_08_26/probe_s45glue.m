// glue E (X1(11) point s0=4/5 over K=Q(sqrt 11)) with its conjugate along
// the unique equivariant 2-torsion iso, via the HLP/BHLS formula of
// Genus2Elliptic2 reimplemented over K; then check Igusa rationality.
SetColumns(0);
SetMemoryLimit(12*10^9);
Q := Rationals();
K<w> := QuadraticField(11);
RK<xk> := PolynomialRing(K);
s0 := K!(4/5);
r0 := Roots(xk^2 - (s0^3-3*s0^2+4*s0)*xk + s0)[1][1];
c0 := s0*(r0-1); b0 := r0*c0;
E := EllipticCurve([1-c0, -b0, -b0, 0, 0]);
sig := hom< K -> K | -w >;
Es := EllipticCurve([sig(a) : a in aInvariants(E)]);
assert Order(E![0,0]) eq 11;

f := HyperellipticPolynomials(WeierstrassModel(E));
g := HyperellipticPolynomials(WeierstrassModel(Es));
printf "f = %o\ng = %o\n", f, g;

// E[2] ~ E^sigma[2] over K  ==>  f and g have the SAME splitting field over K
L := SplittingField(f);
printf "splitting field of f over K: rel degree %o\n", Degree(L);
RL := PolynomialRing(L);
alpha := [ r[1] : r in Roots(RL!f) ];
gamma := [ r[1] : r in Roots(RL!g) ];
printf "roots found: f %o, g %o\n", #alpha, #gamma;
assert #alpha eq 3 and #gamma eq 3;
Df := Discriminant(f); Dg := Discriminant(g);

// Galois group of L over K (glue is over K)
GA := Automorphisms(L);
GK := [ gA : gA in GA | gA(L!w) eq L!w ];
printf "automorphisms: %o total, %o over K\n", #GA, #GK;
Gaction := func< rs | [[ Index(rs, gA(r)) : r in rs ] : gA in GK ] >;

a := [[ alpha[i]-alpha[j] : j in [1..3]] : i in [1..3]];
Rt<t> := PolynomialRing(L);
nglue := 0;
for sgm in SymmetricGroup(3) do
    beta := [ gamma[i^sgm] : i in [1..3] ];
    b := [[ beta[i]-beta[j] : j in [1..3]] : i in [1..3]];
    if alpha[1]*b[3][2] + alpha[2]*b[1][3] + alpha[3]*b[2][1] ne 0 and Gaction(alpha) eq Gaction(beta) then
        nglue +:= 1;
        a1 := a[3][2]^2/b[3][2]+a[2][1]^2/b[2][1]+a[1][3]^2/b[1][3];
        a2 := alpha[1]*b[3][2]+alpha[2]*b[1][3]+alpha[3]*b[2][1];
        b1 := b[3][2]^2/a[3][2]+b[2][1]^2/a[2][1]+b[1][3]^2/a[1][3];
        b2 := beta[1]*a[3][2]+beta[2]*a[1][3]+beta[3]*a[2][1];
        A := Dg*a1/a2;  B := Df*b1/b2;
        sext := -(A*a[2][1]*a[1][3]*t^2+B*b[2][1]*b[1][3])
                *(A*a[3][2]*a[2][1]*t^2+B*b[3][2]*b[2][1])
                *(A*a[1][3]*a[3][2]*t^2+B*b[1][3]*b[3][2]);
        okc := true; sextK := 0;
        try sextK := RK![ K!c : c in Coefficients(sext) ]; catch e; okc := false; end try;
        if not okc then printf "matching %o: sextic NOT over K\n", sgm; continue; end if;
        printf "matching %o: sextic over K found\n", sgm;
        C := HyperellipticCurve(sextK);
        IC := IgusaClebschInvariants(C);
        inQ := &and[ ic in Q : ic in IC ];
        printf "   IgusaClebsch rational: %o\n", inQ;
        if inQ then
            ICQ := [ Q!ic : ic in IC ];
            printf "   IC = %o\n", ICQ;
            // Mestre reconstruction over Q
            C0 := 0; okm := true;
            try C0 := HyperellipticCurveFromIgusaClebsch(ICQ); catch e; okm := false; printf "   Mestre failed: %o\n", e`Object; end try;
            if okm then
                C0 := ReducedMinimalWeierstrassModel(C0);
                printf "   Q-model: %o\n", HyperellipticPolynomials(C0);
            end if;
        end if;
    end if;
end for;
printf "GLUE_DONE matchings=%o\n", nglue;
quit;
