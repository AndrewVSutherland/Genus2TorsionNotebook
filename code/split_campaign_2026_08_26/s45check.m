// two checks:
// (A) my over-K reimplementation vs Genus2Elliptic2 on a known Q-pair with
//     the same 2-torsion field (validation of the transplant);
// (B) the K-glue CK for the s0=4/5 point: is #Jac(CK)(k_p) equal to
//     #E(k_p)*#E^sigma(k_p) at split primes of K? (is CK the right glue?)
SetColumns(0);
SetMemoryLimit(12*10^9);
AttachSpec("/home/claude/Magma/magma.spec");
Attach("/home/claude/torsion_jac/genus2.m");
Q := Rationals();
Px<x> := PolynomialRing(Q);

// ---------- (A) validation over Q ----------
// two curves with full rational 2-torsion (trivially 2-congruent)
fA := (x)*(x-1)*(x-4);
fB := (x-2)*(x+3)*(x+7);
LA := Genus2Elliptic2(fA, fB);
printf "(A) Genus2Elliptic2: %o curves\n", #LA;
invA := { IgusaClebschInvariants(C) : C in LA };

// my transplant, base field Q, splitting field = Q handled via L := Q?
// (roots rational: use L = RationalField-as-number-field trick: work in Q directly)
alpha := [ r[1] : r in Roots(fA) ];
gamma := [ r[1] : r in Roots(fB) ];
Df := Discriminant(fA); Dg := Discriminant(fB);
a := [[ alpha[i]-alpha[j] : j in [1..3]] : i in [1..3]];
Rt<t> := PolynomialRing(Q);
myinv := {};
for sgm in SymmetricGroup(3) do
    beta := [ gamma[i^sgm] : i in [1..3] ];
    b := [[ beta[i]-beta[j] : j in [1..3]] : i in [1..3]];
    if alpha[1]*b[3][2] + alpha[2]*b[1][3] + alpha[3]*b[2][1] ne 0 then
        a1 := a[3][2]^2/b[3][2]+a[2][1]^2/b[2][1]+a[1][3]^2/b[1][3];
        a2 := alpha[1]*b[3][2]+alpha[2]*b[1][3]+alpha[3]*b[2][1];
        b1 := b[3][2]^2/a[3][2]+b[2][1]^2/a[2][1]+b[1][3]^2/a[1][3];
        b2 := beta[1]*a[3][2]+beta[2]*a[1][3]+beta[3]*a[2][1];
        A := Dg*a1/a2;  B := Df*b1/b2;
        sext := -(A*a[2][1]*a[1][3]*t^2+B*b[2][1]*b[1][3])
                *(A*a[3][2]*a[2][1]*t^2+B*b[3][2]*b[2][1])
                *(A*a[1][3]*a[3][2]*t^2+B*b[1][3]*b[3][2]);
        if Discriminant(sext) ne 0 then
            Include(~myinv, IgusaClebschInvariants(HyperellipticCurve(sext)));
        end if;
    end if;
end for;
printf "(A) invariant sets equal: %o (mine %o, theirs %o)\n",
    invA eq myinv, #myinv, #invA;

// ---------- (B) the K-glue at s0 = 4/5 ----------
K<w> := QuadraticField(11);
RK<xk> := PolynomialRing(K);
s0 := K!(4/5);
r0 := Roots(xk^2 - (s0^3-3*s0^2+4*s0)*xk + s0)[1][1];
c0 := s0*(r0-1); b0 := r0*c0;
E := EllipticCurve([1-c0, -b0, -b0, 0, 0]);
sig := hom< K -> K | -w >;
Es := EllipticCurve([sig(a) : a in aInvariants(E)]);
f := HyperellipticPolynomials(WeierstrassModel(E));
g := HyperellipticPolynomials(WeierstrassModel(Es));
L := SplittingField(f);
RL := PolynomialRing(L);
al := [ r[1] : r in Roots(RL!f) ];
gm := [ r[1] : r in Roots(RL!g) ];
DfK := Discriminant(f); DgK := Discriminant(g);
GA := Automorphisms(L);
GK := [ gA : gA in GA | gA(L!w) eq L!w ];
Gaction := func< rs | [[ Index(rs, gA(r)) : r in rs ] : gA in GK ] >;
aK := [[ al[i]-al[j] : j in [1..3]] : i in [1..3]];
RtL<tL> := PolynomialRing(L);
CK := 0; found := false;
for sgm in SymmetricGroup(3) do
    beta := [ gm[i^sgm] : i in [1..3] ];
    b := [[ beta[i]-beta[j] : j in [1..3]] : i in [1..3]];
    if al[1]*b[3][2] + al[2]*b[1][3] + al[3]*b[2][1] ne 0 and Gaction(al) eq Gaction(beta) then
        a1 := aK[3][2]^2/b[3][2]+aK[2][1]^2/b[2][1]+aK[1][3]^2/b[1][3];
        a2 := al[1]*b[3][2]+al[2]*b[1][3]+al[3]*b[2][1];
        b1 := b[3][2]^2/aK[3][2]+b[2][1]^2/aK[2][1]+b[1][3]^2/aK[1][3];
        b2 := beta[1]*aK[3][2]+beta[2]*aK[1][3]+beta[3]*aK[2][1];
        A := DgK*a1/a2;  B := DfK*b1/b2;
        sext := -(A*aK[2][1]*aK[1][3]*tL^2+B*b[2][1]*b[1][3])
                *(A*aK[3][2]*aK[2][1]*tL^2+B*b[3][2]*b[2][1])
                *(A*aK[1][3]*aK[3][2]*tL^2+B*b[1][3]*b[3][2]);
        sextK := RK![ K!c : c in Coefficients(sext) ];
        CK := HyperellipticCurve(sextK);
        found := true;
        break;
    end if;
end for;
error if not found, "no matching";
printf "(B) CK built over K\n";
OK := Integers(K);
ntested := 0;
for p in PrimesInInterval(13, 200) do
    if ntested ge 6 then break; end if;
    dec := Decomposition(OK, p);
    if #dec ne 2 then continue; end if;
    pr := dec[1][1];
    kp, red := ResidueClassField(pr);
    okr := true; np1 := 0; np2 := 0; nj := 0;
    try
        E1p := Reduction(E, pr); np1 := #E1p;
        E2p := Reduction(Es, pr); np2 := #E2p;
        fp := PolynomialRing(kp)![ red(c) : c in Coefficients(HyperellipticPolynomials(CK)) ];
        if Degree(fp) ge 5 and IsSquarefree(fp) then
            nj := #Jacobian(HyperellipticCurve(fp));
        else
            okr := false;
        end if;
    catch e;
        okr := false;
    end try;
    if not okr then continue; end if;
    ntested +:= 1;
    printf "(B) p=%o: #J(CK)=%o  #E*#Es=%o  match=%o  11|J: %o\n",
        p, nj, np1*np2, nj eq np1*np2, nj mod 11 eq 0;
end for;
printf "S45CHECK_DONE\n";
quit;
