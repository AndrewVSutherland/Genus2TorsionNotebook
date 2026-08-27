// claude_ov_b4_c43close.m -- Lane 4 (route B4): CLOSE the (4,3) component.
//
// On the (4,3) component of the quadratic-factor incidence of Flynn's order-11
// family, with the clean parametrization
//        u(w) = (2w^3-2w^2+2w-1)/w^2,  v(w) = ((w^2-w+1)/w)^2,
//        t(w) = -((w-1)(w^2-w+1)/w^2)^2,
// the sextic F_{t(w)} has GENERIC FACTOR TYPE [1,2,3] over Q(w):
//        F = L(x) * G1(x) * C(x),  deg 1 / 2 / 3.
// Hence generic 2-rank 1 (as on the (6,2) stream), and J(Q)[2] has rank >= 2
// for a member iff EITHER
//   (i)  G1 splits   <=>  disc(G1) = -(4w^3-4w^2+4w-1)/w^4 is a square, i.e.
//        a rational point on  E19 : y^2 = -4w^3+4w^2-4w+1   (cond 19, rank 0);
//   (ii) C acquires a rational root  <=>  a rational point on the plane curve
//        C3 : C(w,x) = 0.
// Moreover a Galois-stable (2,2)-kernel EXISTS at all only in case (ii): the
// quartic cofactor L*C has one rational root, so a stable pairing must pair it
// with a Galois-fixed root of C.  So on this component "Richelot is possible"
// == "2-rank is already 2"; the Richelot step buys nothing.
//
// This script certifies both loci.
//
// Run: code/claude_magma_slot.sh -b MemGB:=12 code/claude_ov_b4_c43close.m \
//        > results/claude_ov_b4_c43close.log 2>&1 &

SetColumns(0);
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
if not assigned SearchB then SearchB := 100000; elif Type(SearchB) eq MonStgElt then SearchB := StringToInteger(SearchB); end if;

Q := Rationals();

// ---------------------------------------------------------------- locus (i)
printf "=== LOCUS (i): the quadratic block splits ===\n";
E := EllipticCurve([0,4,0,16,16]);          // Y^2 = X^3+4X^2+16X+16, X=-4w, Y=4y
E := MinimalModel(E);
printf "E19 minimal model: %o  conductor %o\n", aInvariants(E), Conductor(E);
printf "E19 rank bounds: %o\n", [RankBounds(E)];
printf "E19 torsion: %o\n", Invariants(TorsionSubgroup(E));
MW, mp, ok1, ok2 := MordellWeilGroup(E);
printf "E19 MordellWeilGroup: %o  proven: %o %o\n", Invariants(MW), ok1, ok2;
printf "E19 rational points: %o\n", [mp(g) : g in MW];
CE := HyperellipticCurve(-4*PolynomialRing(Q).1^3+4*PolynomialRing(Q).1^2-4*PolynomialRing(Q).1+1);
printf "direct model y^2=-4w^3+4w^2-4w+1 : points to 10^6 = %o\n", Points(CE : Bound := 1000000);

// ------------------------------------------------------------ the parametrization
Fw<w> := FunctionField(Q);
uu := (2*w^3 - 2*w^2 + 2*w - 1)/w^2;
vv := ((w^2 - w + 1)/w)^2;
tS := -((w-1)*(w^2-w+1)/w^2)^2;
PK<X> := PolynomialRing(Fw);
F := X^6 + 2*X^5 + (2*tS+3)*X^4 + 2*X^3 + (tS^2+1)*X^2 + 2*tS*(1-tS)*X + tS^2;
fac := Factorization(F);
Cub := [f[1] : f in fac | Degree(f[1]) eq 3][1];
Lin := [f[1] : f in fac | Degree(f[1]) eq 1][1];
printf "\nLINEAR FACTOR: %o\nCUBIC FACTOR: %o\n", Lin, Cub;

// --------------------------------------------------------------- locus (ii)
printf "\n=== LOCUS (ii): the cubic factor acquires a rational root ===\n";
P2<W,Xv> := PolynomialRing(Q, 2);
den := LCM([Denominator(Coefficient(Cub,i)) : i in [0..3]]);
cc := [P2!Evaluate(Numerator(den*Coefficient(Cub,i)), W) : i in [0..3]];
C3eq := &+[cc[i+1]*Xv^i : i in [0..3]];
printf "C3 equation: %o\n", C3eq;
printf "C3 bidegree (w,x) = (%o,%o)\n", Degree(C3eq, W), Degree(C3eq, Xv);
printf "C3 absolutely irreducible: %o\n", IsAbsolutelyIrreducible(Curve(AffineSpace(P2), C3eq));
A2 := AffineSpace(P2);
C3 := Curve(A2, C3eq);
PC3 := ProjectiveClosure(C3);
g3 := Genus(PC3);
printf "C3 GEOMETRIC GENUS = %o\n", g3;
pts := PointSearch(PC3, 2000);
printf "C3 PointSearch(2000): %o points\n", #pts;
for p in pts do printf "  C3PT %o\n", p; end for;
bh := IsHyperelliptic(PC3);
printf "C3 hyperelliptic: %o\n", bh;

// --- the SYMMETRIC model and the bielliptic quotient ----------------------
printf "\n--- symmetric model  T(w) + T(z) = 0  (z = 1/(1-rho)) ---\n";
P2s<Ws,Zs> := PolynomialRing(Q, 2);
Dsym := (Ws^3-2*Ws^2+2*Ws-1)*Zs^2 + (Zs^3-2*Zs^2+2*Zs-1)*Ws^2;
Csym := Curve(AffineSpace(P2s), Dsym);
printf "Dsym bidegree (%o,%o), absolutely irreducible: %o\n",
    Degree(Dsym,Ws), Degree(Dsym,Zs), IsAbsolutelyIrreducible(Csym);
printf "Dsym GENUS = %o\n", Genus(ProjectiveClosure(Csym));
// quotient by the swap: p = w+z, q = wz  =>  p^2 - p(q^2+2q) + 4q^2-2q = 0
Pq<qv> := PolynomialRing(Q);
Cq := HyperellipticCurve(qv^4 + 4*qv^3 - 12*qv^2 + 8*qv);
printf "QUOTIENT genus-1 curve y^2 = q^4+4q^3-12q^2+8q, genus %o\n", Genus(Cq);
Eq, mq := EllipticCurve(Cq, Cq![0,0,1]);
Eq2, m2 := MinimalModel(Eq);
printf "  MINIMAL MODEL %o  conductor %o\n", aInvariants(Eq2), Conductor(Eq2);
printf "  rank bounds %o  torsion %o\n", [RankBounds(Eq2)], Invariants(TorsionSubgroup(Eq2));
MWq, mpq, o1, o2 := MordellWeilGroup(Eq2);
printf "  MordellWeilGroup %o proven %o %o  generators %o\n",
    Invariants(MWq), o1, o2, [mpq(g) : g in Generators(MWq)];
E62 := EllipticCurve([0,3,0,2,1]);
printf "  SAME CURVE as the (6,2) component's model (92.a1)? %o\n",
    IsIsomorphic(Eq2, MinimalModel(E62));

// ------------------------------------------------- membership + torsion probe
printf "\n=== MEMBERS: exact torsion of y^2 = F_{t(w)} for small w ===\n";
Pol<x> := PolynomialRing(Q);
function Member(wv)
    tv := -((wv-1)*(wv^2-wv+1)/wv^2)^2;
    return x^6 + 2*x^5 + (2*tv+3)*x^4 + 2*x^3 + (tv^2+1)*x^2 + 2*tv*(1-tv)*x + tv^2, tv;
end function;
function IntegralSextic(f)
    // scale x -> X/m, y -> Y/m^3 so the model is integral (TorsionSubgroup trap)
    m := 1;
    cs := Coefficients(f);           // cs[i+1] = coefficient of x^i
    for pr in &join[{p : p in PrimeDivisors(Denominator(c))} : c in cs | c ne 0] do
        need := 0;
        for i in [0..5] do
            e := Valuation(cs[i+1], pr);
            if e lt 0 then need := Max(need, Ceiling(-e/(6-i))); end if;
        end for;
        m := m * pr^need;
    end for;
    g := m^6 * Evaluate(f, x/m);
    assert &and[IsIntegral(c) : c in Coefficients(g)];
    return g, m;
end function;
for wv in [2, 3, -1, 1/2, 4, -2, 5, 3/2, -3, 2/3, 7, 5/2] do
    f, tv := Member(wv);
    if Discriminant(f) eq 0 then printf "w=%o : SINGULAR\n", wv; continue; end if;
    ft := [Degree(q[1]) : q in Factorization(f)];
    g, m := IntegralSextic(f);
    C := HyperellipticCurve(g);
    J := Jacobian(C);
    T := TorsionSubgroup(J);
    printf "w=%o t=%o factortype=%o TORSION=%o cond-disc=%o\n",
        wv, tv, Sort(ft), Invariants(T), Factorization(Integers()!Discriminant(C));
end for;

printf "C43CLOSE_DONE\n";
quit;
