// Lane 8 (overnight 2026-07-25): Howe's order-70 curve f70.
// Exact torsion, factor type, RM/split screen, strict two-prime simplicity certificate.
// Source of f70 in repo: code/claude_c35_validate.m line 47 (literature intake 2026-07-23).
SetColumns(0);
SetMemoryLimit(3*10^9);

QQ := Rationals();
R<x> := PolynomialRing(QQ);

f70 := 22*x^5 + (697/144)*x^4 - (645/4)*x^3 + (1045/4)*x^2 - 162*x + 36;
print "f70 =", f70;
print "degree", Degree(f70);

// integral model: y -> y/L, f -> L^2 f  with L = 12 (denominators 144,4)
L := 12;
g := R!(L^2*f70);
print "integral model g = 144*f70 =", g;
print "squarefree:", IsSquarefree(g);

// factor type (degree-only use of Factorization is safe)
fac := Factorization(g);
unit := g div &*[t[1]^t[2] : t in fac];
assert unit*&*[t[1]^t[2] : t in fac] eq g;
print "leading unit:", unit;
print "factor type:", [<Degree(t[1]),t[2]> : t in fac];
for t in fac do print "  factor:", t[1], "mult", t[2]; end for;

C := HyperellipticCurve(g);
print "genus:", Genus(C);
J := Jacobian(C);
print "conductor-ish disc of g:", Factorization(Integers()!Discriminant(g));

t0 := Cputime();
T, mp := TorsionSubgroup(J);
print "TORSION invariants:", Invariants(T);
print "torsion order:", #T;
print "torsion time:", Cputime(t0);

// 2-rank
T2 := TwoTorsionSubgroup(J);
print "2-torsion invariants:", Invariants(T2);

// generators of the torsion, and their orders
for i in [1..#Generators(T)] do
    gg := T.i;
    print "gen", i, "order", Order(gg), "->", mp(gg);
end for;

// ---------------- simplicity / split screen ----------------
// scan good primes: L-polynomial factorization, root-power strictness, real-subfield disc
Rz<z> := PolynomialRing(QQ);
Rxy<X,Y> := PolynomialRing(QQ,2);

badp := Set(PrimeDivisors(Integers()!Discriminant(g)) cat [2]);
print "bad primes (incl 2):", Sort(Setseq(badp));

sigs := [];
nirred := 0;
for p in PrimesInInterval(3,300) do
    if p in badp then continue; end if;
    Cp := ChangeRing(C, GF(p));
    if not IsNonsingular(Cp) then continue; end if;
    Lp := LPolynomial(Cp);
    cs := Coefficients(Lp);       // L(T) = 1 + a1 T + a2 T^2 + p a1 T^3 + p^2 T^4
    chi := Rz!Reverse(cs);        // charpoly of Frobenius, monic quartic
    c3 := Coefficient(chi,3); c2 := Coefficient(chi,2);
    // real subfield: y^2 - c3 y + (c2 - 2p)
    dd := Integers()!(c3^2 - 4*(c2-2*p));
    if dd eq 0 then sc := 0; else sc := Squarefree(dd); end if;
    ir := IsIrreducible(chi);
    // split screen: factorization type of chi over Q
    ftype := [Degree(t[1]) : t in Factorization(chi)];
    Append(~sigs, <p, sc, ir, ftype>);
    if ir then nirred +:= 1; end if;
end for;
print "primes tested:", #sigs, " irreducible chi count:", nirred;
print "real-subfield squarefree cores:", Sort(Setseq({s[2] : s in sigs}));
print "per-prime <p, sqfree core, chi irreducible, chi factor degrees>:", sigs;

// ---------------- explicit split test: identify the elliptic factors ----------------
// If chi = (x^2-a1 x+p)(x^2-a2 x+p) at every good prime, J ~_Q E1 x E2.
// Identify E1,E2 by matching the a_p sequence against curves of conductor N | 2^a 3^b 11^c 13^d.
print "=== elliptic subcover identification ===";
aps := AssociativeArray();
for s in sigs do
    p := s[1];
    Cp := ChangeRing(C, GF(p));
    chi := Rz!Reverse(Coefficients(LPolynomial(Cp)));
    ff := Factorization(chi);
    if &and[Degree(t[1]) eq 2 : t in ff] and #ff+ (#ff eq 1 select 1 else 0) ge 1 then
        qs := &cat[[t[1] : i in [1..t[2]]] : t in ff];
        aps[p] := Sort([-Coefficient(q,1) : q in qs]);
    end if;
end for;
print "traces (a1,a2) at good primes:", [<p, aps[p]> : p in Sort(Setseq(Keys(aps)))];

D := EllipticCurveDatabase();
print "Cremona database conductor limit:", LargestConductor(D);
cands := [];
for N in [n : n in [1..LargestConductor(D)] | Set(PrimeDivisors(n)) subset {2,3,11,13}] do
    for k in [1..NumberOfIsogenyClasses(D,N)] do
        E := EllipticCurve(D,N,k,1);
        ok := true;
        for p in Sort(Setseq(Keys(aps)))[1..Min(10,#Keys(aps))] do
            ap := TraceOfFrobenius(ChangeRing(E,GF(p)));
            if not ap in aps[p] then ok := false; break; end if;
        end for;
        if ok then Append(~cands, <N,k,Invariants(TorsionSubgroup(E)),aInvariants(E)>); end if;
    end for;
end for;
print "matching elliptic isogeny classes:", cands;

quit;

quit;
