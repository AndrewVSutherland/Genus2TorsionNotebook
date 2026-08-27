// Lane 8 (overnight 2026-07-25): the two HLP "P^0" split models, transcribed from
// Howe-Leprevost-Poonen, Forum Math. 12 (2000) 315-364, Section 3.6:
//   eq (4), p.15 :  Z/63     C : y^2 = 897x^6 - 197570x^4 + 79136353x^2 - 146398496
//   p.15         :  Z/7xZ/7  C : y^2 = x^6 + 3025x^4 + 3232987x^2 + 869675859
// Verify: genus, exact TorsionSubgroup, 2-rank, split/simple screen over good primes,
// and identify the elliptic subcovers (both curves are even in x, so x -> x^2 gives
// the two degree-2 elliptic quotients explicitly).
SetColumns(0);
SetMemoryLimit(3*10^9);
QQ := Rationals();
R<x> := PolynomialRing(QQ);
Rz<z> := PolynomialRing(QQ);

procedure Report(name, f)
    print "==================================================";
    print name;
    print "f =", f;
    print "degree", Degree(f), "squarefree", IsSquarefree(f);
    fac := Factorization(f);
    unit := f div &*[t[1]^t[2] : t in fac];
    assert unit*&*[t[1]^t[2] : t in fac] eq f;
    print "leading unit", unit, "factor degrees", [Degree(t[1]) : t in fac];
    C := HyperellipticCurve(f);
    print "genus", Genus(C);
    J := Jacobian(C);
    dsc := Integers()!Discriminant(f);
    print "disc factorization:", Factorization(dsc);
    t0 := Cputime();
    T, mp := TorsionSubgroup(J);
    print "TORSION invariants:", Invariants(T), " order", #T, " (", Cputime(t0), "s )";
    print "2-torsion invariants:", Invariants(TwoTorsionSubgroup(J));
    for i in [1..#Generators(T)] do
        print "  gen", i, "order", Order(T.i), "->", mp(T.i);
    end for;
    // even-in-x => explicit elliptic quotients y^2 = fE(u), u = x^2, and the twist
    cs := Coefficients(f);
    if &and[cs[i] eq 0 : i in [2..#cs by 2]] then
        c := [cs[1],cs[3],cs[5],cs[7]];
        E1 := EllipticCurve(z^3*c[4] + z^2*c[3] + z*c[2] + c[1]);   // v^2 = f(u), u=x^2
        E2 := EllipticCurve(z^3*c[1] + z^2*c[2] + z*c[3] + c[4]);   // reversed: u=1/x^2
        print "  SUBCOVER E1 (u=x^2):", aInvariants(E1), " cond", Conductor(E1),
              " torsion", Invariants(TorsionSubgroup(E1)), " rank?", Rank(E1);
        print "  SUBCOVER E2 (u=1/x^2):", aInvariants(E2), " cond", Conductor(E2),
              " torsion", Invariants(TorsionSubgroup(E2)), " rank?", Rank(E2);
    end if;
    // simplicity / split screen
    badp := Set(PrimeDivisors(dsc)) join {2};
    sigs := []; nirred := 0;
    for p in PrimesInInterval(3,300) do
        if p in badp then continue; end if;
        Cp := ChangeRing(C,GF(p));
        if not IsNonsingular(Cp) then continue; end if;
        chi := Rz!Reverse(Coefficients(LPolynomial(Cp)));
        c3 := Coefficient(chi,3); c2 := Coefficient(chi,2);
        dd := Integers()!(c3^2-4*(c2-2*p));
        sc := dd eq 0 select 0 else Squarefree(dd);
        ir := IsIrreducible(chi);
        if ir then nirred +:= 1; end if;
        Append(~sigs, <p, sc, ir, [Degree(t[1]) : t in Factorization(chi)]>);
    end for;
    print "good primes tested", #sigs, " chi irreducible count", nirred;
    print "real-subfield squarefree cores:", Sort(Setseq({s[2] : s in sigs}));
    print "sig detail:", sigs;
end procedure;

Report("HLP Z/63  (Forum Math 12 (2000), eq. (4))",
       897*x^6 - 197570*x^4 + 79136353*x^2 - 146398496);
Report("HLP Z/7 x Z/7 (Forum Math 12 (2000), p.15)",
       x^6 + 3025*x^4 + 3232987*x^2 + 869675859);
quit;
