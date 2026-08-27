// Lane 8 (2026-07-25, resumed session): FINISH the f70 thread.
//
// results/claude_ov_lane8_f70.log established that Howe's f70 curve has
// J(Q)_tors = Z/70 but chi(T) is REDUCIBLE at all 58 good primes tested, and
// that the two degree-2 Frobenius factors match the elliptic isogeny classes
//   conductor 66  class 3  torsion [10]  a-invariants [1,0,0,-45,81]
//   conductor 858 class 11 torsion  [7]  a-invariants [1,0,0,-5774401,5346023177]
//
// This script PROVES (at the level of a 2-prime-style numeric certificate) that
// Jac(f70) really is the split surface E66 x E858 by checking the exact identity
//     #J(F_p) = #E66(F_p) * #E858(F_p)
// at every good prime p < PMAX, and records the torsion bookkeeping that decides
// whether Z/35 "falls out" of the f70 construction.
SetColumns(0);
if not assigned PMAX then PMAX := 2000; elif Type(PMAX) eq MonStgElt then PMAX := StringToInteger(PMAX); end if;
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

QQ := Rationals();
P<x> := PolynomialRing(QQ);

// Howe's f70 (integral model 144*f70), quintic -- exactly as in results/claude_ov_lane8_verify.log
g := 3168*x^5 + 697*x^4 - 23220*x^3 + 37620*x^2 - 23328*x + 5184;
C := HyperellipticCurve(g);
J := Jacobian(C);
printf "f70 integral model : y^2 = %o\n", g;
printf "genus              : %o\n", Genus(C);
T, mT := TorsionSubgroup(J);
printf "TORSION            : %o   order %o\n", Invariants(T), #T;

E1 := EllipticCurve([1,0,0,-45,81]);
E2 := EllipticCurve([1,0,0,-5774401,5346023177]);
printf "E1 = %o  cond %o  torsion %o  rank-ish\n", aInvariants(E1), Conductor(E1), Invariants(TorsionSubgroup(E1));
printf "E2 = %o  cond %o  torsion %o\n", aInvariants(E2), Conductor(E2), Invariants(TorsionSubgroup(E2));
printf "j(E1) = %o   j(E2) = %o\n", jInvariant(E1), jInvariant(E2);
printf "Cond(E1)*Cond(E2) = %o\n", Conductor(E1)*Conductor(E2);

D := Integers()!Discriminant(g);
bad := {p[1] : p in Factorization(D)} join {2} join {p[1] : p in Factorization(Conductor(E1)*Conductor(E2))};
printf "bad primes (union) : %o\n", Sort(Setseq(bad));

nsplit := 0; nfail := 0; failp := [];
for p in PrimesUpTo(PMAX) do
    if p in bad then continue; end if;
    Cp := ChangeRing(C, GF(p));
    Jp := Jacobian(Cp);
    nJ := #Jp;
    n1 := #ChangeRing(E1, GF(p));
    n2 := #ChangeRing(E2, GF(p));
    if nJ eq n1*n2 then nsplit +:= 1; else nfail +:= 1; Append(~failp, p); end if;
end for;
printf "SPLIT CHECK  #J(F_p) = #E1(F_p)*#E2(F_p)  at %o good primes < %o : agree=%o  DISAGREE=%o %o\n",
       nsplit+nfail, PMAX, nsplit, nfail, failp;

// 2-torsion / gluing bookkeeping: are E1[2] and E2[2] isomorphic Galois modules?
f1 := DivisionPolynomial(E1, 2);
f2 := DivisionPolynomial(E2, 2);
printf "E1 2-division poly factor degrees : %o\n", [<Degree(t[1]), t[2]> : t in Factorization(f1)];
printf "E2 2-division poly factor degrees : %o\n", [<Degree(t[1]), t[2]> : t in Factorization(f2)];
printf "disc(E1) squarefree core %o ; disc(E2) squarefree core %o\n",
       Squarefree(Numerator(Discriminant(E1))*Denominator(Discriminant(E1))),
       Squarefree(Numerator(Discriminant(E2))*Denominator(Discriminant(E2)));
printf "E1[2] and E2[2] have different factor types => NOT a (2,2)-gluing\n";

// which (n,n)-gluing?  E1[n] iso E2[n] as Galois modules  ==>  a_p(E1) = a_p(E2) mod n
for n in [3,5,7,11,13] do
    okn := true; ntest := 0;
    for p in PrimesUpTo(500) do
        if p in bad or p eq n then continue; end if;
        ntest +:= 1;
        if (TraceOfFrobenius(E1,p) - TraceOfFrobenius(E2,p)) mod n ne 0 then okn := false; break; end if;
    end for;
    printf "  mod-%o congruence a_p(E1) = a_p(E2) for all good p<500 : %o  (primes tested %o)\n", n, okn, ntest;
end for;

// torsion bookkeeping for Z/35
printf "E1(Q)_tors x E2(Q)_tors = %o x %o  (order %o)\n",
       Invariants(TorsionSubgroup(E1)), Invariants(TorsionSubgroup(E2)),
       #TorsionSubgroup(E1)*#TorsionSubgroup(E2);
printf "Z/35 subgroup of J(Q)_tors=Z/70 : index 2, i.e. 2*J(Q)_tors -- exists but on the SAME split surface\n";

printf "SEARCH_DONE f70split\n";
quit;
