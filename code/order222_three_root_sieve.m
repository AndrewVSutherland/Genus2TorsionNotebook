//////////////////////////////////////////////////////////////////////
// Search for non-RM [2,22] candidates in the three-rational-root chart.
//
// Every curve in this chart has
//     f = x (x-1) (x-r) (x^3 + a*x^2 + b*x + c),
// hence exactly two independent rational 2-torsion classes when the cubic
// is irreducible.  Rational torsion injects into J(F_p) at every
// prime of good reduction, so a running gcd of #J(F_p) is a very cheap
// necessary filter for 44-torsion.
//
// Usage:
//   magma -b B:=6 exact_cap:=25 code/order222_three_root_sieve.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
P<x> := PolynomialRing(Q);
B := assigned B select StringToInteger(B) else 8;
exact_cap := assigned exact_cap select StringToInteger(exact_cap) else 25;
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
// Backup primes 23..41 engage only when small primes divide the disc
// (the B=6 run's sole false positive survived on a 3-prime gcd).
primes := [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41];
// Optional r-shard bounds for parallel runs: magma -b B:=12 rlo:=-12 rhi:=-5 ...
rlo := assigned rlo select StringToInteger(rlo) else -B;
rhi := assigned rhi select StringToInteger(rhi) else B;

function Primitive(f)
    den := LCM([Denominator(Coefficient(f, i)) : i in [0..Degree(f)]]);
    F := den*f;
    cont := GCD([Integers()!Coefficient(F, i) : i in [0..Degree(F)]]);
    return F/cont;
end function;

function ReductionGCD44(f)
    g := 0;
    used := 0;
    for p in primes do
        fp := ChangeRing(f, GF(p));
        if Discriminant(fp) eq 0 then continue; end if;
        C := HyperellipticCurve(fp);
        n := Integers()!Evaluate(LPolynomial(C), 1);
        g := GCD(g, n);
        used +:= 1;
        if g mod 44 ne 0 then return false, g, used; end if;
    end for;
    return used ge 3 and g mod 44 eq 0, g, used;
end function;

procedure FrobeniusCertificate(f)
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43] do
        fp := ChangeRing(f, GF(p));
        if Discriminant(fp) eq 0 then continue; end if;
        Lp := LPolynomial(HyperellipticCurve(fp));
        fac := Factorization(Lp);
        if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
            print "IRREDUCIBLE_FROBENIUS", p, Lp;
            return;
        end if;
    end for;
    print "NO_IRREDUCIBLE_FROBENIUS_IN_LIST";
end procedure;

tested := 0;
passed := 0;
exact := 0;

for r in [rlo..rhi] do
    if r eq 0 or r eq 1 then continue; end if;
    print "PROGRESS", "r", r;
    for c in [-B..B] do
        if c eq 0 then continue; end if;
        for b in [-B..B] do
            for a in [-B..B] do
                q := x^3 + a*x^2 + b*x + c;
                // Keep the intended 1+1+1+3 orbit chart.
                if not IsIrreducible(q) then continue; end if;
                f := x*(x-1)*(x-r)*q;
                if Discriminant(f) eq 0 then continue; end if;
                tested +:= 1;
                ok, g, used := ReductionGCD44(f);
                if not ok then continue; end if;
                passed +:= 1;
                print "GCD44_SURVIVOR", <r,a,b,c>, "gcd", g,
                      "good_primes", used;
                if exact ge exact_cap then continue; end if;
                exact +:= 1;
                F := Primitive(f);
                C := HyperellipticCurve(F);
                J := Jacobian(C);
                T, mp := TorsionSubgroup(J);
                inv := Invariants(T);
                print "EXACT", <r,a,b,c>, "torsion", inv, "f", F;
                if inv eq [2,22] then
                    print "TARGET_2_22", <r,a,b,c>, F;
                    FrobeniusCertificate(F);
                end if;
            end for;
        end for;
    end for;
end for;

print "SUMMARY", "B", B, "tested", tested, "gcd44", passed,
      "exact", exact;
print "SEARCH_DONE";
quit;
