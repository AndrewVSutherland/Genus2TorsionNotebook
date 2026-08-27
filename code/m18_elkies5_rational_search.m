//////////////////////////////////////////////////////////////////////
//  Low-height exact search for M_1(8) plus rational 5-torsion.
//
//      C: y^2 = q(x) * (x^4 + q(x)),
//      q = a*x^2 + b*x + c.
//
//  The finite-field sieve gives a nonboundary obstruction modulo 11, so
//  any rational hit with full rational 2-torsion should have bad reduction
//  at 11.  This script nevertheless checks low-height primitive integral
//  triples directly.
//
//  Typical run from torsion_jac:
//      magma -b B:=30 code/m18_elkies5_rational_search.m
//////////////////////////////////////////////////////////////////////

if not assigned B then
    B := 30;
elif Type(B) eq MonStgElt then
    B := StringToInteger(B);
end if;
if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;

Qx<x> := PolynomialRing(Rationals());

function IsPrimitiveTriple(a, b, c)
    return GCD([Integers()!a, Integers()!b, Integers()!c]) eq 1;
end function;

function IsFullSplitSextic(f)
    if Degree(f) ne 6 or Discriminant(f) eq 0 then
        return false;
    end if;
    fac := Factorization(f);
    return #fac eq 6 and &and [ Degree(ff[1]) eq 1 and ff[2] eq 1 : ff in fac ];
end function;

checked := 0;
smooth := 0;
full_split := 0;
bad11 := 0;
torsion_tests := 0;
hits := [];

for a in [-B..B] do
    for b in [-B..B] do
        for c in [-B..B] do
            if a eq 0 and b eq 0 and c eq 0 then
                continue;
            end if;
            if a lt 0 then
                continue;
            end if;
            if a eq 0 and b lt 0 then
                continue;
            end if;
            if a eq 0 and b eq 0 and c lt 0 then
                continue;
            end if;
            if not IsPrimitiveTriple(a,b,c) then
                continue;
            end if;

            checked +:= 1;
            q := a*x^2 + b*x + c;
            f := q * (x^4 + q);
            if Degree(f) ne 6 or Discriminant(f) eq 0 then
                continue;
            end if;
            smooth +:= 1;
            if not IsFullSplitSextic(f) then
                continue;
            end if;
            full_split +:= 1;
            if Integers()!Numerator(Discriminant(f)) mod 11 eq 0 then
                bad11 +:= 1;
            end if;

            C := HyperellipticCurve(f);
            J := Jacobian(C);
            torsion_tests +:= 1;
            G, phi := TorsionSubgroup(J);
            invs := Invariants(G);
            has5 := &or [ n mod 5 eq 0 : n in invs ];
            print "split", [a,b,c], "disc_mod_11",
                  Integers()!Numerator(Discriminant(f)) mod 11,
                  "torsion", invs;
            if has5 then
                Append(~hits, <[a,b,c], invs, f>);
                print "HIT", [a,b,c], invs, f;
                if #hits ge max_hits then
                    break a;
                end if;
            end if;
        end for;
    end for;
end for;

print "DONE B", B;
print "checked", checked, "smooth", smooth, "full_split", full_split,
      "bad11_full_split", bad11, "torsion_tests", torsion_tests, "hits", #hits;

quit;
