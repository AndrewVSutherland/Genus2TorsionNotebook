//////////////////////////////////////////////////////////////////////
//  Root-parametrized exact search for M_1(8) plus rational 5-torsion.
//
//  Write
//
//      x^4 + q(x) = (x-r)(x-s)(x-t)(x+r+s+t),
//
//  and require q = a*x^2 + b*x + c to split over Q.  This directly
//  searches the full-rational-Weierstrass locus inside M_1(8).
//
//  Because the finite-field sieve obstructs nonboundary mod 11, exact
//  torsion is only tested when the sextic discriminant is 0 mod 11.
//
//  Typical run from torsion_jac:
//      magma -b H:=40 code/m18_elkies5_split_root_search.m
//////////////////////////////////////////////////////////////////////

if not assigned H then
    H := 40;
elif Type(H) eq MonStgElt then
    H := StringToInteger(H);
end if;
if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;

Z := Integers();
Qx<x> := PolynomialRing(Rationals());

function IsSquareInteger(n)
    if n lt 0 then
        return false;
    end if;
    ok, root := IsSquare(n);
    return ok;
end function;

seen := {};
root_tuples := 0;
smooth_split := 0;
q_split := 0;
bad11 := 0;
torsion_tests := 0;
hits := [];

for r in [-H..H] do
    for s in [-H..H] do
        for t in [-H..H] do
            u := -r-s-t;
            if Abs(u) gt H then
                continue;
            end if;
            roots := Sort([r,s,t,u]);
            if #(Set(roots)) ne 4 then
                continue;
            end if;
            key := <roots[1], roots[2], roots[3], roots[4]>;
            if key in seen then
                continue;
            end if;
            Include(~seen, key);
            root_tuples +:= 1;

            quartic := &*[ x - Rationals()!rr : rr in roots ];
            coeffs := Coefficients(quartic);
            c := Z!coeffs[1];
            b := Z!coeffs[2];
            a := Z!coeffs[3];
            // coeffs[4] is the x^3 coefficient, forced to be zero.
            if a eq 0 then
                continue;
            end if;
            q := a*x^2 + b*x + c;
            f := q * quartic;
            if Degree(f) ne 6 or Discriminant(f) eq 0 then
                continue;
            end if;
            smooth_split +:= 1;

            delta_q := b^2 - 4*a*c;
            if not IsSquareInteger(delta_q) then
                continue;
            end if;
            q_split +:= 1;

            disc_num := Z!Numerator(Discriminant(f));
            if disc_num mod 11 ne 0 then
                continue;
            end if;
            bad11 +:= 1;

            C := HyperellipticCurve(f);
            J := Jacobian(C);
            torsion_tests +:= 1;
            G, phi := TorsionSubgroup(J);
            invs := Invariants(G);
            print "test roots", roots, "q", [a,b,c], "torsion", invs;
            if &or [ n mod 5 eq 0 : n in invs ] then
                Append(~hits, <roots, [a,b,c], invs, f>);
                print "HIT roots", roots, "q", [a,b,c], "torsion", invs, "f", f;
                if #hits ge max_hits then
                    break r;
                end if;
            end if;
        end for;
    end for;
end for;

print "DONE H", H;
print "root_tuples", root_tuples, "smooth_split", smooth_split,
      "q_split", q_split, "bad11", bad11,
      "torsion_tests", torsion_tests, "hits", #hits;

quit;
