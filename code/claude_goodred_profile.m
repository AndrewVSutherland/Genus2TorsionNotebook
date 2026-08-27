// Item 1, Step 3: chart-free good-reduction profile.
// For p in {5,7,...,47}: enumerate ALL genus 2 curves /F_p with 6 rational
// Weierstrass points (up to PGL2: 6-subsets of P^1(F_p) containing {0,1,oo},
// both quadratic twists), compute J(F_p), and test containment of
// (2,2,4,8), (2,2,2,12), and context groups (2,2,4,4), (2,2,2,8).
// Torsion injects into J(F_p) at odd good p, so "no curve at p contains G"
// proves: every genus-2 Jacobian /Q with torsion >= G has bad reduction at p.

// does the abelian group with invariant factor list inv contain grp t?
contains := function(inv, t)
    // p-part domination test for each prime dividing #t
    for pr in {x[1] : x in Factorization(&*t)} do
        te := Reverse(Sort([Valuation(m, pr) : m in t | Valuation(m, pr) gt 0]));
        ie := Reverse(Sort([Valuation(m, pr) : m in inv | Valuation(m, pr) gt 0]));
        if #te gt #ie then return false; end if;
        for i in [1..#te] do
            if te[i] gt ie[i] then return false; end if;
        end for;
    end for;
    return true;
end function;

targets := [ <"(2,2,4,8)",  [2,2,4,8]>,
             <"(2,2,2,12)", [2,2,2,12]>,
             <"(2,2,4,4)",  [2,2,4,4]>,
             <"(2,2,2,8)",  [2,2,2,8]> ];

for p in [5,7,11,13,17,19,23,29,31,37,41,43,47] do
    Fp := GF(p);
    ns := [x : x in Fp | not IsSquare(x) and x ne 0][1];
    // 6-subsets of P^1(F_p) containing {0,1,oo}: choose 3 more finite pts != 0,1
    rest := [x : x in Fp | x ne 0 and x ne 1];
    counts := AssociativeArray();
    for t in targets do counts[t[1]] := 0; end for;
    frobs := AssociativeArray();
    for t in targets do frobs[t[1]] := {}; end for;
    total := 0;
    P<x> := PolynomialRing(Fp);
    for i in [1..#rest] do
    for j in [i+1..#rest] do
    for k in [j+1..#rest] do
        // Weierstrass pts {0,1,rest[i],rest[j],rest[k],oo}: quintic model
        f0 := x*(x-1)*(x-rest[i])*(x-rest[j])*(x-rest[k]);
        for c in [Fp!1, ns] do
            f := c*f0;
            C := HyperellipticCurve(f);
            J := Jacobian(C);
            A := AbelianGroup(J);
            inv := Invariants(A);
            total +:= 1;
            for t in targets do
                if contains(inv, t[2]) then
                    counts[t[1]] +:= 1;
                    Include(~frobs[t[1]], EulerFactor(J));
                end if;
            end for;
        end for;
    end for;
    end for;
    end for;
    printf "p=%o curves=%o :", p, total;
    for t in targets do
        printf "  %o:%o", t[1], counts[t[1]];
    end for;
    printf "\n";
    // print the surviving Frobenius L-polynomials for the two hard targets at small p
    if p le 13 then
        for t in [targets[1], targets[2]] do
            if counts[t[1]] gt 0 then
                printf "   %o L-polys at p=%o: %o\n", t[1], p, frobs[t[1]];
            end if;
        end for;
    end if;
end for;
print "PROFILE COMPLETE";
quit;
