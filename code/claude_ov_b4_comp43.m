// claude_ov_b4_comp43.m -- Lane 4 (route B4): the SECOND component of the
// quadratic-factor incidence of Flynn's order-11 family.
//
// Res_t(e1,e0) factors over Q[u,v] into exactly TWO irreducible components:
//   (6,2):  (u-2)^4 v^2 + 2(u^5-5u^4+9u^3-8u^2+4u-2) v + u^2(u^2-u+1)^2   [known]
//   (4,3): -4v^3 + (u^2-8u+8) v^2 + (2u^3-6u^2+10u-4) v + (u^2-u+1)^2     [NEW]
// The whole overnight Lane-4 campaign used only the (6,2) stream.  This script
// decides the (4,3) component: geometric genus, a smooth model, rank if
// elliptic, and (if it has points) the seed members it produces.
//
// Run: nohup magma -b code/claude_ov_b4_comp43.m > results/claude_ov_b4_comp43.log 2>&1 &

SetColumns(0);
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
P2<u,v> := PolynomialRing(Q, 2);
C62 := (u-2)^4*v^2 + 2*(u^5-5*u^4+9*u^3-8*u^2+4*u-2)*v + u^2*(u^2-u+1)^2;
C43 := -4*v^3 + (u^2-8*u+8)*v^2 + (2*u^3-6*u^2+10*u-4)*v + (u^2-u+1)^2;

// --- independent re-derivation of the two components ---------------------
Pt<t> := PolynomialRing(P2);
Px<xx> := PolynomialRing(Pt);
Ft := xx^6 + 2*xx^5 + (2*t+3)*xx^4 + 2*xx^3 + (t^2+1)*xx^2 + 2*t*(1-t)*xx + t^2;
_, rem := Quotrem(Ft, xx^2 + (Pt!u)*xx + (Pt!v));
e1 := Coefficient(rem,1); e0 := Coefficient(rem,0);
R := Resultant(e1, e0);
printf "RESULTANT factorization: %o\n", [<Degree(f[1],u), Degree(f[1],v), f[2]> : f in Factorization(R)];
printf "C62 divides R: %o\n", IsDivisibleBy(R, C62);
printf "C43 divides R: %o\n", IsDivisibleBy(R, C43);

A2 := AffineSpace(P2);
for pair in [<"(6,2)", C62>, <"(4,3)", C43>] do
    nm := pair[1]; eqn := pair[2];
    Cc := Curve(A2, eqn);
    printf "\n=== COMPONENT %o ===\n", nm;
    printf "  irreducible: %o  arith/geom genus: ", IsIrreducible(eqn);
    g := Genus(ProjectiveClosure(Cc));
    printf "%o\n", g;
    if g eq 1 then
        b, EE := IsHyperelliptic(ProjectiveClosure(Cc));
        pts := PointSearch(ProjectiveClosure(Cc), 500);
        printf "  PointSearch(500): %o points, e.g. %o\n", #pts, [pts[i] : i in [1..Min(4,#pts)]];
        if #pts gt 0 then
            EJ, mp := EllipticCurve(ProjectiveClosure(Cc), pts[1]);
            EJ := MinimalModel(EJ);
            printf "  ELLIPTIC MODEL: %o  cond=%o\n", aInvariants(EJ), Conductor(EJ);
            rb := RankBounds(EJ);
            printf "  RANK BOUNDS: %o\n", [rb];
            printf "  TORSION: %o\n", Invariants(TorsionSubgroup(EJ));
            if Rank(EJ) gt 0 then
                printf "  GENERATORS: %o\n", Generators(EJ);
            end if;
        end if;
    elif g eq 0 then
        printf "  GENUS 0: rational or a conic; PointSearch: ";
        pts := PointSearch(ProjectiveClosure(Cc), 500);
        printf "%o points\n", #pts;
    else
        pts := PointSearch(ProjectiveClosure(Cc), 2000);
        printf "  GENUS %o (Faltings-finite).  PointSearch(2000): %o points\n", g, #pts;
        for pp in pts do printf "    PT %o\n", pp; end for;
    end if;
end for;

printf "COMP43_DONE\n";
quit;
