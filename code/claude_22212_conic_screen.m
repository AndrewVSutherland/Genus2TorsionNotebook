// Screen the deduped conic-planes on S for nondegeneracy + rational points.
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals();
load "data/claude_22212_conic_planes.m";
PU<u,v,w> := PolynomialRing(Q, 3);
nsurv := 0;
for ip in [1..#PLANES] do
    M := PLANES[ip];
    X := [ M[1][k]*u + M[2][k]*v + M[3][k]*w : k in [1..5] ];
    Q2 := X[1]^2+X[2]^2+X[3]^2-X[4]^2-X[5]^2;
    if Q2 eq 0 then continue; end if;
    // degeneracy forms restricted to the plane
    forms := [ X[i] : i in [1..5] ];
    for i in [1..4] do for j in [i+1..5] do
        Append(~forms, X[i]-X[j]); Append(~forms, X[i]+X[j]);
    end for; end for;
    bad := false;
    for L in forms do
        if L eq 0 then bad := true; break; end if;                 // form vanishes on plane
        if IsDivisibleBy(Q2, L) then bad := true; break; end if;   // form divides conic
    end for;
    if bad then continue; end if;
    // conic must be geometrically irreducible: not a product of two lines
    fc := Factorization(Q2);
    if #fc gt 1 or fc[1][2] gt 1 or TotalDegree(fc[1][1]) lt 2 then continue; end if;
    nsurv +:= 1;
    // rational point?
    Pr2 := ProjectiveSpace(Q, 2);
    Cn := Conic(Pr2, Evaluate(Q2, [Pr2.1, Pr2.2, Pr2.3]));
    hp, pt := HasRationalPoint(Cn);
    printf "SURVIVOR %o: plane %o | conic %o | ratpt %o", ip, M, Q2, hp;
    if hp then printf " %o", pt; end if;
    printf "\n";
end for;
printf "survivors (nondegenerate irreducible conics): %o\n", nsurv;
printf "SCREEN_DONE\n";
quit;
