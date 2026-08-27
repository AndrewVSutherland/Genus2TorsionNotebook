// Round 2: [2,4] derived family; [2,2,4] and [4,4] hits-first design.
SetColumns(0);
SetSeed(1);
SetMemoryLimit(8*10^9);
Q := Rationals();
P<x> := PolynomialRing(Q);
function StrictPrime(C)
    fC, hC := HyperellipticPolynomials(C);
    g := 4*fC + hC^2;
    dsc := Integers()!Numerator(Discriminant(g)/16);
    for p in PrimesInInterval(11, 600) do
        if dsc mod p eq 0 then continue; end if;
        Pp := P!Reverse(Coefficients(LPolynomial(ChangeRing(SimplifiedModel(C), GF(p)))));
        if not IsIrreducible(Pp) then continue; end if;
        K<pi> := NumberField(Pp);
        ok := true;
        for n in [1..12] do
            if Degree(MinimalPolynomial(pi^n)) ne 4 then ok := false; break; end if;
        end for;
        if ok then return p, true; end if;
    end for;
    return 0, false;
end function;
procedure TryFiber(f, tag, target)
    if Degree(f) lt 5 or not IsSquarefree(f) then return; end if;
    den := LCM([Denominator(c) : c in Coefficients(f)]);
    C := HyperellipticCurve(den^2*f);
    T := TorsionSubgroup(Jacobian(SimplifiedModel(C)));
    inv := Sprint(Invariants(T));
    printf "  %o torsion %o", tag, inv;
    if inv eq target then
        sp, oks := StrictPrime(C);
        printf "  STRICT p=%o(%o)", sp, oks;
    end if;
    printf "\n";
end procedure;

// A) [2,4]: f = x(x^2 + 2pq x + u^2)(x^2 - (p^2-u^2)/(1-q^2))
//    (halving conditions for the q2-class solved rationally)
printf "PART A [2,4] family slice p=t, q=2, u=1: f = x(x^2+4tx+1)(x^2-(1-t^2)/3)\n";
for t in [Q!2, Q!3, Q!5, Q!7, Q!-2, Q!1/2] do
    f := x*(x^2 + 4*t*x + 1)*(x^2 - (1-t^2)/3);
    TryFiber(f, Sprintf("t=%o", t), "[ 2, 4 ]");
end for;

// B) [2,2,4]: mini-sweep of x(x-1)(x-L)(x^2+ax+b) for exact [2,2,4]
printf "PART B [2,2,4] hits-first sweep\n";
nb := 0;
for L in [Q!2,Q!3,Q!4,Q!5,Q!-1,Q!-2,Q!1/2,Q!3/2], a in [-4..4], b in [-6..6] do
    if b eq 0 then continue; end if;
    f := x*(x-1)*(x-L)*(x^2+a*x+b);
    if not IsSquarefree(f) then continue; end if;
    // quick prefilter: gcd of #J(F_p) divisible by 16
    ok := true; g := 0;
    for p in [7,11,13] do
        dn := LCM([Denominator(c) : c in Coefficients(f)]);
        fi := dn^2*f;
        if Integers()!Numerator(Discriminant(fi)) mod p eq 0 then continue; end if;
        np := #Jacobian(ChangeRing(HyperellipticCurve(fi), GF(p)));
        g := GCD(g, np);
        if g ne 0 and g mod 16 ne 0 then ok := false; break; end if;
    end for;
    if not ok or g eq 0 then continue; end if;
    den := LCM([Denominator(c) : c in Coefficients(f)]);
    C := HyperellipticCurve(den^2*f);
    T := TorsionSubgroup(Jacobian(SimplifiedModel(C)));
    if Sprint(Invariants(T)) eq "[ 2, 2, 4 ]" then
        nb +:= 1;
        printf "  HIT224 L=%o a=%o b=%o\n", L, a, b;
        if nb ge 8 then break L; end if;
    end if;
end for;

// C) [4,4]: search (beta,q,u) making the final conic square, then fibers
printf "PART C [4,4] condition search\n";
nc := 0;
for be in [Q!2,Q!3,Q!5,Q!2/3,Q!3/2,Q!5/2,Q!4,Q!5/3], qq in [Q!2,Q!3,Q!1/2,Q!5,Q!2/5,Q!4], u in [Q!1,Q!2,Q!3] do
    // alpha^2 * [ be^2/(qq^2(be^2-1)^2(1-qq^2)) + 1 ] = u^2/(1-qq^2) - u^2(1-be^2)
    lhsc := be^2/(qq^2*(be^2-1)^2*(1-qq^2)) + 1;
    rhs := u^2/(1-qq^2) - u^2*(1-be^2);
    if lhsc eq 0 then continue; end if;
    v := rhs/lhsc;
    sq, al := IsSquare(v);
    if not sq or al eq 0 then continue; end if;
    a := 2*al*be/(be^2-1);
    p := al*be/(qq*(be^2-1));
    b := u^2;
    e := (p^2 - u^2)/(1 - qq^2);
    f := x*(x^2 + a*x + b)*(x^2 - e);
    if not IsSquarefree(f) or Degree(f) ne 5 then continue; end if;
    den := LCM([Denominator(c) : c in Coefficients(f)]);
    C := HyperellipticCurve(den^2*f);
    T := TorsionSubgroup(Jacobian(SimplifiedModel(C)));
    printf "  C-cand be=%o q=%o u=%o torsion %o\n", be, qq, u, Invariants(T);
    nc +:= 1;
    if nc ge 12 then break be; end if;
end for;
printf "ROUND2_DONE\n";
quit;
