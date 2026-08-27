// Census F-column upgrades: fiber certificates turning (infty) into infty.
// For each family: find a fiber with EXACT torsion = G, a STRICT fiber
// (irreducible Frobenius charpoly, all 12 root-powers degree 4), and
// >= 2 distinct G2-invariants (nonconstant modulus).
SetColumns(0);
SetSeed(1);
SetMemoryLimit(8*10^9);
Q := Rationals();
P<x> := PolynomialRing(Q);

function StrictPrime(C)
    // returns p, true if some good p < 300 passes the strict root-power test
    fC, hC := HyperellipticPolynomials(C);
    g := 4*fC + hC^2;
    dsc := Integers()!Numerator(Discriminant(g)/16);
    for p in PrimesInInterval(11, 300) do
        if dsc mod p eq 0 then continue; end if;
        Cp := ChangeRing(SimplifiedModel(C), GF(p));
        Lp := LPolynomial(Cp);
        Pp := P!Reverse(Coefficients(Lp));
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

famnames := ["[2]", "[3]", "[2,2]", "[2,2,2]", "[2,2,2,2]", "[2,2,2,4]",
             "[2,2,4,4]", "[2,2,2,8]", "[8]"];
function FamCurve(i, t)
    if i eq 1 then return x*(x^5 + x + t); end if;
    if i eq 2 then return (x^3+x+1)^2 - t*(x^2+1)^3; end if;
    if i eq 3 then return x*(x-1)*(x^4 + x + t); end if;
    if i eq 4 then return x*(x-1)*(x-2)*(x^3 + x + t); end if;
    if i eq 5 then return x*(x-1)*(x-2)*(x-3)*(x-t); end if;
    if i eq 6 then return x*(x+1)*(x+4)*(x+9)*(x+t^2); end if;
    if i eq 7 then
        s := t;
        a := (s-2)*(s^2-s+1/2)*(s^2-3/2)*(s^4+s^2+9/4);
        b := (s^2-2*s+1/2)*(s^2-2*s+9/2)*(s^2-3/2)*(s^2+1);
        c := -(s-2)*(2*s^2-2*s+1)*(s^4+s^2+9/4);
        d := (s-2)*(2*s^2-2*s+1)*(2*s^2-3)*(s^2+1);
        return x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2);
    end if;
    if i eq 8 then
        a := -4*t^2*(t+1)/(t^2+t+1)^2;
        b := -t/(t+1); c := Q!1; d := t;
        return x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2);
    end if;
    if i eq 9 then
        // A(8) chart, slice p=2, r=3, parameter t (verbatim formulas)
        r := Q!3; p := Q!2;
        e := t^2 - 2*p*t/r;
        d := e + 2*p - r^2;
        lambda := r/t;
        aa := r^2 - lambda;
        bb := 2*r*p - 2*lambda*(p + r*t) + 2*r*lambda;
        cc := p^2 + 2*p*r^2 - r^4 - r^3*t - r*p^2/t
            - lambda*(r^2 + e) + 2*lambda*(r*p + r^2*t - 3*p*t + r*t^2);
        Qp := x^2 + d;
        q := aa*x^2 + bb*x + cc;
        return q*(Qp^2 + q);
    end if;
    return P!0;
end function;

tvals := [ Q!2, Q!3, Q!5, Q!7, Q!-2, Q!-3, Q!1/2, Q!-1/3, Q!11, Q!2/5 ];
for i in [1..#famnames] do
    G := famnames[i];
    exactT := 0; strictT := 0; strictP := 0;
    invs := {};
    printf "FAMILY %o %o\n", i, G;
    for t in tvals do
        okf := true;
        try
            f := FamCurve(i, t);
        catch e okf := false; end try;
        if not okf or Degree(f) lt 5 or not IsSquarefree(f) then continue; end if;
        den := LCM([Denominator(c) : c in Coefficients(f)]);
        f := den^2 * f;   // (y -> den*y): same curve, integral model
        C := HyperellipticCurve(f);
        Include(~invs, G2Invariants(C));
        T := TorsionSubgroup(Jacobian(SimplifiedModel(C)));
        tstr := "[" cat (#Invariants(T) eq 0 select "" else
            &cat[ IntegerToString(Invariants(T)[k]) cat (k lt #Invariants(T) select "," else "") : k in [1..#Invariants(T)] ]) cat "]";
        printf "  t=%o torsion %o\n", t, tstr;
        if tstr eq G and Type(exactT) eq RngIntElt then
            exactT := t;
        end if;
        if tstr eq G and strictP eq 0 then
            sp, oks := StrictPrime(C);
            if oks then strictT := t; strictP := sp; end if;
        end if;
        if Type(exactT) ne RngIntElt and strictP ne 0 and #invs ge 2 then break; end if;
    end for;
    printf "CERT %o exact_fiber=%o strict_fiber=%o strict_prime=%o distinct_invs=%o\n",
        G, exactT, strictT, strictP, #invs;
end for;
printf "INFTY_CERTIFY_DONE\n";
quit;
