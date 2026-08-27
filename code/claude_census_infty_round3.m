// Round 3: [10] (contact-5 generic), [2,8] (A(8) with split q),
// [2,10] (contact-5 with [2,2]-split residual quartic).
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
procedure Fiber(f, tag, target)
    if Degree(f) lt 5 or not IsSquarefree(f) then printf "  %o degenerate\n", tag; return; end if;
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

// A) [10]: contact-5 generic slice b = 1: f = (1+ax+x^2)^2 - (2+a)^2 x^5
printf "PART A [10]: f = (1+ax+x^2)^2 - (a+2)^2 x^5\n";
for a in [Q!1, Q!2, Q!3, Q!-3, Q!5, Q!1/2, Q!7, Q!-5] do
    f := (1+a*x+x^2)^2 - (a+2)^2*x^5;
    Fiber(f, Sprintf("a=%o", a), "[ 10 ]");
end for;

// B) [2,8]: A(8) chart with q-disc square.  q = a x^2 + b x + c;
// scan the (r,p,t)-chart for disc(q) square, report torsion of those fibers.
printf "PART B [2,8]: A(8) fibers with split q\n";
nb := 0;
for r in [Q!2,Q!3,Q!5,Q!3/2], p in [Q!1,Q!2,Q!-1,Q!1/2,Q!3], t in [Q!1,Q!2,Q!-1,Q!-2,Q!1/2,Q!3,Q!5,Q!-3] do
    if r eq 0 or t eq 0 then continue; end if;
    e := t^2 - 2*p*t/r;
    d := e + 2*p - r^2;
    lambda := r/t;
    aa := r^2 - lambda;
    bb := 2*r*p - 2*lambda*(p + r*t) + 2*r*lambda;
    cc := p^2 + 2*p*r^2 - r^4 - r^3*t - r*p^2/t
        - lambda*(r^2 + e) + 2*lambda*(r*p + r^2*t - 3*p*t + r*t^2);
    if aa eq 0 then continue; end if;
    disc := bb^2 - 4*aa*cc;
    if not IsSquare(disc) then continue; end if;
    Qp := x^2 + d;
    q := aa*x^2 + bb*x + cc;
    f := q*(Qp^2 + q);
    if Degree(f) ne 6 or not IsSquarefree(f) then continue; end if;
    nb +:= 1;
    Fiber(f, Sprintf("(r,p,t)=(%o,%o,%o)", r, p, t), "[ 2, 8 ]");
    if nb ge 10 then break r; end if;
end for;

// C) [2,10]: contact-5 with residual quartic split into two rational quadratics
printf "PART C [2,10]: contact-5, quartic = q2*q2' locus (search)\n";
nc := 0;
for a in [-8..8], b in [-8..8] do
    if 1+a+b eq 0 then continue; end if;
    f := (1+a*x+b*x^2)^2 - (1+a+b)^2*x^5;
    if Degree(f) lt 5 or not IsSquarefree(f) then continue; end if;
    quart := f div (x-1);
    if Degree(quart) ne 4 then continue; end if;
    fac := Factorization(quart);
    degs := Sort([Degree(t[1]) : t in fac]);
    if degs eq [2,2] then
        nc +:= 1;
        Fiber(f, Sprintf("(a,b)=(%o,%o)", a, b), "[ 2, 10 ]");
        if nc ge 10 then break a; end if;
    end if;
end for;
printf "ROUND3_DONE\n";
quit;
