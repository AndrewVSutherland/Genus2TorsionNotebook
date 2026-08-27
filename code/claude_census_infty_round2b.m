// Round 2b: [4], [2,4], [2,2,4] via the automatic-halving parametrization.
// Core: q_{alpha,beta}(x) = x^2 + ((2ab+1)/b^2) x + a^2/b^2 has
// -theta = ((a + b*theta)/1)^2 in the root field IDENTICALLY, and
// q(0) = (a/b)^2 is a square (norm) -- all halving components trivial.
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
function QAB(al, be)  // over any field
    return (Parent(al)!1)*0 + 1, 0;  // placeholder (unused)
end function;

// [4]: K = Q(sqrt2), alpha = 1 + s*w, beta = 1  (w = sqrt2)
printf "PART A [4]: f = x * Norm_{K/Q}(x^2 + (2*alpha+1)x + alpha^2), K=Q(sqrt2), alpha=1+s*sqrt2\n";
K<w> := QuadraticField(2);
PK<X> := PolynomialRing(K);
for s in [Q!1, Q!2, Q!3, Q!-1, Q!-2, Q!1/2, Q!5, Q!-3] do
    al := 1 + s*w;  be := K!1;
    qK := X^2 + ((2*al*be+1)/be^2)*X + al^2/be^2;
    Q4 := P![ Q!(Norm(Coefficient(qK,0))), 0, 0, 0, 0 ];  // rebuild properly below
    // Norm polynomial: product of qK and its conjugate
    qKc := X^2 + Conjugate(Coefficient(qK,1))*X + Conjugate(Coefficient(qK,0));
    prod := qK*qKc;
    cf := [ Q!c : c in Coefficients(prod) ];   // coefficients are Galois-stable
    Q4 := P!cf;
    if not IsIrreducible(Q4) then printf "  s=%o Q4 reducible\n", s; continue; end if;
    f := x*Evaluate(Q4, x);
    Fiber(f, Sprintf("s=%o", s), "[ 4 ]");
end for;

// [2,4]: f = x * q_{alpha,beta} * q_{gamma,delta} with rational (al,be),(ga,de)
printf "PART B [2,4]: f = x*(x^2+(2a+1)x+a^2)*(x^2+((2g d +1)/d^2)x+g^2/d^2), a=t, (g,d)=(2,3)\n";
for t in [Q!1, Q!2, Q!3, Q!-2, Q!5, Q!1/2, Q!-4, Q!7] do
    al := t;
    q1 := x^2 + (2*al+1)*x + al^2;
    ga := Q!2; de := Q!3;
    q2 := x^2 + ((2*ga*de+1)/de^2)*x + ga^2/de^2;
    f := x*q1*q2;
    Fiber(f, Sprintf("t=%o", t), "[ 2, 4 ]");
end for;

// [2,2,4]: f = x(x+4)(x+9) * q_{alpha,beta}, alpha = t, beta = 1
printf "PART C [2,2,4]: f = x(x+4)(x+9)(x^2+(2t+1)x+t^2)\n";
for t in [Q!1, Q!2, Q!3, Q!-2, Q!5, Q!1/2, Q!-4, Q!7] do
    q1 := x^2 + (2*t+1)*x + t^2;
    f := x*(x+4)*(x+9)*q1;
    Fiber(f, Sprintf("t=%o", t), "[ 2, 2, 4 ]");
end for;
printf "ROUND2B_DONE\n";
quit;
