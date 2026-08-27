//////////////////////////////////////////////////////////////////////
//  Halving diagnostics for Elkies' printed N=32 genus-2 example.
//
//  Elkies gives the curve
//
//      y^2 = (15*x-1)*(1056*x^4 + 156183*x^3
//                      + 26297*x^2 + 649*x - 121)
//
//  with non-Weierstrass points (0,11), (-1,1440), and Weierstrass
//  points infinity and (1/15,0).  The difference between either
//  non-Weierstrass point and either Weierstrass point has order 32.
//
//  This script derives and tests the algebraic condition for halving
//  such a divisor.  If P=(r,y0) and infinity is the chosen Weierstrass
//  point on an odd-degree model y^2=F(X), translate r to 0.  A half of
//  P-infinity exists iff there are
//
//      u(X)   = X^2 + a*X + b,
//      ell(X) = m*X^2 + n*X + k,    k = +/- y0,
//
//  such that
//
//      F(X) - ell(X)^2 = lc(F) * X * u(X)^2.                  (*)
//
//  The top two coefficients solve linearly for a,b:
//
//      a = (F4 - m^2)/(2*F5),
//      b = (F3 - 2*m*n - F5*a^2)/(2*F5).
//
//  The remaining X^2 and X coefficients give two plane equations in
//  m,n.  We factor/result them over Q and count solutions mod good
//  primes.  For the finite Weierstrass point w=1/15, we use the
//  birational change X=1/(x-w), Y=y*X^3, so w becomes infinity.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned output_file then
    output_file := "data/elkies32_halving_conditions.txt";
end if;

Q := Rationals();
Z := Integers();
PX<x> := PolynomialRing(Q);

base_f := (15*x - 1)*(1056*x^4 + 156183*x^3 + 26297*x^2 + 649*x - 121);
wfinite := Q!1/15;

function IntContentOfRationalPolynomial(F)
    coeffs := Coefficients(F);
    if #coeffs eq 0 then
        return Z!1;
    end if;
    den := LCM([Denominator(c) : c in coeffs]);
    ints := [Z!(den*c) : c in coeffs];
    cont := GCD([Abs(c) : c in ints | c ne 0] cat [Z!0]);
    if cont eq 0 then
        return Z!1;
    end if;
    return cont/den;
end function;

function ClearDenoms(F)
    P := Parent(F);
    coeffs := Coefficients(F);
    if #coeffs eq 0 then
        return P!0;
    end if;
    den := LCM([Denominator(c) : c in coeffs]);
    G := P!(den*F);
    ints := [Z!c : c in Coefficients(G)];
    cont := GCD([Abs(c) : c in ints | c ne 0] cat [Z!0]);
    if cont ne 0 then
        G := P!(G/cont);
    end if;
    // Normalize a little, for easier comparison in saved output.
    if #Coefficients(G) gt 0 and Coefficients(G)[#Coefficients(G)] lt 0 then
        G := -G;
    end if;
    return G;
end function;

function ShiftToPoint(F, r)
    return Evaluate(F, x + r);
end function;

function MoveFiniteWeierstrassToInfinity(F, w)
    // If x_old = w + 1/X and Y = y*X^3, then
    // Y^2 = X^6 * F(w + 1/X).
    G := &+[Coefficient(F, i)*(w*x + 1)^i*x^(6-i) : i in [0..Degree(F)]];
    G := PX!G;
    assert Coefficient(G, 6) eq 0;
    return G;
end function;

function GoodPrimeList(F, primes)
    D := Z!Numerator(Discriminant(F));
    bad := Set(PrimeDivisors(2*D));
    return [p : p in primes | not p in bad], Sort(SetToSequence(bad));
end function;

procedure PrintFactorization(out, label, F)
    fprintf out, "%o = %o\n", label, F;
    fprintf out, "%o factors:\n", label;
    fac := Factorization(F);
    if #fac eq 0 then
        fprintf out, "  <unit>\n";
    else
        for pair in fac do
            fprintf out, "  exponent %o: %o\n", pair[2], pair[1];
        end for;
    end if;
end procedure;

function CountModPrime(Fshift, y0, sign, p)
    f5 := Coefficient(Fshift, 5);
    if (Z!(Numerator(f5)) mod p) eq 0 or p eq 2 then
        return -1;
    end if;

    Fp := GF(p);
    Rmn<m,n> := PolynomialRing(Fp, 2);
    f1 := Fp!Coefficient(Fshift, 1);
    f2 := Fp!Coefficient(Fshift, 2);
    f3 := Fp!Coefficient(Fshift, 3);
    f4 := Fp!Coefficient(Fshift, 4);
    f5p := Fp!f5;
    kp := Fp!(sign*y0);

    a := (f4 - m^2)/(2*f5p);
    b := (f3 - 2*m*n - f5p*a^2)/(2*f5p);
    K2 := f2 - n^2 - 2*m*kp - 2*f5p*a*b;
    K1 := f1 - 2*n*kp - f5p*b^2;

    els := [Fp!i : i in [0..p-1]];
    count := 0;
    examples := [];
    for mv in els do
        for nv in els do
            if Evaluate(K1, [mv,nv]) eq 0 and Evaluate(K2, [mv,nv]) eq 0 then
                count +:= 1;
                if #examples lt 5 then
                    Append(~examples, <Z!mv, Z!nv>);
                end if;
            end if;
        end for;
    end for;
    return count, examples;
end function;

procedure AnalyzeDivisor(out, label, Fmodel, r, y0)
    Fshift := ShiftToPoint(Fmodel, r);
    f0 := Coefficient(Fshift, 0);
    f1 := Coefficient(Fshift, 1);
    f2 := Coefficient(Fshift, 2);
    f3 := Coefficient(Fshift, 3);
    f4 := Coefficient(Fshift, 4);
    f5 := Coefficient(Fshift, 5);

    fprintf out, "============================================================\n";
    fprintf out, "%o\n", label;
    fprintf out, "F_shift(X) = %o\n", Fshift;
    fprintf out, "F_shift(0) - y0^2 = %o\n", f0 - y0^2;
    fprintf out, "leading coefficient F5 = %o\n", f5;
    fprintf out, "coefficient vector [F0..F5] = %o\n", [Coefficient(Fshift, i) : i in [0..5]];
    fprintf out, "\n";

    Rmn<m,n> := PolynomialRing(Q, 2, "grevlex");
    f1q := Rmn!f1; f2q := Rmn!f2; f3q := Rmn!f3;
    f4q := Rmn!f4; f5q := Rmn!f5;

    for sign in [-1, 1] do
        k := Rmn!(sign*y0);
        a := (f4q - m^2)/(2*f5q);
        b := (f3q - 2*m*n - f5q*a^2)/(2*f5q);
        K2 := ClearDenoms(f2q - n^2 - 2*m*k - 2*f5q*a*b);
        K1 := ClearDenoms(f1q - 2*n*k - f5q*b^2);

        fprintf out, "sign k=%o\n", sign;
        fprintf out, "a(m,n) = %o\n", a;
        fprintf out, "b(m,n) = %o\n", b;
        PrintFactorization(out, "K2", K2);
        PrintFactorization(out, "K1", K1);

        I := ideal<Rmn | K1, K2>;
        fprintf out, "ideal dimension = %o\n", Dimension(I);

        try
            pts := Variety(I);
            fprintf out, "rational affine solutions in (m,n): %o\n", #pts;
            for P in pts do
                mv := P[1];
                nv := P[2];
                av := Evaluate(a, [mv,nv]);
                bv := Evaluate(b, [mv,nv]);
                fprintf out, "  m=%o n=%o a=%o b=%o\n", mv, nv, av, bv;
            end for;
        catch err
            fprintf out, "Variety(I) failed: %o\n", err`Object;
        end try;

        try
            res_m := ClearDenoms(Resultant(K1, K2, n));
            PrintFactorization(out, "Res_n(K1,K2)", res_m);
        catch err
            fprintf out, "Resultant in n failed: %o\n", err`Object;
        end try;
        fprintf out, "\n";
    end for;

    primes := [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73];
    good, bad := GoodPrimeList(Fmodel, primes);
    fprintf out, "bad primes among tested range from 2*disc(Fmodel): %o\n", [p : p in primes | p in bad];
    fprintf out, "good primes tested: %o\n", good;
    for p in good do
        line := Sprintf("mod %o:", p);
        for sign in [-1, 1] do
            count, examples := CountModPrime(Fshift, y0, sign, p);
            line cat:= Sprintf(" sign %o count %o examples %o;", sign, count, examples);
        end for;
        fprintf out, "%o\n", line;
    end for;
    fprintf out, "\n";
end procedure;

out := Open(output_file, "w");

fprintf out, "# Elkies N=32 printed example: halving diagnostics\n\n";
fprintf out, "base_f = %o\n", base_f;
fprintf out, "expanded base_f = %o\n", base_f;
fprintf out, "disc(base_f) = %o\n\n", Discriminant(base_f);

fprintf out, "Generic translated halving equations for F(X)=sum F_i X^i and P=(0,y0):\n";
fprintf out, "  k = +/- y0\n";
fprintf out, "  a = (F4 - m^2)/(2*F5)\n";
fprintf out, "  b = (F3 - 2*m*n - F5*a^2)/(2*F5)\n";
fprintf out, "  K2 = F2 - n^2 - 2*m*k - 2*F5*a*b = 0\n";
fprintf out, "  K1 = F1 - 2*n*k - F5*b^2 = 0\n\n";

AnalyzeDivisor(out, "D00 = (0,11) - infinity", base_f, Q!0, Q!11);
AnalyzeDivisor(out, "Dm10 = (-1,1440) - infinity", base_f, Q!-1, Q!1440);

gfinite := MoveFiniteWeierstrassToInfinity(base_f, wfinite);
fprintf out, "Finite Weierstrass transformation:\n";
fprintf out, "  w = %o\n", wfinite;
fprintf out, "  X = 1/(x-w), Y = y*X^3\n";
fprintf out, "  g(X) = X^6*f(w+1/X) = %o\n\n", gfinite;

R0 := 1/(Q!0 - wfinite);
Y0 := Q!11 * R0^3;
AnalyzeDivisor(out, "D0w = (0,11) - (1/15,0), after w -> infinity", gfinite, R0, Y0);

Rm1 := 1/(Q!-1 - wfinite);
Ym1 := Q!1440 * Rm1^3;
AnalyzeDivisor(out, "Dm1w = (-1,1440) - (1/15,0), after w -> infinity", gfinite, Rm1, Ym1);

delete out;

print "Wrote", output_file;
quit;
