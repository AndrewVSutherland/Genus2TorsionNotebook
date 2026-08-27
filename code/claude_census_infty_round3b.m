// Locus parametrizations for [2,10] and [2,8]; probe [2,2,10].
SetColumns(0);
SetSeed(1);
SetMemoryLimit(8*10^9);
Q := Rationals();
P<x> := PolynomialRing(Q);

// A) [2,10]: contact-5 quartic-splitting locus.
// f_{a,b} = (1+ax+bx^2)^2-(1+a+b)^2 x^5, quart = f/(x-1).
// Splitting into two rational quadratics <=> resolvent cubic has a rational
// root z AND the associated discriminants are squares.  Explore the surface
// R(a,b,z)=0: for fixed z solve for the (a,b)-curve; test rationality by
// finding many points and attempting a parametrization via Conic/ratl curve.
QF<aa,bb> := FunctionField(Q, 2);
PF<X> := PolynomialRing(QF);
fF := (1+aa*X+bb*X^2)^2 - (1+aa+bb)^2*X^5;
quart := fF div (X-1);
printf "A) quartic coefficients (deg in b):\n";
cf := Coefficients(quart);
for i in [1..#cf] do printf "  x^%o: %o\n", i-1, cf[i]; end for;
// resolvent cubic of monic-normalized quartic q(x) = x^4+c3 x^3+c2 x^2+c1 x+c0:
lc := cf[5];
c3 := cf[4]/lc; c2 := cf[3]/lc; c1 := cf[2]/lc; c0 := cf[1]/lc;
PZ<z> := PolynomialRing(QF);
resv := z^3 - c2*z^2 + (c1*c3 - 4*c0)*z - (c1^2 + c0*c3^2 - 4*c0*c2);
printf "A) resolvent degrees in (a,b): %o\n", [Degree(Numerator(co), i) : i in [1,2], co in Coefficients(resv)];

// B) [2,8]: A(8) disc-square locus on the line (r,p) = (3/2, 3):
// disc(t) = bb^2 - 4 aa cc as a rational function; y^2 = disc(t): genus?
QT<tt> := FunctionField(Q);
r := QT!3/2; p := QT!3;
e := tt^2 - 2*p*tt/r;
lambda := r/tt;
A8a := r^2 - lambda;
A8b := 2*r*p - 2*lambda*(p + r*tt) + 2*r*lambda;
A8c := p^2 + 2*p*r^2 - r^4 - r^3*tt - r*p^2/tt
    - lambda*(r^2 + e) + 2*lambda*(r*p + r^2*tt - 3*p*tt + r*tt^2);
disc := A8b^2 - 4*A8a*A8c;
num := Numerator(disc); den := Denominator(disc);
printf "B) disc numerator (in t): %o\n", num;
printf "B) disc denominator: %o\n", den;
// square condition y^2 = num*den (mod squares)
PT<T> := PolynomialRing(Q);
g := PT!(num*den);
g := g div GCD(g, Derivative(g))^0;   // keep as is
sqfree := &*[ t[1]^(t[2] mod 2) : t in Factorization(g) ];
printf "B) squarefree core: %o (degree %o)\n", sqfree, Degree(sqfree);
if Degree(sqfree) le 2 then
    printf "B) => genus 0 conic-type condition!\n";
elif Degree(sqfree) in [3,4] then
    E := HyperellipticCurve(sqfree);
    printf "B) => genus 1; curve y^2 = %o\n", sqfree;
end if;

// C) [2,2,10] probe: contact-5 with quartic of factor type [1,1,2]
printf "C) [2,2,10] probe\n";
function StrictPrime(C)
    fC, hC := HyperellipticPolynomials(C);
    gg := 4*fC + hC^2;
    dsc := Integers()!Numerator(Discriminant(gg)/16);
    for pq in PrimesInInterval(11, 600) do
        if dsc mod pq eq 0 then continue; end if;
        Pp := P!Reverse(Coefficients(LPolynomial(ChangeRing(SimplifiedModel(C), GF(pq)))));
        if not IsIrreducible(Pp) then continue; end if;
        K<pi> := NumberField(Pp);
        ok := true;
        for n in [1..12] do
            if Degree(MinimalPolynomial(pi^n)) ne 4 then ok := false; break; end if;
        end for;
        if ok then return pq, true; end if;
    end for;
    return 0, false;
end function;
nc := 0;
for a in [-10..10], b in [-10..10] do
    if 1+a+b eq 0 then continue; end if;
    f := (1+a*x+b*x^2)^2 - (1+a+b)^2*x^5;
    if Degree(f) lt 5 or not IsSquarefree(f) then continue; end if;
    quart2 := f div (x-1);
    if Degree(quart2) ne 4 then continue; end if;
    degs := Sort([Degree(t[1]) : t in Factorization(quart2)]);
    if degs eq [1,1,2] then
        C := HyperellipticCurve(f);
        T := TorsionSubgroup(Jacobian(SimplifiedModel(C)));
        inv := Sprint(Invariants(T));
        printf "  (a,b)=(%o,%o) torsion %o", a, b, inv;
        if inv eq "[ 2, 2, 10 ]" then
            sp, oks := StrictPrime(C);
            printf "  STRICT p=%o(%o)", sp, oks;
        end if;
        printf "\n";
        nc +:= 1;
        if nc ge 10 then break a; end if;
    end if;
end for;
printf "ROUND3B_DONE\n";
quit;
