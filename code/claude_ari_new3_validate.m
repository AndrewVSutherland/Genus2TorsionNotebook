// claude_ari_new3_validate.m — full validation of NEWCURVE_3 from surface point (5364,19661,4165)
// per validate-and-record-a-hit: fresh rebuild, exact TorsionSubgroup, simplicity certificates.
SetColumns(0);
SetMemoryLimit(10*10^9);
Q := Rationals(); Px<x> := PolynomialRing(Q);

Av := Q!5364; Bv := Q!19661; Cv := Q!4165;
// verify surface membership from scratch
Fv := (Bv^2-Av^2)*(Bv^2-Cv^2)*(Bv^2-Av^2-Cv^2);
Gv := (Bv^2-Av^2-Cv^2)*(Bv^2*(Av^2+Cv^2) - (Av^4+Av^2*Cv^2+Cv^4));
okF, yv := IsSquare(Fv); okG, zv := IsSquare(Gv);
printf "surface membership: F square %o (y=%o), G square %o (z=%o)\n", okF, yv, okG, zv;
assert okF and okG;

s := (Av^3 - Av*Bv^2 + 2*Av*Cv^2)/(2*Av^2*Bv - 2*Bv^3 + 2*Bv*Cv^2);
m := (-2*Av^2*Bv^2 + 2*Bv^4 + 2*Av^2*Cv^2 - 6*Bv^2*Cv^2 + 4*Cv^4)/(2*Av^3*Bv - 2*Av*Bv^3 + 2*Av*Bv*Cv^2);
n := Av/Bv;
den := LCM([Denominator(s), Denominator(m), Denominator(n)]);
si := Integers()!(s*den); mi := Integers()!(m*den); ni := Integers()!(n*den);
g := GCD([si, mi, ni]); si := si div g; mi := mi div g; ni := ni div g;
printf "(s:m:n) = (%o : %o : %o)\n", si, mi, ni;
s := Q!si; m := Q!mi; n := Q!ni;

a := (m + 2*s)*(n - s)*(n*m - 2*n*s + 4*s^2)*(n*m + n*s - m*s - 2*s^2)*(n*m + 4*n*s - 4*s^2);
b := -8*(m + 2*s)*(n - 2*s)*s^2*(n - s)^2*(n*m + n*s - m*s - 2*s^2);
c := -s*(n - 2*s)*(n*m + n*s - m*s - 2*s^2)*(n*m + 4*n*s - 4*s^2)*(n*m + 4*n*s - 2*m*s - 4*s^2);
d := -s*(n - s)*(m + 2*s)^2*(n - 2*s)^2*(n*m + 4*n*s - 4*s^2);
fq := x*(x+a)*(x+b)*(x+c)*(x+d);
assert Discriminant(fq) ne 0;
printf "quintic built; coefficient size ~ %o digits\n", #IntegerToString(Integers()!Max([Abs(Numerator(cc)) : cc in Coefficients(fq)]));

// clear squares from f by scaling x -> u^2*x, y -> u^5*y won't shrink content; use minimization
Cq := HyperellipticCurve(fq);
tmin := Cputime();
Cmin := ReducedMinimalWeierstrassModel(Cq);
printf "minimization time %o s\n", Cputime(tmin);
fmin, hmin := HyperellipticPolynomials(Cmin);
printf "MINMODEL h = %o\n", hmin;
printf "MINMODEL f = %o\n", fmin;

Csim := SimplifiedModel(Cmin);
fsim, hsim := HyperellipticPolynomials(Csim);
assert hsim eq 0 and forall{cc : cc in Coefficients(fsim) | Denominator(cc) eq 1};
printf "simplified model y^2 = %o\n", fsim;
J := Jacobian(Csim);
tt := Cputime();
Tor := TorsionSubgroup(J);
printf "TORSION: %o (order %o)  [%o s]\n", Invariants(Tor), #Tor, Cputime(tt);

D := Integers()!Discriminant(Cmin);
nc := 0;
for p in PrimesInInterval(29, 400) do
  if nc ge 4 then break; end if;
  if D mod p eq 0 then continue; end if;
  chi := Px!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(Cmin, GF(p))))));
  if Degree(chi) ne 4 or not IsIrreducible(chi) then continue; end if;
  K<aK> := NumberField(chi); c12 := MinimalPolynomial(aK^12);
  if IsIrreducible(c12) and Degree(c12) eq 4 then
    printf "p=%o: simplicity cert OK (chi and chi^(12) irreducible deg 4)\n", p; nc +:= 1;
  end if;
end for;
printf "certificates: %o\n", nc;
printf "G2Invariants: %o\n", G2Invariants(Cmin);
if Invariants(Tor) eq [2,2,2,12] and nc ge 2 then
  printf "CONFIRMED: THIRD geometrically simple (2,2,2,12) curve over Q\n";
end if;
print "ALL_DONE";
quit;
