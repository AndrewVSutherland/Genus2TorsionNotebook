// Analyze Phi38: integral primitive model, F_3 points of projective closure,
// genus over Q, small-height rational points.
S<D,E> := PolynomialRing(RationalField(),2);
// input: raw eval-able Phi38 polynomial (repo-relative; run from the repo
// root, or override with magma -b infile:=/path/to/file)
if not assigned infile then infile := "data/claude_prod_04_35_phi38_poly.txt"; end if;
Phi := eval Read(infile);
Phi := S!Phi;
// primitive integral model
Phi := Phi * LCM([Denominator(c) : c in Coefficients(Phi)]);
Phi := Phi / GCD([Numerator(c) : c in Coefficients(Phi)]);
Phi := S!Phi;
printf "deg %o terms %o\n", TotalDegree(Phi), #Terms(Phi);
printf "content check: gcd of coeffs = %o\n", GCD([Integers()!c : c in Coefficients(Phi)]);

// projective closure over F_3 and small primes: point counts
P2<X,Y,Z> := ProjectiveSpace(RationalField(),2);
PhiH := &+[ MonomialCoefficient(Phi, m) * X^Degree(m,1) * Y^Degree(m,2) * Z^(38-TotalDegree(m)) where m := t : t in Monomials(Phi)];
for p in [3,5,7,9] do
  if p eq 9 then Fq := GF(9); else Fq := GF(p); end if;
  P2p := ProjectiveSpace(Fq,2);
  Rp := CoordinateRing(P2p);
  Phip := Rp!0;
  ok := true;
  for t in Terms(PhiH) do Phip +:= Rp!t; end for;
  if Phip eq 0 then printf "p=%o: reduction vanishes identically!\n", p; continue; end if;
  Cp := Curve(P2p, Phip);
  pts := Points(Cp);
  printf "p=%o: #points of projective closure = %o : %o\n", p, #pts, [Coordinates(q) : q in pts];
end for;
quit;
