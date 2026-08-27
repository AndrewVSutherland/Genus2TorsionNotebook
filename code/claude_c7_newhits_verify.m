SetColumns(0); SetMemoryLimit(3*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q); Z := Integers();
G := func<u | -(u^5 - u^3 - u^2/2)/(u+1)^2>;
trips := [
  [-511/61, -511/625, -1/2],
  [-165/41, -33/16, -165/289],
  [-164/297, -1/2, 164/361],
  [-17/50, -34/189, 34/121],
  [-1/2, -13/49, 13/50],
  [-3,-3/4,-3/5],
  [-10,-10/7,-1/2]
];
for T in trips do
  v1 := T[1]; v2 := T[2];
  c4 := (G(v1)-G(v2))/(v1^2-v2^2); c0 := G(v1) - c4*v1^2;
  b := c4 - 2; a := 9/2 - c0 - c4;
  h := 1 - (7/2)*x + a*x^2 + b*x^3;
  num := h^2 + (x-1)^7;  f := num div x^2;
  if num ne f*x^2 or Degree(f) ne 5 or Discriminant(f) eq 0 then printf "%o DEGENERATE\n",T; continue; end if;
  // integral model: x -> X/m, y -> Y/m^3
  m := LCM([Denominator(c) : c in Coefficients(f)]);
  F := P![ Coefficient(f,i)*m^(6-i) : i in [0..5] ];
  F := P![ Z!c : c in Coefficients(F) ];
  C := HyperellipticCurve(F);  J := Jacobian(C);
  tors := Invariants(TorsionSubgroup(J));
  Pt := C ! [m*1, m^3*Evaluate(h,1)];
  D := J ! (Pt - PointsAtInfinity(C)[1]);
  ordD := Order(D);
  ft := [Degree(t[1]) : t in Factorization(F)];
  D0 := Z!Discriminant(C); sig := {};
  for p in PrimesInInterval(11,140) do
    if D0 mod p ne 0 then
      chi := P!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(C,GF(p))))));
      d := Coefficient(chi,3)^2 - 4*(Coefficient(chi,2) - 2*p);
      if d ne 0 then Include(~sig, Squarefree(Z!Abs(d))); end if;
    end if;
  end for;
  ss := Sort([s : s in sig]);
  printf "v=%o\n  ftype=%o torsion=%o markedord=%o  #discsig=%o sig=%o\n", T, ft, tors, ordD, #ss, ss[1..Min(6,#ss)];
end for;
printf "NEWHITS_DONE\n"; quit;
