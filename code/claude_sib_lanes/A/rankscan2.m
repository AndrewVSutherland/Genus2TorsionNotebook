P<u> := PolynomialRing(Rationals());
q := 4*u^2 - 6*u + 3;
for b in [17..32] do
 for a in [-80..80] do
  if a ne 0 and AbsoluteValue(a) ge b and GCD(AbsoluteValue(a), b) eq 1 then
    g := a/b;
    f2 := q*((g^2+6*g+1) - 8*g*u);
    ok, E2 := IsEllipticCurve(HyperellipticCurve(f2));
    if ok then
      lo, hi := RankBounds(E2 : Effort := 1);
      if lo ge 5 then printf "g=%o/%o: rank [%o,%o]\n", a, b, lo, hi; end if;
    end if;
  end if;
 end for;
 printf "b=%o done\n", b;
end for;
quit;
