// claude_222_qf_geometry.m — [2,22] side geometry (Magma):
// (1) genus + points of the Flynn (4,3) quadratic-factor component
//     u^4+(2v-2)u^3+(v^2-6v+3)u^2+(-8v^2+10v-2)u+(-4v^3+8v^2-4v+1) = 0;
// (2) DS v2 (corrected Daowsud–Schmidt) family: rational root over Q(u),
//     quintic cofactor, and its quadratic-factor incidence curve
//     {(u,A,B): x^2+Ax+B | q5_u} — components, genera (fresh [2,22] targets).
SetColumns(0);
SetMemoryLimit(3*10^9);
Q := Rationals();

// ---------- (1) Flynn (4,3) component ----------
A2<u,v> := AffineSpace(Q, 2);
g43 := u^4+(2*v-2)*u^3+(v^2-6*v+3)*u^2+(-8*v^2+10*v-2)*u+(-4*v^3+8*v^2-4*v+1);
C43 := Curve(A2, g43);
printf "Flynn43: irred=%o absirred=%o\n", IsIrreducible(C43), IsAbsolutelyIrreducible(C43);
try
  printf "Flynn43 GENUS %o\n", Genus(C43);
catch e
  printf "Flynn43 genus failed: %o\n", e`Object;
end try;
pts := PointSearch(C43, 10^4);
printf "Flynn43 small points (%o found): %o\n", #pts, [pts[i] : i in [1..Minimum(#pts,8)]];

// ---------- (2) DS v2 ----------
F<t> := FunctionField(Q);
P<x> := PolynomialRing(F);
den := t^5+8;
g := x^6 - (16*t/den)*x^5
   - ((t^15+40*t^10+128*t^5+512)/(2*t^3*den^2))*x^4
   + ((6*t^15+224*t^10+1536*t^5+3072)/(t^2*den^3))*x^3
   + ((t^30+112*t^25+2880*t^20+25600*t^15+106496*t^10+262144*t^5+262144)/(16*t^6*den^4))*x^2
   - ((t^25+80*t^20+1664*t^15+12288*t^10+36864*t^5+32768)/(2*t^5*den^4))*x
   - (t^25+46*t^20+736*t^15+5248*t^10+18432*t^5+24576)/(2*t^4*den^4);
rts := Roots(g);
printf "DSv2 rational roots over Q(t): %o\n", [<r[1],r[2]> : r in rts];
assert #rts ge 1;
r0 := rts[1][1];
q5 := g div (x - r0);
printf "DSv2 quintic cofactor irreducible over Q(t): %o\n", IsIrreducible(q5);

// quadratic-factor incidence of q5: clear denominators, divide by x^2+A*x+B
R3<tt,A,B> := PolynomialRing(Q, 3);
PR<X> := PolynomialRing(R3);
dd := LCM([Denominator(c) : c in Coefficients(q5)]);
q5i := PR![ Evaluate(Numerator(c*dd), tt) : c in Coefficients(q5) ];
qq, rr := Quotrem(q5i, X^2 + A*X + B);
e1 := Coefficient(rr,1); e0 := Coefficient(rr,0);
A3 := AffineSpace(R3);
X3 := Scheme(A3, [e1, e0]);
comps := IrreducibleComponents(X3);
printf "DSv2 quintic QF incidence: %o components\n", #comps;
for i in [1..#comps] do
  Z := ReducedSubscheme(comps[i]);
  d := Dimension(Z);
  printf "-- comp %o: dim %o deg %o\n", i, d, Degree(Z);
  if d eq 1 then
    try
      CC := Curve(Z);
      ab := IsAbsolutelyIrreducible(CC);
      printf "   absirred %o", ab;
      if ab then printf "  GENUS %o", Genus(CC); end if;
      printf "\n";
    catch e
      printf "   genus failed: %o\n", e`Object;
    end try;
  end if;
end for;

// also: root incidence of q5 (extra rational root => [1,1,4]-type member)
R2<t2,rr2> := PolynomialRing(Q, 2);
q5r := &+[ Evaluate(Numerator(Coefficients(q5)[i]*dd), t2) * rr2^(i-1) : i in [1..#Coefficients(q5)] ];
A2b := AffineSpace(R2);
Xr := Scheme(A2b, [q5r]);
compsr := IrreducibleComponents(Xr);
printf "DSv2 quintic extra-root incidence: %o components\n", #compsr;
for i in [1..#compsr] do
  Z := ReducedSubscheme(compsr[i]);
  if Dimension(Z) eq 1 then
    try
      CC := Curve(Z);
      ab := IsAbsolutelyIrreducible(CC);
      printf "-- comp %o: absirred %o", i, ab;
      if ab then printf "  GENUS %o", Genus(CC); end if;
      printf "\n";
      if ab and Genus(CC) le 1 then
        ptsr := PointSearch(CC, 10^4);
        printf "   small points: %o\n", [ptsr[j] : j in [1..Minimum(#ptsr,8)]];
      end if;
    catch e
      printf "-- comp %o genus failed: %o\n", i, e`Object;
    end try;
  end if;
end for;
print "QF_GEOMETRY_DONE";
quit;
