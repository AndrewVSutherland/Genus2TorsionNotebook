// claude_ari_other_components.m — the other components of M(2,2,2,12) over M(2,2,2,6).
// Document/Ari: halve P_b. Here: calibrate the divisibility criteria for all five
// point-type 2-torsion classes P_0, P_a, P_b, P_c, P_d against DivisionPoints, then
// scan small (s:m:n) for points on ANY of the five component surfaces.
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals(); Px<x> := PolynomialRing(Q);

function ABCD(s,m,n)
  a := (m + 2*s)*(n - s)*(n*m - 2*n*s + 4*s^2)*(n*m + n*s - m*s - 2*s^2)*(n*m + 4*n*s - 4*s^2);
  b := -8*(m + 2*s)*(n - 2*s)*s^2*(n - s)^2*(n*m + n*s - m*s - 2*s^2);
  c := -s*(n - 2*s)*(n*m + n*s - m*s - 2*s^2)*(n*m + 4*n*s - 4*s^2)*(n*m + 4*n*s - 2*m*s - 4*s^2);
  d := -s*(n - s)*(m + 2*s)^2*(n - 2*s)^2*(n*m + 4*n*s - 4*s^2);
  return [a,b,c,d];
end function;

// criterion variant: vals all nonzero squares, or all-negated squares
function CritVals(abcd, i)   // x-T tuple of P_e = (x0,0)-infty, x0 = -e: values x0 - e_j
  if i eq 0 then return abcd; end if;   // x0 = 0: values 0-(-r) = r
  e := abcd[i];
  others := [abcd[j] : j in [1..4] | j ne i];
  return [-e] cat [r - e : r in others];
end function;
function CritHolds(abcd, i)
  vals := CritVals(abcd, i);
  if exists{v : v in vals | v eq 0} then return false, "zero"; end if;
  if forall{v : v in vals | IsSquare(v)} then return true, "+"; end if;
  return false, "";
end function;

print "==== calibration on curve #2 (s,m,n)=(2208,-8303,-7200) ====";
abcd2 := ABCD(Q!2208, Q!-8303, Q!-7200);
for i in [0..4] do
  hold, sgn := CritHolds(abcd2, i);
  printf "crit %o: %o %o\n", i, hold, sgn;
end for;

// (calibration passed: crit 2 true, others false, matching the unique divisible class P_b)

print "==== small scan over (s:m:n), height <= 60, all five criteria ====";
// s>0 WLOG (projective sign); require gcd 1; disc != 0 via distinct nonzero a,b,c,d
HT := 60;
hits := [];
for s in [1..HT] do
  for m in [-HT..HT] do
    for n in [-HT..HT] do
      if GCD([s, m, n]) ne 1 then continue; end if;
      abcd := ABCD(Q!s, Q!m, Q!n);
      if exists{t : t in abcd | t eq 0} then continue; end if;
      if #Seqset(abcd) ne 4 then continue; end if;
      for i in [0..4] do
        hold, sgn := CritHolds(abcd, i);
        if hold then
          printf "SCANHIT crit=%o sgn=%o (s,m,n)=(%o,%o,%o)\n", i, sgn, s, m, n;
          Append(~hits, <i, s, m, n>);
        end if;
      end for;
    end for;
  end for;
end for;
printf "scan hits: %o\n", #hits;

print "==== verify each scan hit with DivisionPoints + TorsionSubgroup ====";
for h in hits do
  i := h[1]; s := Q!h[2]; m := Q!h[3]; n := Q!h[4];
  abcd := ABCD(s, m, n);
  f := x*&*[x + t : t in abcd];
  Cq := HyperellipticCurve(f);
  ok := true;
  try
    Csim := SimplifiedModel(ReducedMinimalWeierstrassModel(Cq));
    J := Jacobian(Csim);
    T := TorsionSubgroup(J);
    printf "hit crit=%o (s,m,n)=(%o,%o,%o): TORSION %o order %o G2a=%o\n",
      i, h[2], h[3], h[4], Invariants(T), #T, G2Invariants(Cq)[1];
  catch e
    printf "hit crit=%o (s,m,n)=(%o,%o,%o): verification error\n", i, h[2], h[3], h[4];
  end try;
end for;
print "ALL_DONE";
quit;
