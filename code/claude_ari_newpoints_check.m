// claude_ari_newpoints_check.m — identify + validate curves from the 7 primitive
// height<=500 points on Ari's surface. Which are new (2,2,2,12) curves?
SetColumns(0);
SetMemoryLimit(6*10^9);
Q := Rationals(); Px<x> := PolynomialRing(Q);

function EmailSMN(Av, Bv, Cv)
  s := (Av^3 - Av*Bv^2 + 2*Av*Cv^2)/(2*Av^2*Bv - 2*Bv^3 + 2*Bv*Cv^2);
  m := (-2*Av^2*Bv^2 + 2*Bv^4 + 2*Av^2*Cv^2 - 6*Bv^2*Cv^2 + 4*Cv^4)/(2*Av^3*Bv - 2*Av*Bv^3 + 2*Av*Bv*Cv^2);
  n := Av/Bv;
  return s, m, n;
end function;

function QuinticModel(s, m, n)
  a := (m + 2*s)*(n - s)*(n*m - 2*n*s + 4*s^2)*(n*m + n*s - m*s - 2*s^2)*(n*m + 4*n*s - 4*s^2);
  b := -8*(m + 2*s)*(n - 2*s)*s^2*(n - s)^2*(n*m + n*s - m*s - 2*s^2);
  c := -s*(n - 2*s)*(n*m + n*s - m*s - 2*s^2)*(n*m + 4*n*s - 4*s^2)*(n*m + 4*n*s - 2*m*s - 4*s^2);
  d := -s*(n - s)*(m + 2*s)^2*(n - 2*s)^2*(n*m + 4*n*s - 4*s^2);
  return x*(x+a)*(x+b)*(x+c)*(x+d), [a,b,c,d];
end function;

// clear square factors from quintic: y^2 = f with f = x*prod(x+ai), rescale x -> lam^2 x
// to shrink coefficients: use minimal model machinery instead
pts := [ [143,218,120], [143,241,120], [133,109,60], [266,241,120], [266,218,143], [266,241,143], [408,437,143] ];
labels := [ "P1=known2", "P2", "P3", "P4", "P5", "P6", "P7=known1" ];

G2list := [];
crvs := [];
for i in [1..#pts] do
  T := pts[i];
  Av := Q!T[1]; Bv := Q!T[2]; Cv := Q!T[3];
  sv, mv, nv := EmailSMN(Av, Bv, Cv);
  // scale (s,m,n) to a primitive integer vector for smaller coefficients
  den := LCM([Denominator(sv), Denominator(mv), Denominator(nv)]);
  si := Integers()!(sv*den); mi := Integers()!(mv*den); ni := Integers()!(nv*den);
  g := GCD([si, mi, ni]); si := si div g; mi := mi div g; ni := ni div g;
  fq, abcd := QuinticModel(Q!si, Q!mi, Q!ni);
  okD := Discriminant(fq) ne 0;
  Cq := HyperellipticCurve(fq);
  g2 := G2Invariants(Cq);
  Append(~G2list, g2);
  Append(~crvs, Cq);
  printf "%o (A:B:C)=(%o:%o:%o): (s:m:n)=(%o:%o:%o) disc_ok=%o\n", labels[i], T[1],T[2],T[3], si, mi, ni, okD;
end for;

print "-- pairwise G2 comparison (which points give the same curve) --";
for i in [1..#pts] do
  for j in [i+1..#pts] do
    if G2list[i] eq G2list[j] then
      printf "SAME CURVE: %o and %o\n", labels[i], labels[j];
    end if;
  end for;
end for;

// distinct curves: pick first representative of each class
seen := [];
reps := [];
for i in [1..#pts] do
  new := true;
  for j in seen do if G2list[j] eq G2list[i] then new := false; break; end if; end for;
  if new then Append(~seen, i); Append(~reps, i); end if;
end for;
printf "distinct curves among the 7 points: %o (reps %o)\n", #reps, [labels[i] : i in reps];

// validate each distinct NEW curve (skip known ones P1, P7)
for i in reps do
  if i eq 1 or i eq 7 then continue; end if;
  printf "==== VALIDATE %o (A:B:C)=(%o:%o:%o) ====\n", labels[i], pts[i][1], pts[i][2], pts[i][3];
  Cq := crvs[i];
  Cmin := ReducedMinimalWeierstrassModel(Cq);
  fmin, hmin := HyperellipticPolynomials(Cmin);
  printf "minimal model: y^2 + (%o)y = %o\n", hmin, fmin;
  J := Jacobian(Cmin);
  tt0 := Cputime();
  Tor := TorsionSubgroup(J);
  printf "TORSION: %o (order %o)  [%o s]\n", Invariants(Tor), #Tor, Cputime(tt0);
  // simplicity certificates
  fdisc, hdisc := HyperellipticPolynomials(Cmin);
  D := Integers()!Discriminant(Cmin);
  nc := 0;
  for p in PrimesInInterval(29, 250) do
    if nc ge 4 then break; end if;
    if D mod p eq 0 then continue; end if;
    chi := Px!Reverse(Coefficients(EulerFactor(Jacobian(ChangeRing(Cmin, GF(p))))));
    if Degree(chi) ne 4 or not IsIrreducible(chi) then continue; end if;
    K<aK> := NumberField(chi); c12 := MinimalPolynomial(aK^12);
    if IsIrreducible(c12) and Degree(c12) eq 4 then
      printf "p=%o: simplicity cert OK\n", p; nc +:= 1;
    end if;
  end for;
  printf "certificates: %o\n", nc;
  if Invariants(Tor) eq [2,2,2,12] and nc ge 2 then
    printf "RESULT %o: NEW geometrically simple (2,2,2,12) curve CONFIRMED\n", labels[i];
  else
    printf "RESULT %o: check by hand (torsion %o, certs %o)\n", labels[i], Invariants(Tor), nc;
  end if;
end for;

print "ALL_DONE";
quit;
