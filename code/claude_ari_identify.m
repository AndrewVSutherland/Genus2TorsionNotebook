// claude_ari_identify.m — classify primitive surface points by the genus-2 curve they induce.
// Reads triples from PTSFILE (lines "A B C"), computes (s,m,n) for both orientations,
// groups by G2Invariants, prints classes vs the two known curves.
SetColumns(0);
SetMemoryLimit(6*10^9);
Q := Rationals(); Px<x> := PolynomialRing(Q);
if not assigned PTSFILE then PTSFILE := "/tmp/claude-1000/-home-claude-torsion-jac/944b35de-89b9-4bb1-ac14-af36f0e69b0a/scratchpad/prims_plain.txt"; end if;

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
  return x*(x+a)*(x+b)*(x+c)*(x+d);
end function;

function CurveOf(Av, Bv, Cv)
  sv, mv, nv := EmailSMN(Q!Av, Q!Bv, Q!Cv);
  den := LCM([Denominator(sv), Denominator(mv), Denominator(nv)]);
  si := Integers()!(sv*den); mi := Integers()!(mv*den); ni := Integers()!(nv*den);
  g := GCD([si, mi, ni]); si := si div g; mi := mi div g; ni := ni div g;
  fq := QuinticModel(Q!si, Q!mi, Q!ni);
  if Discriminant(fq) eq 0 then return false, _, _; end if;
  return true, G2Invariants(HyperellipticCurve(fq)), [si, mi, ni];
end function;

// known curves
f2rec := 36*x^6+36750*x^5-462983772*x^4-301623595823*x^3+1518598238654317*x^2+397058962729817115*x-1282993930035013443975;
C2rec := HyperellipticCurve(f2rec, x^2+x);
f1rec := 756900*x^6 + 737595570*x^5 + 150572203590*x^4 - 15854483576121*x^3 - 530648977741620*x^2 + 32014154874551031*x + 830742747091037849;
C1rec := HyperellipticCurve(f1rec, x^2+1);
KnownG2 := [ G2Invariants(C1rec), G2Invariants(C2rec) ];
KnownNames := [ "curve#1", "curve#2" ];

lines := Split(Read(PTSFILE), "\n");
classG2 := [];   // list of G2 invariant triples of distinct curves
classPts := [];  // parallel: list of point-descriptions
for L in lines do
  w := Split(L, " ");
  if #w lt 3 then continue; end if;
  Av := StringToInteger(w[1]); Bv := StringToInteger(w[2]); Cv := StringToInteger(w[3]);
  for orient in [1, 2] do
    Aa := orient eq 1 select Av else Cv;
    Cc := orient eq 1 select Cv else Av;
    ok, g2, smn := CurveOf(Aa, Bv, Cc);
    if not ok then printf "(%o,%o,%o) orient%o: DEGENERATE quintic\n", Av,Bv,Cv,orient; continue; end if;
    name := "";
    for k in [1..#KnownG2] do if g2 eq KnownG2[k] then name := KnownNames[k]; end if; end for;
    idx := 0;
    for k in [1..#classG2] do if g2 eq classG2[k] then idx := k; end if; end for;
    if idx eq 0 then
      Append(~classG2, g2);
      Append(~classPts, []);
      idx := #classG2;
      if name eq "" then name := Sprintf("NEWCURVE_%o", idx); end if;
      printf "class %o = %o first seen at (%o,%o,%o) orient%o (s:m:n)=(%o:%o:%o)\n", idx, name, Av,Bv,Cv, orient, smn[1],smn[2],smn[3];
    end if;
    classPts[idx] := Append(classPts[idx], Sprintf("(%o,%o,%o)o%o", Av,Bv,Cv,orient));
    if name ne "" and idx le #KnownG2 then ; end if;
  end for;
end for;
print "== classes ==";
for k in [1..#classG2] do
  nm := "";
  for j in [1..#KnownG2] do if classG2[k] eq KnownG2[j] then nm := KnownNames[j]; end if; end for;
  if nm eq "" then nm := Sprintf("NEWCURVE_%o", k); end if;
  printf "class %o [%o]: %o points: %o\n", k, nm, #classPts[k], classPts[k];
end for;
print "ALL_DONE";
quit;
