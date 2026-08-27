//////////////////////////////////////////////////////////////////////
// claude_ov_m612prym_ffprobe.m  (lane 9, overnight 2026-07-25)
//
// Probe: can we do Abel-Prym arithmetic on E8 over F_q?
//   - build the function field of E8 mod q, check Genus = 4
//   - build the Prym involution iota : y -> -y and push places through it
//   - compute the divisor class group and the order of a class
//   - locate the boundary places over x = 0
// Run: magma -b q:=31 claude_ov_m612prym_ffprobe.m
//////////////////////////////////////////////////////////////////////
SetColumns(0);
if not assigned q then q := 31; else q := StringToInteger(q); end if;
printf "=== ffprobe q = %o ===\n", q;

k := GF(q);
Fx<X> := RationalFunctionField(k);
Py<Y> := PolynomialRing(Fx);
Q8 := Y^8 + (216*X^4+72*X^3-24*X^2)*Y^4
    + (-1296*X^6-1728*X^5-432*X^4+64*X^3)*Y^2
    + (-3888*X^8-2592*X^7+432*X^6+288*X^5-48*X^4);
t0 := Cputime();
FF<yy> := FunctionField(Q8);
g := Genus(FF);
printf "genus = %o  (%o s)\n", g, Cputime(t0);

// the Prym involution
ok := true;
try
  iota := hom< FF -> FF | -yy >;
  printf "hom<FF->FF|-y> constructed OK\n";
catch err
  ok := false;
  printf "hom form 1 failed: %o\n", err`Object;
end try;
if not ok then
  iota := hom< FF -> FF | map<Fx->Fx | a :-> a>, -yy >;
  printf "hom form 2 constructed OK\n";
end if;
printf "check iota(y) = %o, iota(iota(y)) = %o\n", iota(yy), iota(iota(yy));

pl1 := Places(FF, 1);
printf "#places of degree 1 = %o\n", #pl1;

// how do we push a place through iota?
P := pl1[1];
printf "test place: x = %o, deg = %o\n", Evaluate(FF!X, P), Degree(P);
worked := "";
try
  D := Pullback(iota, Divisor(P));
  printf "Pullback(iota, Divisor(P)) OK: %o\n", D;
  worked := "pullback-div";
catch err
  printf "Pullback on divisor failed: %o\n", err`Object;
end try;
if worked eq "" then
  try
    D := Pullback(iota, P);
    printf "Pullback(iota, P) OK: %o\n", D;
    worked := "pullback-place";
  catch err
    printf "Pullback on place failed: %o\n", err`Object;
  end try;
end if;

t0 := Cputime();
Cl, mCl := ClassGroup(FF);
printf "ClassGroup = %o   (%o s)\n", Cl, Cputime(t0);
printf "invariants %o\n", Invariants(Cl);

// boundary places over x = 0
z0 := Zeros(Fx.1)[1];
dec := Decomposition(FF, z0);
printf "places over x=0: %o\n", [<Degree(d[1]), d[2]> : d in dec];

// class of P - iota(P) for a few places
for i in [1..Min(4,#pl1)] do
  P := pl1[i];
  Dp := Divisor(P);
  Di := Pullback(iota, Dp);
  c := (Dp - Di) @@ mCl;
  printf "place %o: x=%o  class(P - iota P) = %o, order %o\n",
      i, Evaluate(FF!X, P), c, Order(c);
end for;
print "FFPROBE_DONE";
quit;
