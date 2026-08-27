// claude_222_mw_pullback.m — enumerate Mordell-Weil points of the rank-1
// conductor-92 elliptic quotients of X1 (branch-point cubic-root cover) and
// X2 (Flynn43 extra-root cover), and test each fiber of the double covers
// X_i -> quotient for rational points.  A split fiber = rational point on X_i
// = a [1,1,2,2]-type member of the corresponding [22] family = [2,22] candidate.
// This reaches parameter heights far beyond any sieve (MW point heights grow
// quadratically in n).
SetColumns(0);
SetMemoryLimit(24*10^9);
Q := Rationals();

NMAX := 60;

procedure HuntCover(X, nm : seed := false, seedpt := [])
  Xp := ProjectiveClosure(X);
  G := AutomorphismGroup(Xp);
  Y, qm := CurveQuotient(G);
  printf "%o: quotient computed, genus %o\n", nm, Genus(Y);
  BS := BaseScheme(qm);
  pts := [];
  if seed then
    // transport a known point of X through the quotient map
    xpt := Xp!seedpt;
    ypt := qm(xpt);
    pts := [ ypt ];
    printf "%o: seeded quotient point via qm(%o) = %o\n", nm, seedpt, Eltseq(ypt);
  else
    pts := PointSearch(Y, 10^5);
  end if;
  if #pts eq 0 then printf "%o: no quotient point available\n", nm; return; end if;
  printf "%o: using base point %o on quotient\n", nm, pts[1];
  E, toE := EllipticCurve(Y, pts[1]);
  Em, mmap := MinimalModel(E);
  printf "%o: quotient elliptic = %o, conductor %o\n", nm, aInvariants(Em), Conductor(Em);
  MW, mMW := MordellWeilGroup(Em);
  gens := [ mMW(MW.i) : i in [1..Ngens(MW)] | Order(MW.i) eq 0 ];
  printf "%o: MW rank %o\n", nm, #gens;
  if #gens eq 0 then printf "%o: rank 0!\n", nm; return; end if;
  gen := gens[1];
  nhits := 0;
  for n in [-NMAX..NMAX] do
    if n eq 0 then continue; end if;
    P := n*gen;
    try
      // Em -> E -> Y by scheme pullbacks WITH per-stage image verification
      // (naive Pullback point lists include patch-indeterminacy junk)
      S1 := Pullback(mmap, Cluster(Codomain(mmap)!Eltseq(P)));
      p1 := [];
      for q in RationalPoints(S1) do
        try
          if mmap(Domain(mmap)!Eltseq(q)) eq Codomain(mmap)!Eltseq(P) then Append(~p1, q); end if;
        catch e ; end try;
      end for;
      if #p1 eq 0 then printf "%o n=%o transport E fail\n", nm, n; continue; end if;
      S2 := Pullback(toE, Cluster(Codomain(toE)!Eltseq(p1[1])));
      p2 := [];
      for q in RationalPoints(S2) do
        try
          if toE(Domain(toE)!Eltseq(q)) eq Codomain(toE)!Eltseq(p1[1]) then Append(~p2, q); end if;
        catch e ; end try;
      end for;
      if #p2 eq 0 then printf "%o n=%o transport Y fail\n", nm, n; continue; end if;
      yP := p2[1];
      F := Pullback(qm, Cluster(Codomain(qm)!Eltseq(yP)));
      rp := [];
      for r in RationalPoints(F) do
        if r in BS then continue; end if;
        // verify the candidate genuinely maps to yP (Pullback picks up
        // indeterminacy points of the polynomial patch)
        ok := true;
        try
          im := qm(Domain(qm)!Eltseq(r));
          ok := Eltseq(im) eq Eltseq(Codomain(qm)!Eltseq(yP));
        catch e
          ok := false;
        end try;
        if ok then Append(~rp, r); end if;
      end for;
      if #rp gt 0 then
        nhits +:= 1;
        printf "SPLIT FIBER %o at n=%o: X-points %o\n", nm, n, [Eltseq(r) : r in rp];
      end if;
    catch e
      printf "%o n=%o fiber failed\n", nm, n;
    end try;
  end for;
  printf "%o: split fibers among |n| <= %o: %o\n", nm, NMAX, nhits;
end procedure;

// X1
A2<s,x> := AffineSpace(Q, 2);
c1 := (s-1)^2*x^3 + (2*s^4-4*s^3+5*s^2-4*s+1)*x^2
    + (s^6-2*s^5+3*s^4-2*s^3+2*s^2-2*s+1)*x - s^2*(s^2-s+1)^2;
X1 := Curve(A2, c1);
HuntCover(X1, "X1");

// X2
FF<T> := FunctionField(Q);
uu := (-11/100*T^3 + 919/50*T^2 - 8017/25*T + 40154/25)/(T^3 - 23*T^2 + 168*T - 396);
vv := (2401/1600*T^4 - 12593/200*T^3 + 205647/200*T^2 - 385757/50*T + 2253001/100)/(T^4 - 34*T^3 + 421*T^2 - 2244*T + 4356);
tt := (-21609/40000*T^6 + 435561/10000*T^5 - 2897123/2000*T^4 + 6390763/250*T^3 - 126932427/500*T^2 + 844854361/625*T - 1894773841/625)/(T^6 - 46*T^5 + 865*T^4 - 8520*T^3 + 46440*T^2 - 133056*T + 156816);
PxF<X> := PolynomialRing(FF);
FlT := X^6+2*X^5+(2*tt+3)*X^4+2*X^3+(tt^2+1)*X^2+2*tt*(1-tt)*X+tt^2;
Q4 := FlT div (X^2 + FF!uu*X + FF!vv);
cfs := Coefficients(Q4);
den := LCM([Denominator(cq) : cq in cfs]);
R2<t2,r2> := PolynomialRing(Q, 2);
num := &+[ Evaluate(Numerator(cfs[i]*den), t2) * r2^(i-1) : i in [1..#cfs] ];
A2b<tb,rb> := AffineSpace(Q,2);
comp9 := 0;
for fa in Factorization(Evaluate(num, [tb, rb])) do
  if TotalDegree(fa[1]) ge 8 then comp9 := fa[1]; end if;
end for;
X2 := Curve(A2b, comp9);
// known boundary point of X2: (tau, r) = (58/3, 0)
assert Evaluate(comp9, [58/3, 0]) eq 0;
HuntCover(X2, "X2" : seed := true, seedpt := [58/3, 0, 1]);
print "MW_PULLBACK_DONE";
quit;
