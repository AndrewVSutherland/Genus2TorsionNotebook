// ===========================================================================
// Lane 2 (claude_ov_c7z) -- CHABAUTY on the genus-2 quotient H, and the
// resulting complete determination of X(Q) for the order-112 curve
//     X : u^2 = G1(T) = T(T-3)(T^2-24T+36)
//         v^2 = G2(T) = (T-3)(7T-12)(37T^2-96T+36)        [genus 4]
// (built and verified in code/claude_ov_c7z_xdecide.m).
//
// The (Z/2)^2 quotients are E1 = E2 = y^2 = x^3-12x (rank 1) and
//     H : y^2 = G1*G2/(T-3)^2 = T(7T-12)(T^2-24T+36)(37T^2-96T+36)  [genus 2].
// X(Q) minus the fibre T=3 maps injectively-on-T into H(Q), so H(Q) decides X.
// ===========================================================================
SetColumns(0);
if not assigned MemGB then MemGB := 12; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals();
PP<X> := PolynomialRing(Q);

G1 := X*(X-3)*(X^2-24*X+36);
G2 := (X-3)*(7*X-12)*(37*X^2-96*X+36);
GH := X*(7*X-12)*(X^2-24*X+36)*(37*X^2-96*X+36);
assert G1*G2 eq (X-3)^2*GH;
printf "G1 = %o\nG2 = %o\nGH = %o\n", G1, G2, GH;

CH := HyperellipticCurve(GH);
JH := Jacobian(CH);
printf "H genus = %o ; lc(GH)=%o square? %o (=> no rational points at infinity)\n",
   Genus(CH), LeadingCoefficient(GH), IsSquare(LeadingCoefficient(GH));
tor, tmap := TorsionSubgroup(JH);
printf "Jac(H) torsion = %o\n", Invariants(tor);
lo, hi := RankBounds(JH);
printf "Jac(H) RankBounds = %o .. %o\n", lo, hi;

printf "\n---- rational points of H by naive search ----\n";
ptsH := Points(CH : Bound := 200000);
printf "Points(Bound:=200000) : %o points\n%o\n", #ptsH, ptsH;

printf "\n---- build J(Q) ----\n";
JPts := [];
for i in [1..#ptsH] do
  for j in [1..#ptsH] do
    if i eq j then continue; end if;
    ok := true; D := JH!0;
    try D := JH ! [ptsH[i], ptsH[j]]; catch e ok := false; end try;
    if ok then Append(~JPts, D); end if;
  end for;
end for;
printf "%o divisor classes built from differences of the known points\n", #JPts;
inf := [ D : D in JPts | Order(D) eq 0 ];
printf "of infinite order: %o\n", #inf;
if #inf eq 0 then
  printf "RANK IS 0 (all differences torsion, and RankBounds upper bound 0..%o)\n", hi;
end if;

bas := ReducedBasis(JPts);
printf "ReducedBasis of the subgroup generated: %o element(s)\n", #bas;
for b in bas do printf "   gen %o   height %o\n", b, Height(b); end for;
rk := #bas;
printf "==> rank(J(Q)) = %o  (RankBounds %o..%o, found %o independent)\n", rk, lo, hi, #bas;
assert rk ge lo and rk le hi;

if rk eq 0 then
  printf "\n---- rank 0: H(Q) from the torsion of J ----\n";
  S := Chabauty0(JH);
  printf "CHABAUTY0 H(Q) = %o\n", S;
  HQ := S;
else
  printf "\n---- rank 1: saturate the generator, then Chabauty ----\n";
  P := bas[1];
  sat := bas;
  try
    sat := Saturation(bas, 100);
    printf "Saturation(.,100) gave %o generator(s)\n", #sat;
    for b in sat do printf "   sat gen %o height %o\n", b, Height(b); end for;
  catch e
    printf "Saturation failed: %o\n", e`Object;
  end try;
  P := sat[1];
  printf "using generator P = %o  height %o\n", P, Height(P);
  V, N := Chabauty(P);
  printf "CHABAUTY returned %o points, index-bound N = %o\n", #V, N;
  printf "CHABAUTY H(Q) = %o\n", V;
  HQ := V;
end if;

printf "\n---- pull back to X ----\n";
Tvals := [];
for p in HQ do
  c := Eltseq(p);
  if c[3] eq 0 then printf "  point at infinity on H: %o\n", p; continue; end if;
  Append(~Tvals, c[1]/c[3]);
end for;
Tvals := Sort(Setseq(Seqset(Tvals)));
printf "T-coordinates of H(Q): %o\n", Tvals;
XQ := [];
for tv in Tvals do
  s1, r1 := IsSquare(Evaluate(G1,tv));
  s2, r2 := IsSquare(Evaluate(G2,tv));
  printf "  T=%-8o G1=%-14o sq=%-5o   G2=%-16o sq=%o\n", tv, Evaluate(G1,tv), s1, Evaluate(G2,tv), s2;
  if s1 and s2 then Append(~XQ, tv); end if;
end for;
// the fibre T = 3 (base locus of X -> H) and T = infinity
printf "  T=3      G1=%o G2=%o  -> the point (3,0,0) on X\n", Evaluate(G1,3), Evaluate(G2,3);
printf "  T=oo     lc(G1)=%o (sq %o), lc(G2)=%o (sq %o)  -> %o rational point at infinity\n",
   LeadingCoefficient(G1), IsSquare(LeadingCoefficient(G1)),
   LeadingCoefficient(G2), IsSquare(LeadingCoefficient(G2)),
   (IsSquare(LeadingCoefficient(G1)) and IsSquare(LeadingCoefficient(G2))) select "IS a" else "NO";
XQ := Sort(XQ cat [Q!3]);
printf "X(Q) T-COORDINATES = %o\n", XQ;

printf "\n---- are they degenerate?  (a(T), b(T) and the five z's) ----\n";
KT := FieldOfFractions(PP);
TT := KT!X;
at := (-22/9*TT^2 + 80/3*TT - 40)/(TT^2 + 4*TT - 12);
bt := (7/27*TT^2 + 10/9*TT - 8/3)/(TT^2 - 2*TT);
printf "a(T) = %o  (poles at %o)\n", at, [ r[1] : r in Roots(PP!Denominator(at)) ];
printf "b(T) = %o  (poles at %o ; zeros at %o)\n", bt,
   [ r[1] : r in Roots(PP!Denominator(bt)) ], [ r[1] : r in Roots(PP!Numerator(bt)) ];
Pw<w> := PolynomialRing(Q);
for tv in XQ do
  dena := Evaluate(PP!Denominator(at), tv);
  denb := Evaluate(PP!Denominator(bt), tv);
  numb := Evaluate(PP!Numerator(bt), tv);
  if dena eq 0 or denb eq 0 then
    printf "  T=%-6o : a or b is INFINITE (point at infinity of Cab) -> NOT a configuration\n", tv;
    continue;
  end if;
  av := Evaluate(at,tv); bv := Evaluate(bt,tv);
  if bv eq 0 then
    printf "  T=%-6o : (a,b)=(%o,0) -> some w_i = 0, contradicts E4=1 -> NOT a configuration\n", tv, av;
    continue;
  end if;
  quart := (w^2-av*w+bv)*(w^2-(4-av)*w+1/bv);
  rts := [ r[1] : r in Roots(quart) ];
  mult := [ r[2] : r in Roots(quart) ];
  zs := [Q!2] cat &cat[ [rts[i] : k in [1..mult[i]]] : i in [1..#rts] ];
  printf "  T=%-6o : (a,b)=(%o,%o)  quartic=%o  roots %o (mult %o)\n", tv, av, bv, quart, rts, mult;
  printf "            five z = %o  -> distinct: %o\n", zs, #Seqset(zs) eq 5;
  if #Seqset(zs) eq 5 then
    printf "            *** NON-DEGENERATE ORDER-112 CANDIDATE AT T=%o ***\n", tv;
  end if;
end for;
printf "CHABAUTY_DONE\n";
quit;
