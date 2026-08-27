// ===========================================================================
// Lane 2 (claude_ov_c7z) -- Aut(X) for the genus-4 order-112 curve
//     X : u^2 = G1(T) = T(T-3)(T^2-24T+36)
//         v^2 = G2(T) = (T-3)(7T-12)(37T^2-96T+36)
// and the elliptic splitting of Jac(X).
//
// Magma's Automorphisms() on the genus-4 function field does not terminate in
// reasonable time.  Instead we determine Aut(X) structurally:
//   X -> P^1_T is Galois with group V = (Z/2)^2 (deck transformations
//   s1 : u -> -u,  s2 : v -> -v), so V <= Aut(X).
//   V is the unique subgroup with X/V = P^1 of that shape; any automorphism of X
//   normalises V (V is the kernel of the action on the 4-dim'l space H^0(K) ... more
//   robustly: we bound Aut(X) by the Mobius transformations preserving the BRANCH
//   DATA of X -> P^1, which every element of the normaliser of V must induce).
// Branch data: 7 points of P^1(Qbar),
//     type-1 (only u ramifies): roots of G1 other than 3      = {0} u {12 +- 6 sqrt3}
//     type-2 (only v ramifies): roots of G2 other than 3      = {12/7} u roots(37T^2-96T+36)
//     type-3 (both u,v ramify): T = 3
// A Mobius map induced by an automorphism must fix the type-3 point T=3 and permute
// {type-1} and {type-2} (possibly swapping the two classes, since s1<->s2 is allowed).
// ===========================================================================
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals();
PP<X> := PolynomialRing(Q);
G1 := X*(X-3)*(X^2-24*X+36);
G2 := (X-3)*(7*X-12)*(37*X^2-96*X+36);
GH := X*(7*X-12)*(X^2-24*X+36)*(37*X^2-96*X+36);

B1 := ExactQuotient(G1, X-3);            // type-1 branch points
B2 := ExactQuotient(G2, X-3);            // type-2 branch points
printf "type-1 branch poly = %o\ntype-2 branch poly = %o\ntype-3 branch point: T = 3\n", B1, B2;
K := SplittingField(B1*B2);
printf "splitting field of the branch locus: degree %o over Q, disc %o\n", Degree(K), Discriminant(MaximalOrder(K));
PK<Y> := PolynomialRing(K);
S1 := [ r[1] : r in Roots(PK!B1) ];
S2 := [ r[1] : r in Roots(PK!B2) ];
printf "#type-1 = %o, #type-2 = %o\n", #S1, #S2;
assert #S1 eq 3 and #S2 eq 3;

// A Mobius transformation is determined by the images of 3 points.  It must send
// 3 |-> 3, and either S1 -> S1 & S2 -> S2, or S1 -> S2 & S2 -> S1.
function MobiusFrom(src, dst)
  // the unique Mobius map sending src[i] |-> dst[i] (i=1..3), as a 2x2 matrix, or false
  M := Matrix(K,3,3,[]);
  // solve (a x + b) - y (c x + d) = 0 for the 3 pairs
  rows := [];
  for i in [1..3] do
    x := src[i]; y := dst[i];
    Append(~rows, [x, K!1, -y*x, -y]);
  end for;
  A := Matrix(K, 3, 4, &cat rows);
  ns := NullSpace(Transpose(A));
  if Dimension(ns) ne 1 then return false, _; end if;
  w := Eltseq(Basis(ns)[1]);
  a,b,c,d := Explode(w);
  if a*d - b*c eq 0 then return false, _; end if;
  return true, [a,b,c,d];
end function;

auts := [];
for swp in [false, true] do
  A0 := S1; B0 := swp select S2 else S1;
  Bo := swp select S1 else S2;
  for pm in Permutations({1,2,3}) do
    src := [ S1[1], S1[2], S1[3] ];
    dst := [ (swp select S2 else S1)[i] : i in pm ];
    ok, M := MobiusFrom(src, dst);
    if not ok then continue; end if;
    a,b,c,d := Explode(M);
    f := func< t | (c*t+d) eq 0 select Infinity() else (a*t+b)/(c*t+d) >;
    // must fix T = 3 and map the other class correctly
    if c*3+d eq 0 or (a*3+b)/(c*3+d) ne 3 then continue; end if;
    im2 := { f(t) : t in S2 };
    tgt2 := swp select Seqset(S1) else Seqset(S2);
    if im2 ne tgt2 then continue; end if;
    Append(~auts, <swp, M>);
  end for;
end for;
printf "Mobius transformations preserving the branch data: %o\n", #auts;
for t in auts do printf "   swap=%o  matrix %o\n", t[1], t[2]; end for;
printf "=> reduced automorphism group Aut(X)/V has order <= %o\n", #auts;
printf "=> |Aut(X)| <= 4 * %o = %o ; and Aut(X) >= V = (Z/2)^2 (the deck group)\n", #auts, 4*#auts;

printf "\n---- upper bound from reduction mod good primes ----\n";
for p in [5,11,13,17,19,23,29,31] do
  if (Integers()!Discriminant(HyperellipticCurve(GH))) mod p eq 0 then continue; end if;
  Fp := GF(p);
  KA := FunctionField(Fp); TA := KA.1;
  PA<YA> := PolynomialRing(KA);
  ok := true;
  try
    F1 := FunctionField(YA^2 - Evaluate(PolynomialRing(Fp)!G1, TA));
    PB<YB> := PolynomialRing(F1);
    F2 := FunctionField(YB^2 - (F1!Evaluate(PolynomialRing(Fp)!G2, TA)));
    printf "  p=%o : genus(X mod p) = %o, #Aut(X_Fp) = %o\n", p, Genus(F2), #Automorphisms(F2);
  catch e
    printf "  p=%o : %o\n", p, e`Object;
  end try;
end for;

printf "\n---- the elliptic splitting of Jac(X) ----\n";
CH := HyperellipticCurve(GH);
printf "Aut(H)/Q order = %o ; geometric = %o\n", #AutomorphismGroup(CH), IdentifyGroup(GeometricAutomorphismGroup(CH));
// bielliptic: find the involution T -> (aT+b)/(cT+d) preserving GH up to squares
found := false;
for cand in CartesianPower([-400..400],2) do
  // involutions of P^1 fixing the pair structure have the form T -> (m T + n)/(T - m)
  m := cand[1]; n := cand[2];
  if m^2 + n eq 0 then continue; end if;
  // pull back GH: (T-m)^6 * GH((mT+n)/(T-m)) must equal lambda * GH(T) for some lambda
  num := &+[ Coefficient(GH,i) * (m*X+n)^i * (X-m)^(6-i) : i in [0..6] ];
  if num eq 0 then continue; end if;
  qq := num / GH;
  if Denominator(qq) eq 1 and Degree(Numerator(qq)) eq 0 then
    lam := Q ! Numerator(qq);
    printf "BIELLIPTIC involution T -> (%o T + %o)/(T - %o)  with GH pullback factor %o (square %o)\n",
      m, n, m, lam, IsSquare(lam);
    found := true;
    break;
  end if;
end for;
if not found then printf "no small-height bielliptic involution found by this search\n"; end if;
printf "AUT_DONE\n";
quit;
