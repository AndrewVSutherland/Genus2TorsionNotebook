// Lane 2 (claude_ov_c7z): walk the RANK-1 Mordell-Weil group of the c=2 slice's
// three-root curve and test the ORDER-112 (split-all) condition on every point.
// This is the INDEPENDENT confirmation of the Chabauty verdict in
// code/claude_ov_c7z_chabauty.m: it reaches T of astronomically large height
// (h(nP) ~ 0.2506 n^2), far beyond any sieve.
//
//   Slice z_1 = 2: the other four z are the roots of  P_m(w) = w^4-4w^3+m w^2-(2/3)m w+1.
//   Cab (genus 0) parametrized by T; the square classes (claude_ov_c7z_xdecide.m) are
//        d1 ~ G1(T) = T(T-3)(T^2-24T+36)        [ >=2 rational w  ->  [2,2,14] ]
//        d2 ~ G2(T) = (T-3)(7T-12)(37T^2-96T+36)[ the other two   ->  [2,2,2,14] ]
//   E1 : u^2 = G1(T) is genus 1.  Explicit reduction (root T=0, s = 1/T, w = u/T^2):
//        w^2 = -108 s^3 + 108 s^2 - 27 s + 1 ,  then  Y = 108 w, Xc = -108 s  gives
//        E : Y^2 = Xc^3 + 108 Xc^2 + 2916 Xc + 11664 ,   T = -108/Xc ,  u = 108 Y / Xc^2.
//
// FIXED 2026-07-25: the original script built the (3,3) plane correspondence curve and
// called EllipticCurve/Inverse on it, which exhausted 8 GB and left G uninitialised.
// The explicit map above is closed-form and costs nothing.
SetColumns(0);
if not assigned NMAX then NMAX := 250; elif Type(NMAX) eq MonStgElt then NMAX := StringToInteger(NMAX); end if;
SetMemoryLimit(4*10^9);
Q := Rationals();
PP<X> := PolynomialRing(Q);
G1 := X*(X-3)*(X^2-24*X+36);
G2 := (X-3)*(7*X-12)*(37*X^2-96*X+36);

E := EllipticCurve([0,108,0,2916,11664]);
Emin, iso := MinimalModel(E);
printf "E (three-root curve of the c=2 slice) = %o\n", E;
printf "minimal model %o   conductor %o   j = %o\n", Emin, Conductor(Emin), jInvariant(Emin);
G, mg, pr1, pr2 := MordellWeilGroup(Emin);
printf "MW = %o proven=(%o,%o)\n", Invariants(G), pr1, pr2;
gens := [ mg(g) : g in Generators(G) ];
printf "gens = %o  heights %o\n", gens, [Height(g) : g in gens];
Pinf := [ g : g in gens | Order(g) eq 0 ][1];
Tors := [ Emin!0 ] cat [ g : g in gens | Order(g) eq 2 ];
printf "free generator %o (height %o); torsion reps %o\n", Pinf, Height(Pinf), Tors;
isoi := Inverse(iso);

// verify the explicit birational map on the known points of X
printf "-- map check --\n";
for tv in [Q| 2, -6, 12/7, 3 ] do
  Xc := -108/tv;
  yy2 := Xc^3 + 108*Xc^2 + 2916*Xc + 11664;
  s, Yv := IsSquare(yy2);
  printf "   T=%-6o Xc=%-8o Y^2=%-10o square=%-5o", tv, Xc, yy2, s;
  if s then
    uu := (Xc eq 0) select 0 else 108*Yv/Xc^2;
    printf "  u=%o  u^2=G1(T)? %o", uu, uu^2 eq Evaluate(G1,tv);
  end if;
  printf "\n";
end for;

printf "-- walking n*P + t, |n| <= %o --\n", NMAX;
n3 := 0; n5 := 0; hits := [];
for n in [-NMAX..NMAX] do
  for t in Tors do
    P := n*Pinf + t;
    if P eq Emin!0 then continue; end if;
    Pe := isoi(P);
    c := Eltseq(Pe);
    if c[3] eq 0 then continue; end if;
    Xc := c[1]/c[3]; Yv := c[2]/c[3];
    if Xc eq 0 then continue; end if;          // T = infinity
    tv := -108/Xc;
    // by construction G1(tv) is a square; sanity-check it cheaply for small n
    if Abs(n) le 3 then
      assert IsSquare(Evaluate(G1,tv));
    end if;
    n3 +:= 1;
    g2v := Evaluate(G2,tv);
    if g2v eq 0 then Append(~hits, <n, t eq Emin!0, tv, "G2=0">); n5 +:= 1; continue; end if;
    nu := Numerator(g2v); de := Denominator(g2v);
    if nu*de lt 0 then continue; end if;
    if IsSquare(nu*de) then
      n5 +:= 1;
      Append(~hits, <n, t eq Emin!0, tv, "square">);
      printf "SPLITALL_HIT n=%o tors=%o T=%o\n", n, t, tv;
    end if;
  end for;
end for;
printf "MWSCAN_DONE NMAX=%o  three-root points walked=%o  split-all=%o\n", NMAX, n3, n5;
printf "split-all T-values: %o\n", [ h[3] : h in hits ];
printf "(these must be exactly the degenerate T in {-6,0,12/7,2,3}; see claude_ov_c7z_chabauty.m)\n";
quit;
