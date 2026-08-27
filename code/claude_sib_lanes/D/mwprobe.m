// mwprobe.m — MW-lattice probe of C_rho' for the near-miss members.
// Enumerate P = c1 g1 + c2 g2 + c3 g3 + t (|ci| <= BND, t in torsion) on
// E = Jac(C_rho'), pull back to (u,y) on C, exact-test X1,X2,X3 squares.
// Reaches u-heights ~ exp(h) far beyond hyperellratpoints.
Q := Rationals(); P<u> := PolynomialRing(Q);
BND := StringToInteger(bnd);
rn := StringToInteger(rrn); rd := StringToInteger(rrd);
u0n := StringToInteger(uu0n); u0d := StringToInteger(uu0d);   // a known C-point
q := 4*u^2-6*u+3;
Pq := (q*rn - rd)*(q*rn - (2*u-1)*rd);
C := HyperellipticCurve(Pq);
u0 := u0n/u0d;
ok, y0 := IsSquare(Evaluate(Pq, u0));
if not ok then printf "BAD BASE POINT\n"; quit; end if;
pt0 := C![u0, y0];
E, m := EllipticCurve(C, pt0);
Em := MinimalModel(E);
_, toEm := IsIsomorphic(E, Em);
printf "member %o/%o: E = %o\n", rn, rd, aInvariants(Em);
G, f := MordellWeilGroup(Em);
gens := [f(G.i) : i in [1..Ngens(G)]];
T, mt := TorsionSubgroup(Em);
tors := [mt(t) : t in T];
gl := [g : g in gens | Order(g) eq 0];
r := #gl;
printf "rank gens %o, torsion %o\n", r, Invariants(T);
mi := Inverse(toEm) * Inverse(m);   // Em -> E -> C
function TestU(uu)
    // exact X1,X2,X3 test at u = uu (rational)
    p := Numerator(uu); qq := Denominator(uu);
    QQ := 4*p^2-6*p*qq+3*qq^2;
    C3 := 16*p^4-40*p^3*qq+40*p^2*qq^2-18*p*qq^3+3*qq^4;
    C2 := -16*p^4+32*p^3*qq-28*p^2*qq^2+10*p*qq^3-qq^4;
    C1 := (8*p^3-12*p^2*qq+10*p*qq^2-3*qq^3)*qq;
    C0 := (-2*p+qq)*qq^3;
    X1 := rn*(QQ*rn-(2*p-qq)*qq*rd);
    X2 := QQ*rn^2-(4*p^2-4*p*qq+2*qq^2)*rn*rd+(2*p-qq)*qq*rd^2;
    X3 := rn*(C3*rn^3+C2*rn^2*rd+C1*rn*rd^2+C0*rd^3);
    X4 := C3*rn^2+(-16*p^3+28*p^2*qq-18*p*qq^2+4*qq^3)*qq*rn*rd+(2*p-qq)^2*qq^2*rd^2;
    if X3 eq 0 or X4 eq 0 then return false, 0; end if;
    s1 := X1 gt 0 and IsSquare(X1);
    s2 := X2 gt 0 and IsSquare(X2);
    s3 := X3 gt 0 and IsSquare(X3);
    n := (s1 select 1 else 0) + (s2 select 1 else 0) + (s3 select 1 else 0);
    return n eq 3, n;
end function;
cnt := 0; best := 0;
for c1 in [-BND..BND] do
 for c2 in [-BND..BND] do
  for c3 in [-BND..BND] do
   for t in tors do
    P0 := c1*gl[1] + (r ge 2 select c2*gl[2] else Em!0)
          + (r ge 3 select c3*gl[3] else Em!0) + t;
    if P0 eq Em!0 then continue; end if;
    ptC := mi(P0);
    if ptC[3] eq 0 then continue; end if;
    uu := ptC[1]/ptC[3];
    if uu eq 0 or uu eq 1 or uu eq 1/2 then continue; end if;
    cnt +:= 1;
    hit, n := TestU(uu);
    if n gt best then best := n; end if;
    if hit then
        printf "MW-HIT member %o/%o u = %o (c=[%o,%o,%o])\n", rn, rd, uu, c1, c2, c3;
    end if;
   end for;
   if r lt 3 then break; end if;
  end for;
  if r lt 2 then break; end if;
 end for;
end for;
printf "member %o/%o BND=%o: tested %o MW points, best partial %o/3, done\n",
    rn, rd, BND, cnt, best;
quit;
