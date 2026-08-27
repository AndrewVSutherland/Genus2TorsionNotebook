// Lane 2: the c=2 slice of the contact-7 three-root surface is the rank-1 elliptic curve
//   E : y^2 = x^3 - 12x  (conductor 288),  MW = Z/2 x Z.
// Every rational point gives z = (2, x, y) with three rational z's, i.e. a genus-2 curve
// with quintic factor type [1,1,1,2] -> 2-rank 3 -> torsion contains [2,2,14].
// ORDER 112 ([2,2,2,14]) <=> the residual quadratic splits <=> D := s^2 - 4p is a SQUARE,
//   s = 4 - (x+y),  p = 1/(x*y).
// This script (a) enumerates E(Q) far beyond any brute-force height and tests D = square,
//            (b) computes the divisor of D on E, hence the genus of the double cover w^2 = D.
SetColumns(0);
Q := Rationals();
A2<x,y> := AffineSpace(Q,2);
F := -x^3*y^2 + 2/3*x^3*y - x^2*y^3 + 14/3*x^2*y^2 - 8/3*x^2*y + 2/3*x*y^3 - 8/3*x*y^2 + x + y - 2/3;
C := Curve(A2,F);
PC := ProjectiveClosure(C);
P0 := PC ! [-1/9,-7/3,1];
E, mp := EllipticCurve(PC, P0);
Emin, iso := MinimalModel(E);
isoi := Inverse(iso);
printf "E (minimal) = %o   conductor %o\n", Emin, Conductor(Emin);
G, mg := MordellWeilGroup(Emin);
printf "MW invariants = %o\n", Invariants(G);
gens := [ mg(g) : g in Generators(G) ];
printf "gens = %o\n", gens;
// identify the free generator and the torsion generator
Tor := [ P : P in gens | Order(P) ne 0 ];
Fre := [ P : P in gens | Order(P) eq 0 ];
printf "torsion gens = %o ; free gens = %o\n", Tor, Fre;
PG := Fre[1];
tors := [ Emin!0 ] cat Tor;

mi := Inverse(mp);
NMAX := 60;
printf "enumerating n*PG + t for |n| <= %o, t in torsion (%o points)\n", NMAX, (2*NMAX+1)*#tors;
nhit := 0; ndeg := 0; ngood := 0;
sqfreecount := AssociativeArray();
for n in [-NMAX..NMAX] do
  for t in tors do
    P := n*PG + t;
    if P eq Emin!0 then continue; end if;
    ok := true;
    pt := 0;
    try pt := mi(isoi(P)); catch e ok := false; end try;
    if not ok then continue; end if;
    cc := Eltseq(pt);
    if cc[3] eq 0 then ndeg +:= 1; continue; end if;
    xx := cc[1]/cc[3]; yy := cc[2]/cc[3];
    if xx eq 0 or yy eq 0 or xx eq 1 or yy eq 1 then ndeg +:= 1; continue; end if;
    s := 4 - xx - yy; p := 1/(xx*yy);
    D := s^2 - 4*p;
    ngood +:= 1;
    if D ne 0 and IsSquare(D) then
      nhit +:= 1;
      printf "SPLITALL_HIT n=%o t=%o  (x,y)=(%o,%o)  D=%o\n", n, t, xx, yy, D;
    end if;
    if Abs(n) le 6 then
      nu := Numerator(D); de := Denominator(D);
      printf "  n=%o t=%o  x=%o y=%o  D=%o  sqfree=%o\n", n, t, xx, yy, D, Squarefree(nu*de);
    end if;
  end for;
end for;
printf "ENUM_DONE NMAX=%o good=%o degenerate=%o splitall_hits=%o\n", NMAX, ngood, ndeg, nhit;
quit;
