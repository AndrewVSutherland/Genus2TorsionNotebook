// Lane 2 (claude_ov_c7z): COMPLETE determination of the special slices of the contact-7
// three-root surface.  Slice z_1 = c gives the quartic pencil
//     P_{m}(w) = w^4 - e1 w^3 + m w^2 - (A m + B) w + e4,   e1=6-c, e4=2/c,
//     A=(T-2)/(T-1), B=(T(c-5+2/c)-3c+14)/(T-1),  T = 3-1/(1-c).
// The three-root slice curve C_c = {(x,y): P has x and y as roots, x != y} is a (3,3) curve
// in P1xP1 of arithmetic genus 4; it is SMOOTH (genus 4) unless the branch divisor of
// M: w -> m degenerates, i.e. unless  disc_m( disc_w( P_m ) )  vanishes.
// This computes that double discriminant as a polynomial in c and factors it: its rational
// roots are exactly the candidate special slices.
SetColumns(0);
SetMemoryLimit(20*10^9);
Q := Rationals();
Fc<c> := FunctionField(Q);
Rm<m> := PolynomialRing(Fc);
Rw<w> := PolynomialRing(Rm);
T := 3 - 1/(1-c);
e1 := 6-c; e4 := 2/c;
A := (T-2)/(T-1);
B := (T*(c-5+2/c) - 3*c + 14)/(T-1);
P := w^4 - e1*w^3 + m*w^2 - (A*m+B)*w + e4;
D1 := Discriminant(P);          // in Rm = Q(c)[m]
printf "deg_m disc_w(P) = %o\n", Degree(D1);
D2 := Discriminant(D1);         // in Q(c)
printf "double discriminant computed\n";
Nn := Numerator(D2); Dd := Denominator(D2);
Pc<C> := PolynomialRing(Q);
NN := Pc![ Coefficient(Nn,i) : i in [0..Degree(Nn)] ];
printf "numerator degree %o\n", Degree(NN);
fa := Factorisation(NN);
printf "SPECIALC factorisation of the double discriminant numerator:\n";
for g in fa do
  rts := [ r[1] : r in Roots(g[1]) ];
  printf "  factor deg %o mult %o : %o   rational roots: %o\n", Degree(g[1]), g[2], g[1], rts;
end for;
allrts := &cat[ [ r[1] : r in Roots(g[1]) ] : g in fa ];
printf "SPECIALC rational roots (candidate special slices): %o\n", Sort(Setseq(Seqset(allrts)));
printf "SPECIALC leading/denominator info: den=%o\n", Dd;
printf "SPECIALC_DONE\n";
quit;
