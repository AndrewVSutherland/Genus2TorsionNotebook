// claude_ov_88chi_w.m -- Lane 5 (2026-07-25).  Solve the Lambda_334 halving system
// above T2 = [L2] SYMBOLICALLY over Q(s,t,v): the four order-4 Mumford polynomials
//     u(x) = x^2 + w x - (w t^2 + b c)
// are the roots w of  disc( (k u^2 + d2 L1 L3)/L2 ) = 0,  k = -d2 (bc-a)^2/(e2 t^2 w^2).
// Prints the factorization of that discriminant in w, the explicit rational roots,
// and, for each, the norm characters  Res(L1,u), Res(L2,u), Res(L3,u)  as square
// classes over Q(s,t,v) -- i.e. the symbolic lift layer.
SetColumns(0);
SetMemoryLimit(24*10^9);
Qq := Rationals();
R<s,t,v> := PolynomialRing(Qq, 3);
F := FieldOfFractions(R);
Px<x> := PolynomialRing(F);

A  := s^2 - t^4 + t^2;
uu := (-s^2*A*v^2 - 2*A*v - 1)/(-s^2*t*A*v^2 + t);
a  := A/(1 - t^2);
b  := A/(uu^2*s^2 + 1 - t^2);
c  := t^2;
d2 := A*(s^2*uu^2 + t^4 - 2*t^2 + 1)
        *(s^4*uu^2 - s^2*t^2*uu^2 + s^2*uu^2 - t^6 + 3*t^4 - 3*t^2 + 1);
e2 := -a + b + c - 1;
e1 := 2*a - 2*b*c;
e0 := a*b*c - a*b - a*c + b*c;
L1 := x^2 + (e1/e2)*x + (e0/e2);
L2 := x^2 - b*c;
L3 := x^2 - a;
Delta := s^2*(A*v+1)^2 + t^2*(1-t^2);
b0 := A*t*(s^2*A*v^2 - 1)/Delta;
Psi := v*(1 + s^2*v)*(1 + A*v);

// squarefree part of an element of F, as (rational content class, polynomial part)
function SqClass(el)
  if el eq 0 then return 0, R!1; end if;
  nm := Numerator(el); dn := Denominator(el);
  pr := nm*dn;
  fa := Factorization(pr);
  poly := R!1;
  for g in fa do if IsOdd(g[2]) then poly *:= g[1]; end if; end for;
  // rational content
  cont := pr / (&*[R| g[1]^g[2] : g in fa ]);
  cq := Qq ! cont;
  nn := Numerator(cq)*Denominator(cq);
  sg := Sign(nn); nn := AbsoluteValue(nn);
  cc := sg;
  for pf in Factorization(nn) do if IsOdd(pf[2]) then cc *:= pf[1]; end if; end for;
  return cc, poly;
end function;

procedure ShowSq(name, el)
  cc, poly := SqClass(el);
  printf "%o  square class:  %o * (%o)\n", name, cc, poly;
end procedure;

Fw<w> := FunctionField(F);
Pw<X> := PolynomialRing(Fw);
L1w := Pw ! [ Fw!cc : cc in Coefficients(L1) ];
L2w := Pw ! [ Fw!cc : cc in Coefficients(L2) ];
L3w := Pw ! [ Fw!cc : cc in Coefficients(L3) ];
uw  := X^2 + w*X - (w*(Fw!t)^2 + Fw!(b*c));
kw  := -(Fw!d2)*(Fw!(b*c-a))^2/((Fw!e2)*(Fw!t)^2*w^2);
Rw  := kw*uw^2 + (Fw!d2)*L1w*L3w;
qq  := Rw div L2w;
dsc := Coefficient(qq,1)^2 - 4*Coefficient(qq,0)*Coefficient(qq,2);
dn2 := Numerator(dsc);
PolW := Parent(dn2);
printf "=== disc(quotient) numerator: factorization in w over Q(s,t,v) ===\n";
lins := [];
for g in Factorization(dn2) do
  printf "--- factor: degree %o in w, multiplicity %o\n", Degree(g[1]), g[2];
  printf "%o\n", g[1];
  if Degree(g[1]) eq 1 then
    root := -Coefficient(g[1],0)/Coefficient(g[1],1);
    Append(~lins, F!root);
  end if;
end for;
printf "=== %o rational roots w in Q(s,t,v) ===\n", #lins;

for i in [1..#lins] do
  w0 := lins[i];
  printf "\n=== ROOT %o ===\n", i;
  ShowSq("w0", w0);
  u0 := x^2 + w0*x - (w0*t^2 + b*c);
  // verify the halving identity fully
  k0 := -d2*(b*c-a)^2/(e2*t^2*w0^2);
  Rz := k0*u0^2 + d2*L1*L3;
  qz := Rz div L2;
  printf "  exact-division check: %o\n", Rz - qz*L2 eq 0;
  dz := Coefficient(qz,1)^2 - 4*Coefficient(qz,0)*Coefficient(qz,2);
  printf "  disc(quotient) = %o\n", dz;
  okB, Bv := IsSquare(k0 + d2);
  printf "  k0+d2 is a square in Q(s,t,v): %o\n", okB;
  N1 := Resultant(L1, u0);  N2 := Resultant(L2, u0);  N3 := Resultant(L3, u0);
  ShowSq("  Res(L1,u)", N1);
  ShowSq("  Res(L2,u)", N2);
  ShowSq("  Res(L3,u)", N3);
  ShowSq("  Psi", Psi);
  ok1 := IsSquare(N1*Psi);  ok2 := IsSquare(N2*Psi);
  printf "  Res(L1,u)*Psi is a square: %o ;  Res(L2,u)*Psi is a square: %o ;  Res(L3,u) is a square: %o\n",
         ok1, ok2, IsSquare(N3);
  // the split-component representative and the K3 representative
  ShowSq("  u(r)  = w0*t*(b0-t)", w0*t*(b0-t));
  ShowSq("  u(-r) = -w0*t*(b0+t)", -w0*t*(b0+t));
end for;

printf "W_DONE\n";
quit;
