//////////////////////////////////////////////////////////////////////
//  Symbolic setup for the p=0 boundary chart in the M_1(8,2,2)
//  halving problem.
//
//  At p=0 the marked point Q=(-1,p(s+1)) collapses to the double
//  branch point X=-1.  Away from s=1,-2, the square-quartic condition
//  forces a(-1)=0, so on the boundary we write
//
//      a(X) = (X+1)(X+B),
//      root(X) = (X+1)(M X+N).
//
//  For the transverse chart we use original split variables
//
//      u=t,  v=s+tV,
//
//  and allow first-order variation in A,B,C,D,E.
//////////////////////////////////////////////////////////////////////

Q := Rationals();

R0<s,B,M,N> := PolynomialRing(Q, 4, "grevlex");
P0<X> := PolynomialRing(R0);

f0 := (X+1)^2*(X+s^2+s+1)*((s-1)*X-1)*((s+2)*X+1);
L0 := s^2+s-2;
boundary_res := ExactQuotient(f0, (X+1)^2)
                - L0*(X+1)*(X+B)^2
                - (M*X+N)^2;
boundary_eqs := [Coefficient(boundary_res, i) : i in [0..2]];

print "p=0 factorization";
print Factorization(f0);
print "boundary square equations in s,B,M,N";
for i in [1..#boundary_eqs] do
    print "g", i-1, boundary_eqs[i];
end for;

if assigned do_boundary_elimination then
    if Type(do_boundary_elimination) eq MonStgElt then
        do_boundary_elimination := do_boundary_elimination in {"true", "True", "1", "yes"};
    end if;
end if;
if assigned do_boundary_elimination and do_boundary_elimination then
    Rlex<M1,N1,s1,B1> := PolynomialRing(Q, 4, "lex");
    phi := hom<R0 -> Rlex | s1, B1, M1, N1>;
    I := ideal<Rlex | [phi(g) : g in boundary_eqs]>;
    print "Computing boundary elimination...";
    time G := GroebnerBasis(I);
    for g in G do
        if Degree(g, M1) eq 0 and Degree(g, N1) eq 0 then
            print "boundary elimination", Factorization(g);
        end if;
    end for;
end if;

R<sT,t,V,BT,UA,UB,M0,N0,UC,UD,UE> := PolynomialRing(Q, 11, "grevlex");
P<X> := PolynomialRing(R);
u := t;
v := sT + t*V;
ss := u+v;
pp := u*v;
qtilde := -X^2 + (pp*ss - ss^2 + 2*pp - ss - 2)*X
          - (ss^2 - pp + ss + 1);
f := ((pp-ss+1)*X^2 + (2-ss)*X + 1)*((ss+2)*X + 1)*qtilde;
L := Coefficient(f, 5);

Afull := BT + 1 + t*UA;
Bfull := BT + t*UB;
Cfull := M0 + t*UC;
Dfull := M0 + N0 + t*UD;
Efull := N0 + t*UE;

a := X^2 + Afull*X + Bfull;
root := Cfull*X^2 + Dfull*X + Efull;
res := f - L*(X+1)*a^2 - root^2;
eqs := [Coefficient(res, i) : i in [0..4]];
e0 := [Evaluate(e, [sT,0,V,BT,UA,UB,M0,N0,UC,UD,UE]) : e in eqs];
e1 := [Evaluate(Coefficient(e, t, 1),
        [sT,0,V,BT,UA,UB,M0,N0,UC,UD,UE]) : e in eqs];

print "corrected first-order p=0 chart";
print "variables: s,t,V,B,UA,UB,M,N,UC,UD,UE";
print "constant equation term counts", [#Terms(g) : g in e0];
print "first-order equation term counts", [#Terms(g) : g in e1];
for i in [1..#e1] do
    print "e1", i-1, e1[i];
end for;

quit;
