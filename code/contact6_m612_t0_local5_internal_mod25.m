//////////////////////////////////////////////////////////////////////
// Exact mod-25 lift test for the smooth compensated omega-pole chart.
//
// This treats pole_order=1:
//   s=25*S, omega=O/5, m^2-c4=5*D.
// All rank-5 mod-5 residues have U+2*nu=0.  We solve the exact affine
// first-lift equations and ask whether that last contact divisor is left at
// order one.  A rank-5 escaping lift is a genuine smooth Q_5 cone.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Z:=Integers(); k:=GF(5);
R<S,O,m,D,L,U,nu>:=PolynomialRing(Z,7,"grevlex");
a:=S*O^2-3; b:=25*S-3; M:=L^2;
c5:=50*S; c4:=b^2+2*a-15; c3:=2*a*b+22;
c2:=a^2+2*b-15; c1:=2*(a+3);
G:=m^2-c4-5*D;
H:=(8*S*(c3-20*S*O)-D^2)^2
   -256*S^2*m^2*(c2+O*D);
B3:=c5*M+3*U;
Delta3:=4*c4*M+12*(U^2+nu^2)-B3^2;
F3:=B3*Delta3+16*nu^3-8*c3*M-8*U^3-48*U*nu^2;
F2:=Delta3^2+64*B3*nu^3-64*c2*M
    -192*(U^2*nu^2+nu^4);
F1:=Delta3*nu^3-4*c1*M-12*U*nu^4;
eqs:=[G,H,F1,F2,F3];
derivs:=[[Derivative(f,j):j in [1..7]]:f in eqs];

rank5_residues:=0; affine_lifts:=0; escaping:=0;
by_SO:=AssociativeArray(); samples:=[];
for s0 in [0..4] do for o0 in [0..4] do for m0 in [0..4] do
for d0 in [0..4] do for l0 in [0..4] do for u0 in [0..4] do
for n0 in [0..4] do
    x0:=[s0,o0,m0,d0,l0,u0,n0];
    if s0*o0*m0*l0*n0*(u0-2*n0) eq 0 then continue; end if;
    f0:=[Z!Evaluate(f,x0):f in eqs];
    if not &and[z mod 5 eq 0:z in f0] then continue; end if;
    A:=Matrix(k,5,7,
        &cat[[k!((Z!Evaluate(derivs[i][j],x0)) mod 5):j in [1..7]]
              :i in [1..5]]);
    if Rank(A) ne 5 then continue; end if;
    rank5_residues+:=1;
    rhs:=Vector(k,[k!((-(f0[i] div 5)) mod 5):i in [1..5]]);
    for s1 in [0..4] do for o1 in [0..4] do for m1 in [0..4] do
    for d1 in [0..4] do for l1 in [0..4] do for u1 in [0..4] do
    for n1 in [0..4] do
        x1:=[s1,o1,m1,d1,l1,u1,n1];
        if Vector(k,x1)*Transpose(A) ne rhs then continue; end if;
        affine_lifts+:=1;
        x:=[x0[i]+5*x1[i]:i in [1..7]];
        assert &and[(Z!Evaluate(f,x)) mod 25 eq 0:f in eqs];
        if (x[6]+2*x[7]) mod 25 eq 0 then continue; end if;
        escaping+:=1;
        key:=Sprint(<s0,o0>);
        if IsDefined(by_SO,key) then by_SO[key]+:=1; else by_SO[key]:=1; end if;
        if #samples lt 32 then
            Append(~samples,<x0,x1,x,
                [(Z!Evaluate(f,x)) mod 25:f in eqs],
                (x[6]+2*x[7]) mod 25>);
        end if;
    end for; end for; end for; end for; end for; end for; end for;
end for; end for; end for; end for; end for; end for; end for;

print "CONTACT6_M612_T0_LOCAL5_INTERNAL_MOD25";
print "RANK5_RESIDUES",rank5_residues,"AFFINE_MOD25_LIFTS",affine_lifts,
      "ESCAPING_U_PLUS_2NU",escaping;
print "ESCAPING_COUNTS_BY_S_O_RESIDUE",
      Sort([<x,by_SO[x]>:x in Keys(by_SO)]);
print "SAMPLES_x0_x1_xmod25_residuals_Um",samples;
print "CONTACT6_M612_T0_LOCAL5_INTERNAL_MOD25_DONE";
quit;
