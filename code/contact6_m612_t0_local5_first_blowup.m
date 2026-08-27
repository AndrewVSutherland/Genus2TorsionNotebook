//////////////////////////////////////////////////////////////////////
// Deflated first blow-up of the unique rank-3 escaping p=5 cone of the
// integral T0-halved contact-6 core pullback.
//
// The center is
//   (s,omega,m,L,U,nu) = (3,1,0,2,2,-1) mod 5,
// on DB=DC=m=(U+2*nu)=(omega^2-1)=0.
// We put x=x0+5*y.  The four leading equations have rank three.  An
// integral row combination spanning the left kernel is divisible by an
// additional 5; replacing one equation by that deflated combination is an
// invertible Z_5 row operation and exposes the true first exceptional fiber.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Z:=Integers(); k:=GF(5);
R<s,omega,m,L,U,nu>:=PolynomialRing(Z,6,"grevlex");

a:=s*omega^2-3; b:=s-3; W:=m^2; M:=L^2;
c5:=2*s; c4:=b^2+2*a-15; c3:=2*a*b+22;
c2:=a^2+2*b-15; c1:=2*(a+3);
N:=c2-omega*c4+omega*W;
R8:=8*s*(c3-4*s*omega)-(c4-W)^2;
Fhalf:=R8^2-256*s^2*W*N;
B3:=c5*M+3*U;
Delta3:=4*c4*M+12*(U^2+nu^2)-B3^2;
F3:=B3*Delta3+16*nu^3-8*c3*M-8*U^3-48*U*nu^2;
F2:=Delta3^2+64*B3*nu^3-64*c2*M
    -192*(U^2*nu^2+nu^4);
F1:=Delta3*nu^3-4*c1*M-12*U*nu^4;
eqs:=[Fhalf,F1,F2,F3];
derivs:=[[Derivative(f,j):j in [1..6]]:f in eqs];
DB:=(a-3)^2-8*(b+3);
DC:=(b-3)^2-8*(a+3);

x0:=[3,1,0,2,2,-1];
J:=Matrix(k,4,6,
    &cat[[k!((Z!Evaluate(derivs[i][j],x0)) mod 5):j in [1..6]]
          :i in [1..4]]);
assert Rank(J) eq 3;
// Magma's Kernel(J) is the row/left kernel {lambda:lambda*J=0}.
Kleft:=Kernel(J);
assert Dimension(Kleft) eq 1;
lambdaF:=Basis(Kleft)[1];
lambda:=[Z!lambdaF[i]:i in [1..4]];
replace:=[i:i in [1..4]|lambda[i] mod 5 ne 0][1];
keep:=[i:i in [1..4]|i ne replace];

T<S,O,mu,E,V,Nn>:=PolynomialRing(Z,6,"grevlex");
subs:=[3+5*S,1+5*O,5*mu,2+5*E,2+5*V,-1+5*Nn];

function Strip5(f)
    ord:=0;
    while f ne 0 and &and[(Z!c) mod 5 eq 0:c in Coefficients(f)] do
        f:=ExactQuotient(f,T!5); ord+:=1;
    end while;
    return f,ord;
end function;

transformed:=[]; orders:=[];
for i in keep do
    g:=T!Evaluate(eqs[i],subs);
    gg,ord:=Strip5(g);
    Append(~transformed,gg); Append(~orders,ord);
end for;
H:=&+[lambda[i]*eqs[i]:i in [1..4]];
HH,ordH:=Strip5(T!Evaluate(H,subs));
Append(~transformed,HH); Append(~orders,ordH);

Tk<SS,OO,mm,EE,VV,NN>:=PolynomialRing(k,6,"grevlex");
red:=hom<T->Tk|SS,OO,mm,EE,VV,NN>;
geqs:=[Tk!red(g):g in transformed];

DB1,ordDB:=Strip5(T!Evaluate(DB,subs));
DC1,ordDC:=Strip5(T!Evaluate(DC,subs));
db1:=Tk!red(DB1); dc1:=Tk!red(DC1);
open_normal:=OO*mm*(VV+2*NN)*db1*dc1;
deriv_new:=[[Derivative(g,j):j in [1..6]]:g in geqs];

total:=0; open_count:=0; rank_counts:=AssociativeArray();
smooth_open:=0; samples:=[];
for s1 in k do for o1 in k do for m1 in k do
for e1 in k do for v1 in k do for n1 in k do
    pt:=[s1,o1,m1,e1,v1,n1];
    if not &and[Evaluate(g,pt) eq 0:g in geqs] then continue; end if;
    total+:=1;
    if Evaluate(open_normal,pt) eq 0 then continue; end if;
    open_count+:=1;
    JJ:=Matrix(k,4,6,
        &cat[[Evaluate(deriv_new[i][j],pt):j in [1..6]]
              :i in [1..4]]);
    rk:=Rank(JJ);
    if IsDefined(rank_counts,rk) then rank_counts[rk]+:=1;
    else rank_counts[rk]:=1; end if;
    if rk eq 4 then
        smooth_open+:=1;
        if #samples lt 24 then Append(~samples,<Z!s1,Z!o1,Z!m1,
                                       Z!e1,Z!v1,Z!n1>); end if;
    end if;
end for; end for; end for; end for; end for; end for;

test_y:=[0,1,1,3,3,4];
test_x25:=[x0[i]+5*test_y[i]:i in [1..6]];

print "CONTACT6_M612_T0_LOCAL5_FIRST_BLOWUP";
print "CENTER_s_omega_m_L_U_nu",x0;
print "ORIGINAL_JACOBIAN_MOD5",J,"RANK",Rank(J);
print "LEFT_KERNEL_ROW",lambda,"REPLACED_EQUATION",replace,
      "KEPT_EQUATIONS",keep;
print "TRANSFORMED_5_ORDERS",orders;
print "DB_DC_5_ORDERS",ordDB,ordDC;
print "TRANSFORMED_SHAPES",[<TotalDegree(g),#Terms(g)>:g in geqs];
print "TRANSFORMED_EQUATIONS_MOD5",geqs;
print "NORMALIZED_DB_DC_MOD5",db1,dc1;
print "TEST_Y",test_y,"TEST_X_MOD25",[z mod 25:z in test_x25];
print "TEST_ORIGINAL_EQUATIONS_MOD25",
      [(Z!Evaluate(f,test_x25)) mod 25:f in eqs];
print "TEST_TRANSFORMED_EQUATIONS_MOD5",
      [Z!Evaluate(g,[k!z:z in test_y]):g in geqs];
print "EXCEPTIONAL_POINTS",total,"OPEN_POINTS",open_count,
      "OPEN_RANK_COUNTS",Sort([<r,rank_counts[r]>:r in Keys(rank_counts)]),
      "SMOOTH_OPEN",smooth_open;
print "SMOOTH_OPEN_SAMPLES_S_O_mu_E_V_N",samples;
print "CONTACT6_M612_T0_LOCAL5_FIRST_BLOWUP_DONE";
quit;
