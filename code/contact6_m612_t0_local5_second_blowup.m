//////////////////////////////////////////////////////////////////////
// Second deflated blow-up of the rank-3 T0 p=5 cone.
//
// The first exceptional open fiber is
//   O=r, mu=q, E=3r, S=2r+3, V=1-3r, N=-r,
// with r*q != 0.  For all 16 (r,q), substitute y=y0+5*z into
// the exact first-deflated equations, perform the next invertible row
// deflation, and enumerate the mod-5 second exceptional fiber.  A point in
// this fiber would solve the first-deflated system modulo 25; in original
// normalization this means F1,F2,F3 modulo 125 and Fhalf modulo 625.
// These unequal precisions are intentional: every exact Q_5 point must
// satisfy them.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Z:=Integers(); k:=GF(5);
R<s,omega,m,L,U,nu>:=PolynomialRing(Z,6,"grevlex");

a:=s*omega^2-3; b:=s-3; W:=m^2; M:=L^2;
c5:=2*s; c4:=b^2+2*a-15; c3:=2*a*b+22;
c2:=a^2+2*b-15; c1:=2*(a+3);
N0:=c2-omega*c4+omega*W;
R8:=8*s*(c3-4*s*omega)-(c4-W)^2;
Fhalf:=R8^2-256*s^2*W*N0;
B3:=c5*M+3*U;
Delta3:=4*c4*M+12*(U^2+nu^2)-B3^2;
F3:=B3*Delta3+16*nu^3-8*c3*M-8*U^3-48*U*nu^2;
F2:=Delta3^2+64*B3*nu^3-64*c2*M
    -192*(U^2*nu^2+nu^4);
F1:=Delta3*nu^3-4*c1*M-12*U*nu^4;
original:=[Fhalf,F1,F2,F3];

T<S,O,mu,E,V,Nn>:=PolynomialRing(Z,6,"grevlex");
xcenter:=[3,1,0,2,2,-1];
subs1:=[3+5*S,1+5*O,5*mu,2+5*E,2+5*V,-1+5*Nn];

function StripInT(f)
    ord:=0;
    while f ne 0 and &and[(Z!c) mod 5 eq 0:c in Coefficients(f)] do
        f:=ExactQuotient(f,T!5); ord+:=1;
    end while;
    return f,ord;
end function;

// The original first row has zero differential at the center; it is the
// left-kernel equation and receives the extra division by 5.
first:=[];
for i in [2,3,4] do
    g,ord:=StripInT(T!Evaluate(original[i],subs1));
    assert ord eq 1; Append(~first,g);
end for;
g,ord:=StripInT(T!Evaluate(original[1],subs1));
assert ord eq 2; Append(~first,g);

Zr<A,B,C,D,G,H>:=PolynomialRing(Z,6,"grevlex");
Kr<AA,BB,CC,DD,GG,HH>:=PolynomialRing(k,6,"grevlex");
redZ:=hom<Zr->Kr|AA,BB,CC,DD,GG,HH>;
Tr<SS,OO,MM,EE,VV,NN>:=PolynomialRing(k,6,"grevlex");
redT:=hom<T->Tr|SS,OO,MM,EE,VV,NN>;
firstmod:=[Tr!redT(f):f in first];
dfirst:=[[Derivative(f,j):j in [1..6]]:f in firstmod];

function StripInZr(f)
    ord:=0;
    while f ne 0 and &and[(Z!c) mod 5 eq 0:c in Coefficients(f)] do
        f:=ExactQuotient(f,Zr!5); ord+:=1;
    end while;
    return f,ord;
end function;

print "CONTACT6_M612_T0_LOCAL5_SECOND_BLOWUP";
grand_points:=0; grand_rank4:=0; branch_data:=[]; smooth_samples:=[];
for r in [1..4] do for q in [1..4] do
    y0:=[(2*r+3) mod 5,r,q,(3*r) mod 5,(1-3*r) mod 5,(-r) mod 5];
    y0k:=[k!z:z in y0];
    assert &and[Evaluate(f,y0k) eq 0:f in firstmod];
    J1:=Matrix(k,4,6,
        &cat[[Evaluate(dfirst[i][j],y0k):j in [1..6]]:i in [1..4]]);
    assert Rank(J1) eq 3;
    Kl:=Kernel(J1); assert Dimension(Kl) eq 1;
    lamF:=Basis(Kl)[1]; lam:=[Z!lamF[i]:i in [1..4]];
    replace:=[i:i in [1..4]|lam[i] mod 5 ne 0][1];
    keep:=[i:i in [1..4]|i ne replace];
    subs2:=[y0[1]+5*A,y0[2]+5*B,y0[3]+5*C,
            y0[4]+5*D,y0[5]+5*G,y0[6]+5*H];
    second:=[]; orders:=[];
    for i in keep do
        ff,oo:=StripInZr(Zr!Evaluate(first[i],subs2));
        Append(~second,ff); Append(~orders,oo);
    end for;
    combo:=&+[lam[i]*first[i]:i in [1..4]];
    ff,oo:=StripInZr(Zr!Evaluate(combo,subs2));
    Append(~second,ff); Append(~orders,oo);
    secondmod:=[Kr!redZ(f):f in second];
    dsecond:=[[Derivative(f,j):j in [1..6]]:f in secondmod];
    count:=0; rank4:=0; ranks:=AssociativeArray();
    for a1 in k do for b1 in k do for c1v in k do
    for d1 in k do for g1 in k do for h1 in k do
        zpt:=[a1,b1,c1v,d1,g1,h1];
        if not &and[Evaluate(f,zpt) eq 0:f in secondmod] then continue; end if;
        count+:=1;
        JJ:=Matrix(k,4,6,
            &cat[[Evaluate(dsecond[i][j],zpt):j in [1..6]]
                  :i in [1..4]]);
        rk:=Rank(JJ);
        if IsDefined(ranks,rk) then ranks[rk]+:=1; else ranks[rk]:=1; end if;
        if rk eq 4 then
            rank4+:=1;
            if #smooth_samples lt 32 then
                zZ:=[Z!z:z in zpt];
                y25:=[y0[i]+5*zZ[i]:i in [1..6]];
                x125:=[xcenter[i]+5*y25[i]:i in [1..6]];
                Append(~smooth_samples,<r,q,zZ,[z mod 125:z in x125]>);
            end if;
        end if;
    end for; end for; end for; end for; end for; end for;
    grand_points+:=count; grand_rank4+:=rank4;
    Append(~branch_data,<r,q,lam,replace,orders,count,
                         Sort([<rk,ranks[rk]>:rk in Keys(ranks)]),rank4>);
end for; end for;

print "BRANCH_DATA_r_q_lambda_replaced_orders_points_ranks_rank4";
for rec in branch_data do print rec; end for;
print "TOTAL_SECOND_EXCEPTIONAL_POINTS",grand_points,
      "TOTAL_RANK4",grand_rank4;
print "SMOOTH_SAMPLES_r_q_z_xmod125",smooth_samples;
print "CONTACT6_M612_T0_LOCAL5_SECOND_BLOWUP_DONE";
quit;
