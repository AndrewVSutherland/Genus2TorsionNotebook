//////////////////////////////////////////////////////////////////////
// Compensated internal omega-pole chart for the T0 pullback at p=5.
//
// For t>=1 put
//   s=5^(2t)*S, omega=O/5^t, a+3=S*O^2,
// and write m^2-c4=5^t*D.  The reduced halving equation becomes
//
// (8*S*(c3-4*5^t*S*O)-D^2)^2
//       = 256*S^2*m^2*(c2+O*D).
//
// This is the generic compensated pole cone v(s)=-2v(omega), in which
// a+3 stays a unit.  Its mod-5 fiber is independent of t>=1.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Z:=Integers(); k:=GF(5);
if not assigned pole_order then pole_order:=1;
elif Type(pole_order) eq MonStgElt then
    pole_order:=StringToInteger(pole_order);
end if;
require pole_order ge 1: "pole_order must be positive";
pt:=5^pole_order; p2t:=pt^2;

R<S,O,m,D,L,U,nu>:=PolynomialRing(Z,7,"grevlex");
a:=S*O^2-3; b:=p2t*S-3; M:=L^2;
c5:=2*p2t*S; c4:=b^2+2*a-15; c3:=2*a*b+22;
c2:=a^2+2*b-15; c1:=2*(a+3);
G:=m^2-c4-pt*D;
H:=(8*S*(c3-4*pt*S*O)-D^2)^2
   -256*S^2*m^2*(c2+O*D);
B3:=c5*M+3*U;
Delta3:=4*c4*M+12*(U^2+nu^2)-B3^2;
F3:=B3*Delta3+16*nu^3-8*c3*M-8*U^3-48*U*nu^2;
F2:=Delta3^2+64*B3*nu^3-64*c2*M
    -192*(U^2*nu^2+nu^4);
F1:=Delta3*nu^3-4*c1*M-12*U*nu^4;
eqs:=[G,H,F1,F2,F3];
derivs:=[[Derivative(f,j):j in [1..7]]:f in eqs];

total:=0; base_units:=0; m_units:=0; contact_open:=0;
rank_counts:=AssociativeArray(); rank5:=0; rank5_contact_open:=0;
rank5_samples:=[]; deep_counts:=AssociativeArray();
for s0 in k do for o0 in k do for m0 in k do for d0 in k do
for l0 in k do for u0 in k do for n0 in k do
    z:=[s0,o0,m0,d0,l0,u0,n0];
    if not &and[(Z!Evaluate(f,z)) mod 5 eq 0:f in eqs] then continue; end if;
    total+:=1;
    if s0*o0 eq 0 then continue; end if;
    base_units+:=1;
    if m0 ne 0 then m_units+:=1; end if;
    cop:=m0*l0*n0*(u0^2-4*n0^2) ne 0;
    if cop then contact_open+:=1; end if;
    A:=Matrix(k,5,7,
        &cat[[k!((Z!Evaluate(derivs[i][j],z)) mod 5):j in [1..7]]
              :i in [1..5]]);
    rk:=Rank(A);
    if IsDefined(rank_counts,rk) then rank_counts[rk]+:=1;
    else rank_counts[rk]:=1; end if;
    if rk eq 5 then
        rank5+:=1;
        if cop then rank5_contact_open+:=1; end if;
        if #rank5_samples lt 32 then
            Append(~rank5_samples,<Z!s0,Z!o0,Z!m0,Z!d0,
                                   Z!l0,Z!u0,Z!n0,cop>);
        end if;
    end if;
    names:=["m","L","nu","Up","Um"];
    vals:=[m0,l0,n0,u0-2*n0,u0+2*n0];
    for i in [1..5] do if vals[i] eq 0 then
        if IsDefined(deep_counts,names[i]) then deep_counts[names[i]]+:=1;
        else deep_counts[names[i]]:=1; end if;
    end if; end for;
end for; end for; end for; end for; end for; end for; end for;

print "CONTACT6_M612_T0_LOCAL5_INTERNAL_OMEGA";
print "POLE_ORDER",pole_order,"EQUATION_SHAPES",
      [<TotalDegree(f),#Terms(f)>:f in eqs];
print "TOTAL",total,"S_O_UNITS",base_units,"M_UNITS",m_units,
      "FULL_CONTACT_OPEN",contact_open;
print "RANK_COUNTS_ON_S_O_UNITS",
      Sort([<r,rank_counts[r]>:r in Keys(rank_counts)]);
print "RANK5",rank5,"RANK5_CONTACT_OPEN",rank5_contact_open;
print "DEEP_VANISH_COUNTS",Sort([<x,deep_counts[x]>:x in Keys(deep_counts)]);
print "RANK5_SAMPLES_S_O_m_D_L_U_nu_open",rank5_samples;
print "CONTACT6_M612_T0_LOCAL5_INTERNAL_OMEGA_DONE";
quit;
