//////////////////////////////////////////////////////////////////////
// Exact mod-7 slice masks for the T0-halved contact-6 core pullback.
//
// The global solver fixes (omega,nu), so this script projects every
// relevant mod-7 point to those two coordinates.  The priority set consists
// of points on the forced base boundary for which all auxiliary T0/contact
// factors are units and the four-equation pullback has Jacobian rank four.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Z:=Integers(); k:=GF(7);
R<s,omega,m,L,U,nu>:=PolynomialRing(k,6,"grevlex");

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
RR:=a+b+2;
base_forms:=[s,s*omega^2,RR,DB,DC];
// These are the non-base factors required to stay in the genuine T0 and
// cubic-contact charts.  A base factor may vanish in the special fiber and
// is then required to be escapable transversely.
aux_open:=m*L*nu*(omega^2-1)*(U^2-4*nu^2);

total:=0; aux_open_total:=0; boundary_aux_open:=0;
boundary_rank4:=0; boundary_rank4_escape:=0;
pairs:={@ @}; full_points:=[]; by_pair:=AssociativeArray();

for s0 in k do for o0 in k do for m0 in k do
for l0 in k do for u0 in k do for n0 in k do
    pt:=[s0,o0,m0,l0,u0,n0];
    if not &and[Evaluate(f,pt) eq 0:f in eqs] then continue; end if;
    total+:=1;
    if Evaluate(aux_open,pt) eq 0 then continue; end if;
    aux_open_total+:=1;
    bzero:=[i:i in [1..#base_forms]|Evaluate(base_forms[i],pt) eq 0];
    if #bzero eq 0 then continue; end if;
    boundary_aux_open+:=1;
    J:=Matrix(k,4,6,&cat[[Evaluate(derivs[i][j],pt):j in [1..6]]
                          :i in [1..4]]);
    rk:=Rank(J);
    if rk ne 4 then continue; end if;
    boundary_rank4+:=1;
    // Each vanishing base factor must be nonconstant on the smooth local
    // surface.  Then a common Q_7 point can be chosen off their union.
    escapes:=true;
    for i in bzero do
        g:=Matrix(k,1,6,[Evaluate(Derivative(base_forms[i],j),pt)
                         :j in [1..6]]);
        if Rank(VerticalJoin(J,g)) eq 4 then escapes:=false; end if;
    end for;
    if not escapes then continue; end if;
    boundary_rank4_escape+:=1;
    pair:=<Z!o0,Z!n0>;
    Include(~pairs,pair);
    key:=Sprint(pair);
    if IsDefined(by_pair,key) then by_pair[key]+:=1; else by_pair[key]:=1; end if;
    if #full_points lt 80 then Append(~full_points,<Z!s0,Z!o0,Z!m0,
                                       Z!l0,Z!u0,Z!n0,bzero>); end if;
end for; end for; end for; end for; end for; end for;

print "CONTACT6_M612_T0_LOCAL7_MASKS";
print "TOTAL",total,"AUX_OPEN",aux_open_total,
      "BOUNDARY_AUX_OPEN",boundary_aux_open;
print "BOUNDARY_RANK4",boundary_rank4,
      "BOUNDARY_RANK4_ESCAPE_EACH_BASE_FACTOR",boundary_rank4_escape;
print "OMEGA_NU_PAIRS_MOD7",Sort(Setseq(pairs));
print "PAIR_COUNTS",Sort([<x,by_pair[x]>:x in Keys(by_pair)]);
print "SAMPLE_POINTS_s_omega_m_L_U_nu_base_indices",full_points;
print "CONTACT6_M612_T0_LOCAL7_MASKS_DONE";
quit;
