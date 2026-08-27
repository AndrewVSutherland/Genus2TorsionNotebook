//////////////////////////////////////////////////////////////////////
// Integral mod-5 diagnostic for the T0-halved contact-6 [6,12] cover.
//
// Variables are (s,omega,m,L,U,nu), with
//   b+3=s, a+3=s*omega^2, W=m^2, M=L^2.
// This enumerates the complete special fiber, records the forced base and
// deeper contact strata, and tests smoothness of the four-equation cover.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
k:=GF(5); Z:=Integers();
R<s,omega,m,L,U,nu>:=PolynomialRing(k,6,"grevlex");

a:=s*omega^2-3; b:=s-3; W:=m^2; M:=L^2;
c5:=2*s;
c4:=b^2+2*a-15;
c3:=2*a*b+22;
c2:=a^2+2*b-15;
c1:=2*(a+3);
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

DB:=(a-3)^2-8*(b+3);
DC:=(b-3)^2-8*(a+3);
RRred:=a+b+2;
base_names:=["b3","a3","RR","DB","DC"];
base_forms:=[s,s*omega^2,RRred,DB,DC];
deep_names:=["m","L","nu","Up","Um","omega_pm_1"];
deep_forms:=[m,L,nu,U-2*nu,U+2*nu,omega^2-1];

derivs:=[[Derivative(f,j):j in [1..6]]:f in eqs];
base_counts:=AssociativeArray(); deep_counts:=AssociativeArray();
pair_counts:=AssociativeArray(); rank_counts:=AssociativeArray();
total:=0; forced:=0; cover_open:=0; smooth:=0; smooth_open:=0;
smooth_can_escape_all:=0; smooth_can_escape_deep:=0;
samples:=[];
for s0 in k do for o0 in k do for m0 in k do
for l0 in k do for u0 in k do for n0 in k do
    pt:=[s0,o0,m0,l0,u0,n0];
    if not &and[Evaluate(f,pt) eq 0:f in eqs] then continue; end if;
    total+:=1;
    bset:=[i:i in [1..#base_forms]|Evaluate(base_forms[i],pt) eq 0];
    dset:=[i:i in [1..#deep_forms]|Evaluate(deep_forms[i],pt) eq 0];
    if #bset gt 0 then forced+:=1; end if;
    isopen:=s0*m0*l0*n0*(o0^2-1)*(u0^2-4*n0^2)
            *Evaluate(RRred,pt) ne 0;
    if isopen then cover_open+:=1; end if;
    J:=Matrix(k,4,6,&cat[[Evaluate(derivs[i][j],pt):j in [1..6]]
                          :i in [1..4]]);
    rk:=Rank(J);
    if IsDefined(rank_counts,rk) then rank_counts[rk]+:=1;
    else rank_counts[rk]:=1; end if;
    if rk eq 4 then
        smooth+:=1;
        if isopen then smooth_open+:=1; end if;
        // A vanishing factor can be left to first order precisely when its
        // gradient is nonzero on the tangent space, equivalently adjoining
        // its gradient increases the Jacobian rank.
        can_escape_deep:=true;
        for j in dset do
            gd:=Matrix(k,1,6,[Evaluate(Derivative(deep_forms[j],h),pt)
                              :h in [1..6]]);
            if Rank(VerticalJoin(J,gd)) eq rk then
                can_escape_deep:=false;
            end if;
        end for;
        if can_escape_deep then smooth_can_escape_deep+:=1; end if;
        can_escape_base:=true;
        for i in bset do
            gb:=Matrix(k,1,6,[Evaluate(Derivative(base_forms[i],h),pt)
                              :h in [1..6]]);
            if Rank(VerticalJoin(J,gb)) eq rk then
                can_escape_base:=false;
            end if;
        end for;
        if can_escape_deep and can_escape_base then
            smooth_can_escape_all+:=1;
        end if;
        if #samples lt 24 then
            Append(~samples,Sprint(<[Z!x:x in pt],
                [base_names[i]:i in bset],[deep_names[i]:i in dset],
                can_escape_base,can_escape_deep>));
        end if;
    end if;
    for i in bset do
        bn:=base_names[i];
        if IsDefined(base_counts,bn) then base_counts[bn]+:=1;
        else base_counts[bn]:=1; end if;
        for j in dset do
            key:=bn cat ":" cat deep_names[j];
            if IsDefined(pair_counts,key) then pair_counts[key]+:=1;
            else pair_counts[key]:=1; end if;
        end for;
    end for;
    for j in dset do
        dn:=deep_names[j];
        if IsDefined(deep_counts,dn) then deep_counts[dn]+:=1;
        else deep_counts[dn]:=1; end if;
    end for;
end for; end for; end for; end for; end for; end for;

print "CONTACT6_M612_T0_LOCAL5_INTEGRAL";
print "EQUATION_SHAPES",[<TotalDegree(f),#Terms(f)>:f in eqs];
print "TOTAL",total,"ON_FORCED_BASE",forced,"FULL_COVER_OPEN",cover_open;
print "JACOBIAN_RANK_COUNTS",Sort([<x,rank_counts[x]>:x in Keys(rank_counts)]);
print "SMOOTH",smooth,"SMOOTH_FULL_OPEN",smooth_open;
print "SMOOTH_CAN_ESCAPE_DEEP_FIRST_ORDER",smooth_can_escape_deep,
      "SMOOTH_CAN_ESCAPE_ALL_VANISHING_FIRST_ORDER",smooth_can_escape_all;
print "BASE_COUNTS",Sort([<x,base_counts[x]>:x in Keys(base_counts)]);
print "DEEP_COUNTS",Sort([<x,deep_counts[x]>:x in Keys(deep_counts)]);
print "BASE_DEEP_INCIDENCES",Sort([<x,pair_counts[x]>:x in Keys(pair_counts)]);
print "SMOOTH_SAMPLES_s_omega_m_L_U_nu__base__deep__escape_base__escape_deep",
      samples;
print "CONTACT6_M612_T0_LOCAL5_INTEGRAL_DONE";
quit;
