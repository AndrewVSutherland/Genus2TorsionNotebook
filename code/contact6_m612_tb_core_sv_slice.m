//////////////////////////////////////////////////////////////////////
// Discovery-oriented algebraic slice of the full T_B-halving + [3,3]
// fiber product.
//
// Fix rational (s,v).  The T_B chart has b=2*s^2-3 and recovers a
// linearly from (r,K).  Substitute that recovery into the three cubic-
// contact equations.  The resulting four equations are in (r,K,M,U):
// one T_B cubic and three contact equations.  This removes any height
// bound on r and is the natural complement to the fast fixed-(s,r)
// verifier.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned s_num then s_num:=1; end if;
if not assigned s_den then s_den:=1; end if;
if not assigned v_num then v_num:=-1; end if;
if not assigned v_den then v_den:=1; end if;
if Type(s_num) eq MonStgElt then s_num:=StringToInteger(s_num); end if;
if Type(s_den) eq MonStgElt then s_den:=StringToInteger(s_den); end if;
if Type(v_num) eq MonStgElt then v_num:=StringToInteger(v_num); end if;
if Type(v_den) eq MonStgElt then v_den:=StringToInteger(v_den); end if;

Q:=Rationals(); Z:=Integers();
s0:=Q!s_num/s_den; v0:=Q!v_num/v_den;
if s0 eq 0 or v0 eq 0 then error "require s*v != 0"; end if;

function PrimitiveA(f)
    if f eq 0 then return f; end if;
    den:=LCM([Denominator(c):c in Coefficients(f)]);
    g:=Parent(f)!(den*f);
    nums:=[Z!c:c in Coefficients(g)];
    cont:=GCD([Abs(n):n in nums|n ne 0]);
    return cont gt 1 select Parent(f)!(g/cont) else g;
end function;

A<a,r,K,M,U>:=PolynomialRing(Q,5,"grevlex");
b0:=2*s0^2-3;
A3:=a-3+4*s0^2*r-2*K;
H1:=s0*((a-3)*r^2+4*r-K*(a+3))-r*A3;
H2:=8*s0^2*(2*s0^2*r^2+2*(a-3)*r+2-K*(2*s0^2-6))
    -A3^2-32*s0^3*r;
D:=Coefficient(H1,a,1);
N:=-Evaluate(H1,[A!0,r,K,M,U]);
assert H1 eq D*a-N;

function SubA(P)
    d:=Degree(P,a);
    return &+[Coefficient(P,a,i)*N^i*D^(d-i):i in [0..d]];
end function;

HTB:=SubA(H2);
if IsDivisibleBy(HTB,K) then HTB:=ExactQuotient(HTB,K); end if;

c1:=2*a+6;
c2:=a^2+2*b0-15;
c3:=2*a*b0+22;
c4:=2*a+b0^2-15;
c5:=2*b0+6;
B3:=c5*M+3*U;
Delta3:=4*c4*M+12*(U^2+v0^2)-B3^2;
F3:=B3*Delta3+16*v0^3-8*c3*M-8*U^3-48*U*v0^2;
F2:=Delta3^2+64*B3*v0^3-64*c2*M-192*(U^2*v0^2+v0^4);
F1:=Delta3*v0^3-4*c1*M-12*U*v0^4;

R<rr,KK,MM,UU>:=PolynomialRing(Q,4,"grevlex");
toR:=hom<A -> R | R!0,rr,KK,MM,UU>;
G0:=PrimitiveA(R!toR(HTB));
G1:=PrimitiveA(R!toR(SubA(F1)));
G2:=PrimitiveA(R!toR(SubA(F2)));
G3:=PrimitiveA(R!toR(SubA(F3)));
D4:=R!toR(D);
N4:=R!toR(N);

print "CONTACT6_M612_TB_CORE_SV_SLICE";
print "s",s0,"v",v0,"b",b0;
for rec in [<"TB",G0>,<"F1",G1>,<"F2",G2>,<"F3",G3>] do
    print rec[1],"total_degree",TotalDegree(rec[2]),
          "degrees",[Degree(rec[2],i):i in [1..4]],
          "terms",#Terms(rec[2]);
end for;

t0:=Cputime();
I:=ideal<R|G0,G1,G2,G3>;
dim0:=Dimension(I);
print "raw_dimension",dim0,"seconds",Cputime(t0);

boundary:=rr*KK*MM*D4*(UU^2-4*v0^2)*(N4+(b0+2)*D4);
sat_ok:=true;
t1:=Cputime();
try I:=Saturation(I,ideal<R|boundary>);
catch e sat_ok:=false; end try;
dim:=Dimension(I);
print "saturated",sat_ok,"dimension",dim,"seconds",Cputime(t1);

sols:=[]; variety_ok:=false;
if dim eq 0 then
    t2:=Cputime();
    try sols:=Variety(I); variety_ok:=true; catch e variety_ok:=false; end try;
    print "variety_ok",variety_ok,"rational_points",#sols,
          "seconds",Cputime(t2);
    for pt in sols do
        r0:=Q!pt[1]; K0:=Q!pt[2]; M0:=Q!pt[3]; U0:=Q!pt[4];
        d0:=Q!Evaluate(D4,<r0,K0,M0,U0>);
        if d0 eq 0 then continue; end if;
        a0:=Q!Evaluate(N4,<r0,K0,M0,U0>)/d0;
        sqK,m0:=IsSquare(K0); sqM,L0:=IsSquare(M0);
        print "POINT","r",r0,"K",K0,"M",M0,"U",U0,"a",a0,
              "Ksquare",sqK,"Msquare",sqM;
    end for;
end if;

print "DONE","total_seconds",Cputime(t0);
quit;
