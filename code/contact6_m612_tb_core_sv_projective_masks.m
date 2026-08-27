//////////////////////////////////////////////////////////////////////
// Projective finite masks on the fixed slice parameters (s,v) for the
// T_B-halving + [3,3] fiber product, with K=m^2 and M=L^2 built in.
//
// For every (s,v) over F_p we compare:
//   raw_projective: projective closure of the cleared equations;
//   sat_projective: first saturate by the generic open product, then close.
//
// The raw mask is conservative for saturation-boundary reductions.  The
// saturated mask is much sharper; their difference is printed explicitly
// and must be treated as boundary, not discarded silently.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Z:=Integers();
if assigned primes then
    if Type(primes) eq MonStgElt then
        prime_list:=[StringToInteger(t):t in Split(primes,",")|#t gt 0];
    else prime_list:=primes; end if;
else prime_list:=[5,7,11]; end if;

function SliceProjectiveStatus(F,s0,v0)
    A<a,r,m,L,U>:=PolynomialRing(F,5,"grevlex");
    K:=m^2; M:=L^2;
    b0:=2*s0^2-3;
    A3:=a-3+4*s0^2*r-2*K;
    H1:=s0*((a-3)*r^2+4*r-K*(a+3))-r*A3;
    H2:=8*s0^2*(2*s0^2*r^2+2*(a-3)*r+2-K*(2*s0^2-6))
        -A3^2-32*s0^3*r;
    D:=Coefficient(H1,a,1);
    N:=-Evaluate(H1,[A!0,r,m,L,U]);
    if D eq 0 then return true,true,-2,-2,"D_identically_zero"; end if;
    function SubA(poly)
        d:=Degree(poly,a);
        return &+[Coefficient(poly,a,i)*N^i*D^(d-i):i in [0..d]];
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

    R<rr,mm,LL,UU>:=PolynomialRing(F,4,"grevlex");
    toR:=hom<A -> R | R!0,rr,mm,LL,UU>;
    G0:=R!toR(HTB); G1:=R!toR(SubA(F1));
    G2:=R!toR(SubA(F2)); G3:=R!toR(SubA(F3));
    D4:=R!toR(D); N4:=R!toR(N);
    Iraw:=ideal<R|G0,G1,G2,G3>;

    raw_has:=true; raw_dim:=-2;
    try
        Xraw:=Scheme(AffineSpace(R),Basis(Iraw));
        Praw:=ProjectiveClosure(Xraw);
        raw_dim:=Dimension(Praw);
        raw_has:=#Points(Praw) gt 0;
    catch e
        raw_has:=true;
    end try;

    boundary:=rr*mm*LL*D4*(UU^2-4*v0^2)*(N4+(b0+2)*D4);
    if boundary eq 0 then return raw_has,true,raw_dim,-2,"boundary_zero"; end if;
    Isat:=Iraw; sat_has:=true; sat_dim:=-2;
    try
        Isat:=Saturation(Isat,ideal<R|boundary>);
        Xsat:=Scheme(AffineSpace(R),Basis(Isat));
        Psat:=ProjectiveClosure(Xsat);
        sat_dim:=Dimension(Psat);
        sat_has:=#Points(Psat) gt 0;
    catch e
        sat_has:=true;
    end try;
    return raw_has,sat_has,raw_dim,sat_dim,"ok";
end function;

print "CONTACT6_M612_TB_CORE_SV_PROJECTIVE_MASKS";
for p in prime_list do
    F:=GF(p);
    raw_allowed:=[]; sat_allowed:=[]; boundary_only:=[]; failures:=[];
    dim_stats:=AssociativeArray();
    t0:=Cputime();
    for s0 in F do
        for v0 in F do
            raw,sat,dr,ds,status:=SliceProjectiveStatus(F,s0,v0);
            key:=<Z!s0,Z!v0>;
            if raw then Append(~raw_allowed,key); end if;
            if sat then Append(~sat_allowed,key); end if;
            if raw and not sat then Append(~boundary_only,key); end if;
            if status ne "ok" then Append(~failures,<key,status>); end if;
            dkey:=<dr,ds>;
            if IsDefined(dim_stats,dkey) then dim_stats[dkey]+:=1;
            else dim_stats[dkey]:=1; end if;
        end for;
    end for;
    print "p",p,"raw_allowed",#raw_allowed,"sat_allowed",#sat_allowed,
          "boundary_only",#boundary_only,"failures",#failures,
          "seconds",Cputime(t0);
    print "RAW",raw_allowed;
    print "SAT",sat_allowed;
    print "BOUNDARY_ONLY",boundary_only;
    print "FAILURES",failures;
    print "DIM_STATS",Sort([<k,dim_stats[k]>:k in Keys(dim_stats)]);
end for;
quit;
