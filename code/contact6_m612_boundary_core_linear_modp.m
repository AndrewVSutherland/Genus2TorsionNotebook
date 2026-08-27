//////////////////////////////////////////////////////////////////////
// Lightweight first-order p=5,7 boundary test for the full contact-6
// cubic-contact core, with M=L^2 imposed.
//
// For every raw point z0=(t,L,U,v) on the special fiber e=0, solve
//
//   (dF/dz)(z0) z1 = -(dF/de)(z0).
//
// Solvability is an exact necessary condition for a deformation
// transverse to the boundary.  Unlike naive point enumeration, it detects
// which points survive to first order off the excess singular fiber.  It is
// intentionally much cheaper than a full Groebner saturation.  The two
// infinity charts keep L,U,v affine, as in boundary_core_modp.m.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Z:=Integers();
if assigned primes then
    prime_list:=Type(primes) eq MonStgElt
        select [StringToInteger(s):s in Split(primes,",")|#s gt 0]
        else primes;
else prime_list:=[5,7]; end if;

all_charts:=["b3","a3","rr","DB","DC","infB","infA"];
if assigned charts then
    chart_list:=Type(charts) eq MonStgElt
        select [s:s in Split(charts,",")|#s gt 0]
        else charts;
    if "all" in chart_list then chart_list:=all_charts; end if;
else chart_list:=all_charts; end if;
if not assigned geometric_saturation then
    geometric_saturation:=false;
elif Type(geometric_saturation) eq MonStgElt then
    geometric_saturation:=geometric_saturation in {"true","True","1","yes"};
end if;

function CoreEquations(a,b,L,U,v)
    M:=L^2;
    c1:=2*a+6; c2:=a^2+2*b-15; c3:=2*a*b+22;
    c4:=2*a+b^2-15; c5:=2*b+6;
    B3:=c5*M+3*U;
    Delta3:=4*c4*M+12*(U^2+v^2)-B3^2;
    F3:=B3*Delta3+16*v^3-8*c3*M-8*U^3-48*U*v^2;
    F2:=Delta3^2+64*B3*v^3-64*c2*M
        -192*(U^2*v^2+v^4);
    F1:=Delta3*v^3-4*c1*M-12*U*v^4;
    return [F1,F2,F3];
end function;

print "CONTACT6_M612_BOUNDARY_CORE_LINEAR_MODP";
print "primes",prime_list,"charts",chart_list,
      "geometric_saturation",geometric_saturation;

for p0 in prime_list do
    p:=Z!p0; require p in {5,7}: "only p=5,7 are supported";
    k:=GF(p); R<e,t,L,U,v>:=PolynomialRing(k,5,"grevlex");
    inv8:=(k!8)^-1;
    for label in chart_list do
        if label eq "b3" then
            eqs:=CoreEquations(t,-3+e,L,U,v);
            description:="a=t,b=-3+e";
        elif label eq "a3" then
            eqs:=CoreEquations(-3+e,t,L,U,v);
            description:="a=-3+e,b=t";
        elif label eq "rr" then
            eqs:=CoreEquations(t,-t-2+e,L,U,v);
            description:="a=t,b=-t-2+e";
        elif label eq "DB" then
            eqs:=CoreEquations(t,inv8*(t-3)^2-3+e,L,U,v);
            description:="a=t,b=(t-3)^2/8-3+e";
        elif label eq "DC" then
            eqs:=CoreEquations(inv8*(t-3)^2-3+e,t,L,U,v);
            description:="b=t,a=(t-3)^2/8-3+e";
        elif label in {"infB","infA"} then
            K:=FieldOfFractions(R);
            if label eq "infB" then
                aa:=(K!t)/(K!e); bb:=(K!1)/(K!e);
                description:="[A:B:T]=[t:1:e]";
            else
                aa:=(K!1)/(K!e); bb:=(K!t)/(K!e);
                description:="[A:B:T]=[1:t:e]";
            end if;
            raw:=CoreEquations(aa,bb,K!L,K!U,K!v);
            degs:=[2,4,3]; eqs:=[];
            for i in [1..3] do
                cleared:=(K!e)^degs[i]*raw[i];
                assert Denominator(cleared) eq 1;
                Append(~eqs,R!Numerator(cleared));
            end for;
        else error "unknown chart",label;
        end if;

        dnormal:=[Derivative(g,1):g in eqs];
        dz:=[[Derivative(g,j):j in [2..5]]:g in eqs];
        if geometric_saturation then
            S<tt,LL,UU,vv>:=PolynomialRing(k,4,"grevlex");
            toS:=hom<R -> S | S!0,tt,LL,UU,vv>;
            I0:=ideal<S|[S!toS(g):g in eqs]>;
            I0open:=Saturation(I0,ideal<S|LL*vv*(UU^2-4*vv^2)>);
            open_empty_geom:=S!1 in I0open;
            open_geom_dim:=open_empty_geom select -1 else Dimension(I0open);
            open_basis_size:=#Basis(I0open);
        else
            open_empty_geom:=false; open_geom_dim:=-2; open_basis_size:=-1;
        end if;
        raw:=0; openraw:=0; transverse:=0; open_transverse:=0;
        fullrank:=0; fullrank_open:=0;
        rank_counts:=AssociativeArray(); by_t:=AssociativeArray();
        samples:=[];
        for t0 in k do for l0 in k do for u0 in k do for v0 in k do
            pt:=[k!0,t0,l0,u0,v0];
            if not &and[Evaluate(g,pt) eq 0:g in eqs] then continue; end if;
            raw+:=1;
            openpt:=l0*v0*(u0^2-4*v0^2) ne 0;
            if openpt then openraw+:=1; end if;
            A:=Matrix(k,3,4,
                &cat[[Evaluate(dz[i][j],pt):j in [1..4]]:i in [1..3]]);
            rhs:=Matrix(k,3,1,[-Evaluate(dnormal[i],pt):i in [1..3]]);
            rk:=Rank(A); augrk:=Rank(HorizontalJoin(A,rhs));
            key:=<rk,augrk>;
            if IsDefined(rank_counts,key) then rank_counts[key]+:=1;
            else rank_counts[key]:=1; end if;
            if rk eq 3 then
                fullrank+:=1; if openpt then fullrank_open+:=1; end if;
            end if;
            if rk eq augrk then
                transverse+:=1;
                tk:=Z!t0;
                if IsDefined(by_t,tk) then by_t[tk]+:=1; else by_t[tk]:=1; end if;
                if openpt then
                    open_transverse+:=1;
                    if #samples lt 16 then
                        Append(~samples,<Z!t0,Z!l0,Z!u0,Z!v0,rk>);
                    end if;
                end if;
            end if;
        end for; end for; end for; end for;
        print "PRIME",p,"CHART",label,description;
        print " RAW_POINTS",raw,"CONTACT_OPEN_RAW",openraw,
              "FIRST_ORDER_TRANSVERSE",transverse,
              "CONTACT_OPEN_TRANSVERSE",open_transverse;
        if geometric_saturation then
            print " GEOMETRIC_CONTACT_OPEN_EMPTY",open_empty_geom,
                  "GEOMETRIC_CONTACT_OPEN_DIM",open_geom_dim,
                  "OPEN_SAT_BASIS_SIZE",open_basis_size;
        else
            print " GEOMETRIC_CONTACT_OPEN_SATURATION SKIPPED";
        end if;
        print " FULL_FIBER_JACOBIAN_RANK",fullrank,
              "FULL_RANK_CONTACT_OPEN",fullrank_open;
        print " RANK_AUGRANK_COUNTS",
              Sort([<key,rank_counts[key]>:key in Keys(rank_counts)]);
        print " TRANSVERSE_COUNTS_BY_T",
              Sort([<key,by_t[key]>:key in Keys(by_t)]);
        print " OPEN_TRANSVERSE_SAMPLES_t_L_U_v_rank",samples;
    end for;
end for;

print "CONTACT6_M612_BOUNDARY_CORE_LINEAR_MODP_DONE";
quit;
