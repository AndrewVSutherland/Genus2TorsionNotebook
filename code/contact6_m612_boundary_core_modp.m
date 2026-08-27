//////////////////////////////////////////////////////////////////////
// Flat-limit p=5,7 boundary charts of the full contact-6 cubic-contact
// [6,6] core (including the square M=L^2).
//
// A naive substitution of a discriminant divisor into F1,F2,F3 creates
// excess-dimensional components.  For each local equation e=0 below we
// instead compute
//
//       ( <F1,F2,F3> : e^infinity ) + <e>.
//
// Thus every reported special fiber is the closure of the e != 0 core,
// not the naive singular fiber.  The two infinity charts use the base
// projective coordinates [A:B:T] and first clear the exact base degrees
// (2,4,3) of (F1,F2,F3).
//
// Examples:
//   magma -b primes:=5 charts:=b3,a3,rr code/contact6_m612_boundary_core_modp.m
//   magma -b primes:=7 charts:=all code/contact6_m612_boundary_core_modp.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Z:=Integers();

if assigned primes then
    if Type(primes) eq MonStgElt then
        prime_list:=[StringToInteger(s):s in Split(primes,",")|#s gt 0];
    else prime_list:=primes; end if;
else
    prime_list:=[5,7];
end if;
if not assigned sat_mode then sat_mode:="contact"; end if;

all_charts:=["b3","a3","rr","DB","DC","infB","infA"];
if assigned charts then
    chart_list:=Type(charts) eq MonStgElt
        select [s:s in Split(charts,",")|#s gt 0]
        else charts;
    if "all" in chart_list then chart_list:=all_charts; end if;
else
    chart_list:=all_charts;
end if;

function CoreEquations(a,b,L,U,v)
    M:=L^2;
    c1:=2*a+6;
    c2:=a^2+2*b-15;
    c3:=2*a*b+22;
    c4:=2*a+b^2-15;
    c5:=2*b+6;
    B3:=c5*M+3*U;
    Delta3:=4*c4*M+12*(U^2+v^2)-B3^2;
    F3:=B3*Delta3+16*v^3-8*c3*M-8*U^3-48*U*v^2;
    F2:=Delta3^2+64*B3*v^3-64*c2*M
        -192*(U^2*v^2+v^4);
    F1:=Delta3*v^3-4*c1*M-12*U*v^4;
    return [F1,F2,F3];
end function;

print "CONTACT6_M612_BOUNDARY_CORE_MODP";
print "primes",prime_list,"charts",chart_list,"sat_mode",sat_mode;
print "chart variables are (e,t,L,U,v); e is the normal/base-infinity parameter";

for p0 in prime_list do
    p:=Z!p0;
    require p in {5,7}: "this diagnostic is intended for p=5,7";
    k:=GF(p);
    R<e,t,L,U,v>:=PolynomialRing(k,5,"grevlex");
    inv8:=(k!8)^-1;
    for label in chart_list do
        is_infinity:=label in {"infB","infA"};
        if label eq "b3" then
            // a=t, b=-3+e; normal b+3=e.
            eqs:=CoreEquations(t,-3+e,L,U,v);
            aa:=t; bb:=-3+e;
            description:="a=t, b=-3+e";
        elif label eq "a3" then
            // b=t, a=-3+e; normal a+3=e.
            eqs:=CoreEquations(-3+e,t,L,U,v);
            aa:=-3+e; bb:=t;
            description:="a=-3+e, b=t";
        elif label eq "rr" then
            // RR=2(a+b+2)^3, so its reduced support is this line.
            eqs:=CoreEquations(t,-t-2+e,L,U,v);
            aa:=t; bb:=-t-2+e;
            description:="a=t, b=-t-2+e (reduced RR support)";
        elif label eq "DB" then
            // DB=-8e in this transverse chart.
            eqs:=CoreEquations(t,inv8*(t-3)^2-3+e,L,U,v);
            aa:=t; bb:=inv8*(t-3)^2-3+e;
            description:="a=t, b=(t-3)^2/8-3+e";
        elif label eq "DC" then
            // DC=-8e in this transverse chart.
            eqs:=CoreEquations(inv8*(t-3)^2-3+e,t,L,U,v);
            aa:=inv8*(t-3)^2-3+e; bb:=t;
            description:="b=t, a=(t-3)^2/8-3+e";
        elif label eq "infB" then
            // B=1 chart on [A:B:T]: a=t/e, b=1/e.
            K:=FieldOfFractions(R);
            raw:=CoreEquations((K!t)/(K!e),(K!1)/(K!e),
                               K!L,K!U,K!v);
            degs:=[2,4,3];
            eqs:=[];
            for i in [1..3] do
                cleared:=(K!e)^degs[i]*raw[i];
                assert Denominator(cleared) eq 1;
                Append(~eqs,R!Numerator(cleared));
            end for;
            description:="[A:B:T]=[t:1:e]";
            aa:=(K!t)/(K!e); bb:=(K!1)/(K!e);
        elif label eq "infA" then
            // A=1 chart: a=1/e, b=t/e.  Its t=0 point is [1:0:0],
            // the only point missed by infB.
            K:=FieldOfFractions(R);
            raw:=CoreEquations((K!1)/(K!e),(K!t)/(K!e),
                               K!L,K!U,K!v);
            degs:=[2,4,3];
            eqs:=[];
            for i in [1..3] do
                cleared:=(K!e)^degs[i]*raw[i];
                assert Denominator(cleared) eq 1;
                Append(~eqs,R!Numerator(cleared));
            end for;
            description:="[A:B:T]=[1:t:e]";
            aa:=(K!1)/(K!e); bb:=(K!t)/(K!e);
        else
            error "unknown chart",label;
        end if;

        // Close the genuine core open, not the large degenerate contact
        // components (L=0, v=0, or disc(q)=0), and not a component trapped
        // in another base-discriminant support.  Saturating by the product
        // still retains intersections approached from the simultaneous open.
        DBexpr:=(aa-3)^2-8*(bb+3);
        DCexpr:=(bb-3)^2-8*(aa+3);
        Iraw:=ideal<R|eqs>;
        raw_special:=Iraw+ideal<R|e>;
        raw_dim:=Dimension(raw_special);
        tstart:=Cputime();
        // Sequential saturation is much smaller than one saturation by the
        // high-degree product.  The default closes the nondegenerate
        // cubic-contact chart transversely to e=0 while retaining crossings
        // with the other base supports.  sat_mode:="full" also removes
        // components trapped in any other base-boundary support.
        sat_factors:=[e,L,v,U^2-4*v^2];
        if sat_mode eq "full" then
            sat_factors cat:= [R!Numerator(aa+3),R!Numerator(bb+3),
                               R!Numerator(aa+bb+2),
                               R!Numerator(DBexpr),R!Numerator(DCexpr)];
        end if;
        Isat:=Iraw;
        for sf in sat_factors do
            if sf ne 0 then Isat:=Saturation(Isat,ideal<R|sf>); end if;
        end for;
        sat_seconds:=Cputime(tstart);
        J:=Isat+ideal<R|e>;
        dim_total:=Dimension(Isat);
        dim_special:=Dimension(J);
        G:=Basis(J);

        total:=0; contact_open:=0; smooth:=0; smooth_open:=0;
        tangent_dims:=AssociativeArray();
        by_t:=AssociativeArray();
        samples:=[];
        derivs:=[[Derivative(g,j):j in [1..5]]:g in G];
        for t0 in k do for l0 in k do for u0 in k do for v0 in k do
            pt:=[k!0,t0,l0,u0,v0];
            if not &and[Evaluate(g,pt) eq 0:g in G] then continue; end if;
            total+:=1;
            key:=Z!t0;
            if IsDefined(by_t,key) then by_t[key]+:=1; else by_t[key]:=1; end if;
            openpt:=l0*v0*(u0^2-4*v0^2) ne 0;
            if openpt then contact_open+:=1; end if;
            JM:=Matrix(k,#G,5,
                &cat[[Evaluate(derivs[i][j],pt):j in [1..5]]
                      :i in [1..#G]]);
            td:=5-Rank(JM);
            if IsDefined(tangent_dims,td) then tangent_dims[td]+:=1;
            else tangent_dims[td]:=1; end if;
            smoothpt:=td eq dim_special;
            if smoothpt then smooth+:=1; end if;
            if smoothpt and openpt then
                smooth_open+:=1;
                if #samples lt 12 then
                    Append(~samples,<Z!t0,Z!l0,Z!u0,Z!v0>);
                end if;
            end if;
        end for; end for; end for; end for;

        print "PRIME",p,"CHART",label,description;
        print " RAW_SPECIAL_DIM",raw_dim,
              "SAT_TOTAL_DIM",dim_total,"FLAT_SPECIAL_DIM",dim_special,
              "SAT_SECONDS",sat_seconds;
        print " BASIS_SIZE",#G,
              "BASIS_SHAPES",[<TotalDegree(g),#Terms(g)>:g in G];
        print " F_POINTS",total,"CONTACT_OPEN",contact_open,
              "SMOOTH",smooth,"SMOOTH_CONTACT_OPEN",smooth_open;
        print " TANGENT_DIM_COUNTS",
              Sort([<key,tangent_dims[key]>:key in Keys(tangent_dims)]);
        print " COUNTS_BY_T",Sort([<key,by_t[key]>:key in Keys(by_t)]);
        print " SMOOTH_OPEN_SAMPLES_t_L_U_v",samples;
        if is_infinity then
            print " INFINITY_WARNING this chart keeps L,U,v affine; simultaneous";
            print " poles in cubic-contact variables require additional weighted charts.";
        end if;
    end for;
end for;

print "CONTACT6_M612_BOUNDARY_CORE_MODP_DONE";
quit;
