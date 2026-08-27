//////////////////////////////////////////////////////////////////////
// Weighted contact-variable poles above base-projective infinity for
// the contact-6 [6,12] core.
//
// Keep the two auxiliary coefficients N,R instead of eliminating them:
//
//   q = x^2 + U*x + v^2,
//   H = x^3 + N*x^2 + R*x + v^3,
//   H^2 - q^3 = L^2*f.
//
// Coefficient comparison gives five much smaller equations.  Their
// weighted initial forms detect simultaneous poles in (L,U,v) without
// the excess components introduced by the eliminated F1,F2,F3 model.
//
// The generic B != 0 infinity chart has only the pure H^2=q^3 initial
// form and one x^5-forced form.  The latter is checked over F_5,F_7.
// At the exceptional endpoint [A:B:T]=[1:0:0], ten exact initial-form
// signatures found by the valuation fan scan are checked.  The leading
// unit of e is retained as E; setting E=1 loses rational residue orbits.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Z:=Integers();

if assigned primes then
    if Type(primes) eq MonStgElt then
        prime_list:=[StringToInteger(s):s in Split(primes,",")|#s gt 0];
    else
        prime_list:=primes;
    end if;
else
    prime_list:=[5,7];
end if;

// Every weight is (v(e),v(L),v(U),v(v),v(N),v(R)).
generic_weights:=[
    <1,-1,-4,-4,-4,-8>,
    <1,-1,-3,-3,-3,-6>,
    <1,-1,-3,-3,-2,-6>
];

endpoint_weights:=[
    <1,-2,-5,-5,-5,-10>,
    <4,-1,-3,-3,-3,-6>,
    <1,-2,-4,-4,-4,-8>,
    <2,-1,-2,-2,-2,-4>,
    <6,-3,-2,-4,-6,-6>,
    <8,-4,-7,-6,-8,-13>,
    <6,-3,-5,-4,-6,-9>,
    <6,-3,-5,-1,-6,-9>,
    <1,-2,-4,-4,-3,-8>,
    <2,-1,-2,-2,-1,-4>
];

function MinimumIndices(vals)
    m:=Min(vals);
    return [i:i in [1..#vals]|vals[i] eq m];
end function;

function GenericSignature(w)
    ee,l,u,v,n,r:=Explode([w[i]:i in [1..6]]);
    vals:=[
        [n,u,-ee+2*l],
        [2*n,r,2*u,2*v,-2*ee+2*l],
        [3*v,n+r,3*u,u+2*v,-2*ee+2*l],
        [2*r,n+3*v,2*u+2*v,4*v,-2*ee+2*l],
        [r+3*v,u+4*v,-ee+2*l]
    ];
    return [MinimumIndices(z):z in vals];
end function;

function EndpointSignature(w)
    ee,l,u,v,n,r:=Explode([w[i]:i in [1..6]]);
    // At t=0 in [A:B:T]=[1:t:e], the coefficient valuations of
    // (c5,c4,c3,c2,c1) are (0,-ee,0,-2ee,-ee).
    vals:=[
        [n,u,2*l],
        [2*n,r,2*u,2*v,-ee+2*l],
        [3*v,n+r,3*u,u+2*v,2*l],
        [2*r,n+3*v,2*u+2*v,4*v,-2*ee+2*l],
        [r+3*v,u+4*v,-ee+2*l]
    ];
    return [MinimumIndices(z):z in vals];
end function;

function SumSelected(row,inds)
    return &+[row[i]:i in inds];
end function;

// Equations are multiplied by E^2.  This does not alter their common
// zero set on the E != 0 torus and avoids Laurent-polynomial syntax.
function GenericInitialEquations(S,sig)
    E,L,U,V,N,R:=Explode([S.i:i in [1..6]]);
    rows:=[
        [2*N,-3*U,-2*L^2/E],
        [N^2,2*R,-3*U^2,-3*V^2,-L^2/E^2],
        [2*V^3,2*N*R,-U^3,-6*U*V^2,-2*L^2/E^2],
        [R^2,2*N*V^3,-3*U^2*V^2,-3*V^4,-L^2/E^2],
        [2*R*V^3,-3*U*V^4,-2*L^2/E]
    ];
    K:=FieldOfFractions(S);
    ans:=[];
    for i in [1..5] do
        g:=(K!E)^2*(&+[K!rows[i][j]:j in sig[i]]);
        assert Denominator(g) eq 1;
        Append(~ans,S!Numerator(g));
    end for;
    return ans;
end function;

function EndpointInitialEquations(S,sig)
    E,L,U,V,N,R:=Explode([S.i:i in [1..6]]);
    K:=FieldOfFractions(S);
    rows:=[
        [2*N,-3*U,-6*L^2],
        [N^2,2*R,-3*U^2,-3*V^2,-2*L^2/E],
        [2*V^3,2*N*R,-U^3,-6*U*V^2,-22*L^2],
        [R^2,2*N*V^3,-3*U^2*V^2,-3*V^4,-L^2/E^2],
        [2*R*V^3,-3*U*V^4,-2*L^2/E]
    ];
    ans:=[];
    for i in [1..5] do
        g:=(K!E)^2*(&+[K!rows[i][j]:j in sig[i]]);
        assert Denominator(g) eq 1;
        Append(~ans,S!Numerator(g));
    end for;
    return ans;
end function;

function OpenTorusData(S,eqs,p)
    k:=BaseRing(S);
    E,L,U,V,N,R:=Explode([S.i:i in [1..6]]);
    open_product:=E*L*U*V*N*R*(U^2-4*V^2);
    I:=ideal<S|eqs>;
    Iopen:=Saturation(I,ideal<S|open_product>);
    empty:=S!1 in Iopen;
    dim:=empty select -1 else Dimension(Iopen);

    total:=0;
    ranks:=AssociativeArray();
    samples:=[];
    square_data:=AssociativeArray();
    derivs:=[[Derivative(eqs[i],j):j in [1..6]]:i in [1..5]];
    nz:=[x:x in k|x ne 0];
    for e0 in nz do for l0 in nz do for u0 in nz do
    for v0 in nz do for n0 in nz do for r0 in nz do
        pt:=[e0,l0,u0,v0,n0,r0];
        if u0^2-4*v0^2 eq 0 then continue; end if;
        if not &and[Evaluate(g,pt) eq 0:g in eqs] then continue; end if;
        total+:=1;
        J:=Matrix(k,5,6,
            &cat[[Evaluate(derivs[i][j],pt):j in [1..6]]
                  :i in [1..5]]);
        rk:=Rank(J);
        if IsDefined(ranks,rk) then ranks[rk]+:=1; else ranks[rk]:=1; end if;
        // At [1:0:0], DB is always a local square.  When v(e) is even,
        // the leading unit of DC is -8/E.  Record its residue squareclass.
        dc_square:=IsSquare(k!(-8)/e0);
        key:=<Z!e0,dc_square>;
        if IsDefined(square_data,key) then square_data[key]+:=1;
        else square_data[key]:=1;
        end if;
        if #samples lt 24 then
            Append(~samples,<Z!e0,Z!l0,Z!u0,Z!v0,Z!n0,Z!r0>);
        end if;
    end for; end for; end for; end for; end for; end for;
    return total,ranks,samples,empty,dim,square_data;
end function;

print "CONTACT6_M612_WEIGHTED_INFINITY_LEADING";
print "VARIABLE_ORDER E,L,U,V,N,R";
print "IDENTITY H^2-q^3=L^2*f";

for p0 in prime_list do
    p:=Z!p0;
    require p in {5,7}: "this diagnostic is intended for p=5,7";
    k:=GF(p);
    S<E,L,U,V,N,R>:=PolynomialRing(k,6,"grevlex");
    print "PRIME",p;

    print " GENERIC_B_NONZERO";
    generic_counts:=[];
    for j in [1..#generic_weights] do
        w:=generic_weights[j];
        sig:=GenericSignature(w);
        eqs:=GenericInitialEquations(S,sig);
        total,ranks,samples,empty,dim,sq:=OpenTorusData(S,eqs,p);
        print "  G_SIGNATURE",j-1,"WEIGHT",w,"MIN_INDICES",sig;
        print "   OPEN_EMPTY",empty,"OPEN_DIM",dim,"OPEN_POINTS",total,
              "RANK_COUNTS",Sort([<q,ranks[q]>:q in Keys(ranks)]);
        if total gt 0 then print "   SAMPLES_E_L_U_V_N_R",samples; end if;
        Append(~generic_counts,total);
    end for;
    assert generic_counts eq [0,0,0];

    print " ENDPOINT_[1:0:0]";
    endpoint_counts:=[];
    for j in [1..#endpoint_weights] do
        w:=endpoint_weights[j];
        sig:=EndpointSignature(w);
        eqs:=EndpointInitialEquations(S,sig);
        total,ranks,samples,empty,dim,sq:=OpenTorusData(S,eqs,p);
        print "  E_SIGNATURE",j-1,"WEIGHT",w,"MIN_INDICES",sig;
        print "   OPEN_EMPTY",empty,"OPEN_DIM",dim,"OPEN_POINTS",total,
              "RANK_COUNTS",Sort([<q,ranks[q]>:q in Keys(ranks)]);
        if total gt 0 then
            assert endpoint_weights[j][1] mod 2 eq 0;
            assert &and[not key[2]:key in Keys(sq)];
            print "   SAMPLES_E_L_U_V_N_R",samples;
            print "   DC_LEADING_SQUARE_BY_E",
                  Sort([<key,sq[key]>:key in Keys(sq)]);
            print "   NOTE DB is square; R3 passes its square cover.";
            print "   NOTE R2 and R1 pass only when DC leading square is true";
            print "        and v(e) is even.";
        end if;
        Append(~endpoint_counts,total);
    end for;
    if p eq 5 then
        assert endpoint_counts eq [0,0,0,0,0,0,0,0,0,16];
    else
        assert endpoint_counts eq [0,0,0,6,4,0,8,24,0,0];
    end if;
end for;

print "CONTACT6_M612_WEIGHTED_INFINITY_LEADING_DONE";
quit;
