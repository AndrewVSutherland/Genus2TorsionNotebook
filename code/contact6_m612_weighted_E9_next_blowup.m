//////////////////////////////////////////////////////////////////////
// A regular first blowup chart for each of the four surviving p=5 E9
// points of the contact-6 [6,12] core.
//
// Put x=x0+5*d and G_i(d)=F_i(x)/5.  At a surviving x0, the five
// affine-linear reductions G_i mod 5 have rank four and one left
// relation lambda.  The polynomial sum lambda_i*G_i is identically
// divisible by 5, so adjoin the divided relation
//
//              H=(sum lambda_i G_i)/5.
//
// Drop one G_i whose lambda coefficient is a unit.  The remaining four
// G_i together with H are exactly equivalent over Z_5 to all five G_i.
// Their mod-5 zeros are the corrections that continue two more raw
// congruence layers.  A rank-five Jacobian at such a zero proves a
// smooth one-dimensional Z_5 branch by multivariate Hensel.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Z:=Integers();
Q:=Rationals();
p:=5;
k:=GF(p);

A<E,l,u,w,n,r>:=PolynomialRing(Q,6,"grevlex");
F5:=2*p*n-3*u-6*l^2;
F4:=E*(p^2*n^2+2*r-3*u^2-3*w^2+15*p^2*l^2)-2*l^2;
F3:=2*w^3+2*p*n*r-u^3-6*u*w^2-22*p^4*l^2;
F2:=E^2*(r^2+2*p*n*w^3-3*u^2*w^2-3*w^4
          +15*p^6*l^2)-p^2*l^2;
F1:=E*(2*r*w^3-3*u*w^4-6*p^8*l^2)-2*p^6*l^2;
Fs:=[F5,F4,F3,F2,F1];

survivors:=[
    [1,2,2,2,3,1], [1,3,2,2,3,1],
    [4,1,3,3,2,1], [4,4,3,3,2,1]
];

D<dE,dl,du,dw,dn,dr>:=PolynomialRing(Q,6,"grevlex");
dsvars:=[dE,dl,du,dw,dn,dr];
Db<cE,cl,cu,cw,cn,cr>:=PolynomialRing(k,6,"grevlex");
reduceD:=hom<D -> Db | [cE,cl,cu,cw,cn,cr]>;

function EvalZ(g,pt)
    z:=Evaluate(g,[Q!a:a in pt]);
    assert Denominator(z) eq 1;
    return Z!z;
end function;

function EvalDZ(g,pt)
    z:=Evaluate(g,[Q!a:a in pt]);
    assert Denominator(z) eq 1;
    return Z!z;
end function;

function JacobianMod5AtA(pt)
    return Matrix(k,5,6,
        &cat[[k!EvalZ(Derivative(Fs[i],j),pt):j in [1..6]]
             :i in [1..5]]);
end function;

function LeftRelation(J)
    zero:=Vector(k,6,[0,0,0,0,0,0]);
    for a1 in [0..p-1] do for a2 in [0..p-1] do
    for a3 in [0..p-1] do for a4 in [0..p-1] do
    for a5 in [0..p-1] do
        aa:=[a1,a2,a3,a4,a5];
        if aa eq [0,0,0,0,0] then continue; end if;
        v:=Vector(k,[k!x:x in aa]);
        if v*J ne zero then continue; end if;
        first:=Min([i:i in [1..5]|aa[i] ne 0]);
        inv:=Z!(1/(k!aa[first]));
        return [(inv*aa[i]) mod p:i in [1..5]];
    end for; end for; end for; end for; end for;
    error "left relation not found";
end function;

function DropRowAndCheck(J,lambda)
    for drop in [1..5] do
        if lambda[drop] eq 0 then continue; end if;
        keep:=[i:i in [1..5]|i ne drop];
        M:=Matrix(k,4,6,&cat[[J[i,j]:j in [1..6]]:i in keep]);
        if Rank(M) eq 4 then return drop,keep; end if;
    end for;
    error "no unit-coefficient dependent row can be dropped";
end function;

function CoefficientsIntegral(g)
    return &and[Denominator(c) eq 1:c in Coefficients(g)];
end function;

print "CONTACT6_M612_WEIGHTED_E9_NEXT_BLOWUP";
print "F_ORDER",["F5","F4","F3","F2","F1"];
print "CORRECTION_ORDER",["dE","dl","du","dw","dn","dr"];

total_smooth:=0;
for base_index in [1..#survivors] do
    base:=survivors[base_index];
    J:=JacobianMod5AtA(base);
    assert Rank(J) eq 4;
    lambda:=LeftRelation(J);
    drop,keep:=DropRowAndCheck(J,lambda);

    phi:=hom<A -> D |
        [D!(base[j]+p*dsvars[j]):j in [1..6]]>;
    Gs:=[];
    for f in Fs do
        numerator:=phi(f);
        G:=numerator/p;
        assert CoefficientsIntegral(G);
        Append(~Gs,G);
    end for;

    Hnumerator:=&+[Z!lambda[i]*Gs[i]:i in [1..5]];
    H:=Hnumerator/p;
    assert CoefficientsIntegral(H);
    proper:=[Gs[i]:i in keep] cat [H];

    smooth_points:=[];
    tangent_vectors:=[];
    for a1 in [0..p-1] do for a2 in [0..p-1] do
    for a3 in [0..p-1] do for a4 in [0..p-1] do
    for a5 in [0..p-1] do for a6 in [0..p-1] do
        pt:=[a1,a2,a3,a4,a5,a6];
        if not &and[(EvalDZ(g,pt) mod p) eq 0:g in proper] then
            continue;
        end if;
        JP:=Matrix(k,5,6,
            &cat[[k!EvalDZ(Derivative(proper[i],j),pt):j in [1..6]]
                 :i in [1..5]]);
        assert Rank(JP) eq 5;
        // Record the unique projective tangent direction by brute force,
        // normalized at its first nonzero coordinate.
        tangent:=[];
        for t1 in [0..p-1] do for t2 in [0..p-1] do
        for t3 in [0..p-1] do for t4 in [0..p-1] do
        for t5 in [0..p-1] do for t6 in [0..p-1] do
            tt:=[t1,t2,t3,t4,t5,t6];
            if tt eq [0,0,0,0,0,0] then continue; end if;
            tv:=Vector(k,[k!x:x in tt]);
            if tv*Transpose(JP) ne Vector(k,5,[0,0,0,0,0]) then
                continue;
            end if;
            first:=Min([i:i in [1..6]|tt[i] ne 0]);
            inv:=Z!(1/(k!tt[first]));
            normalized:=[(inv*tt[i]) mod p:i in [1..6]];
            if #tangent eq 0 then
                tangent:=normalized;
            else
                assert tangent eq normalized;
            end if;
        end for; end for; end for; end for; end for; end for;
        assert #tangent eq 6;
        Append(~smooth_points,pt);
        Append(~tangent_vectors,tangent);
    end for; end for; end for; end for; end for; end for;

    assert #smooth_points eq 5;
    total_smooth+:=#smooth_points;
    print "BASE_INDEX",base_index,"BASE",base;
    print " LEFT_RELATION_lambda",lambda;
    print " DROP_ROW",drop,"KEEP_ROWS",keep;
    print " DIVIDED_RELATION_H_SHAPE",<TotalDegree(H),#Terms(H)>;
    print " DIVIDED_RELATION_H_MOD5",reduceD(H);
    print " SMOOTH_POINT_COUNT",#smooth_points;
    print " SMOOTH_CORRECTIONS_AND_TANGENTS",
          [<smooth_points[i],tangent_vectors[i]>:i in [1..#smooth_points]];
    print " PROPER_JACOBIAN_RANKS",[5:i in [1..#smooth_points]];
end for;

assert total_smooth eq 20;
print "TOTAL_SMOOTH_MOD5_POINTS",total_smooth;
print "VERDICT",
      "each of the four E9 survivors has five smooth one-dimensional Z_5 branches";
print "CONTACT6_M612_WEIGHTED_E9_NEXT_BLOWUP_DONE";
quit;
