//////////////////////////////////////////////////////////////////////
// Exact continuation of the four surviving p=5 E9 core points.
//
// Starting from the integral normalized endpoint equations in
// contact6_m612_weighted_E9_lift.m, enumerate correction digits through
// mod 125 and mod 625.  The rank-4 Jacobian is constant modulo 5 over
// each residue class, so its 5^6 possible correction digits are grouped
// once by their five-dimensional image; every compatible right-hand
// side then has exactly 25 preimages.
//
// The weighted-scaling tangent at x=(E,l,u,w,n,r) is
//
//      sigma(x)=(2E,-l,-2u,-2w,-n,-4r).
//
// It is checked to lie in the special-fiber kernel.  First correction
// digits are also grouped modulo this line, exposing the transverse
// correction parameter independently of a choice of kernel basis.
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

function EvalZ(g,pt)
    z:=Evaluate(g,[Q!a:a in pt]);
    assert Denominator(z) eq 1;
    return Z!z;
end function;

function JacobianMod5(pt)
    return Matrix(k,5,6,
        &cat[[k!EvalZ(Derivative(Fs[i],j),pt):j in [1..6]]
             :i in [1..5]]);
end function;

function ImageKey(vals)
    return &+[(Z!vals[i])*p^(i-1):i in [1..#vals]];
end function;

// buckets[key+1] is the complete inverse image of key under J.
function CorrectionBuckets(J)
    buckets:=[[]:i in [0..p^5-1]];
    for dE in [0..p-1] do for dl in [0..p-1] do
    for du in [0..p-1] do for dw in [0..p-1] do
    for dn in [0..p-1] do for dr in [0..p-1] do
        ds:=[dE,dl,du,dw,dn,dr];
        im:=[&+[J[i,j]*k!ds[j]:j in [1..6]]:i in [1..5]];
        key:=ImageKey(im);
        Append(~buckets[key+1],ds);
    end for; end for; end for; end for; end for; end for;
    return buckets;
end function;

function LiftOne(pt,modulus,buckets)
    rhs:=[];
    for f in Fs do
        value:=EvalZ(f,pt);
        assert value mod modulus eq 0;
        Append(~rhs,k!((-value div modulus) mod p));
    end for;
    dslist:=buckets[ImageKey(rhs)+1];
    lifts:=[];
    nextmod:=p*modulus;
    for ds in dslist do
        vals:=[(pt[j]+modulus*ds[j]) mod nextmod:j in [1..6]];
        assert &and[(EvalZ(f,vals) mod nextmod) eq 0:f in Fs];
        Append(~lifts,vals);
    end for;
    return lifts,rhs,dslist;
end function;

function ScalingVector(pt)
    exponents:=[2,-1,-2,-2,-1,-4];
    return [(exponents[j]*pt[j]) mod p:j in [1..6]];
end function;

function AddScaling(ds,sigma,c)
    return [(ds[j]+c*sigma[j]) mod p:j in [1..6]];
end function;

function LexLess(a,b)
    for i in [1..#a] do
        if a[i] lt b[i] then return true; end if;
        if a[i] gt b[i] then return false; end if;
    end for;
    return false;
end function;

function ScalingRepresentative(ds,sigma)
    rep:=ds;
    for c in [1..p-1] do
        q:=AddScaling(ds,sigma,c);
        if LexLess(q,rep) then rep:=q; end if;
    end for;
    return rep;
end function;

print "CONTACT6_M612_WEIGHTED_E9_NEXT_LIFTS";
print "VARIABLE_ORDER",["E","l","u","w","n","r"];
print "MODULI",[5,25,125,625];

grand_counts:=[#survivors];
all_mod25:=[];
all_mod125:=[];
all_mod625:=[];

for base_index in [1..#survivors] do
    base:=survivors[base_index];
    J:=JacobianMod5(base);
    assert Rank(J) eq 4;
    buckets:=CorrectionBuckets(J);
    nonempty:=[i-1:i in [1..#buckets]|#buckets[i] gt 0];
    assert #nonempty eq p^4;
    assert &and[#buckets[i+1] eq p^2:i in nonempty];

    sigma:=ScalingVector(base);
    sigmak:=Vector(k,[k!x:x in sigma]);
    assert sigmak*Transpose(J) eq Vector(k,5,[0,0,0,0,0]);

    lifts25,rhs5,ds25:=LiftOne(base,p,buckets);
    assert #lifts25 eq 25;
    all_mod25 cat:=lifts25;

    // Quotient the 25 first corrections by the scaling line.  Each of
    // the five representatives is one transverse correction class.
    reps:=Setseq({ Sprint(ScalingRepresentative(ds,sigma)):ds in ds25 });
    Sort(~reps);
    assert #reps eq 5;

    lifts125:=[];
    compatible25:=[];
    compatible_by_scaling_rep:=AssociativeArray();
    for pt25 in lifts25 do
        q125,rhs25,ds125:=LiftOne(pt25,p^2,buckets);
        if #q125 gt 0 then
            Append(~compatible25,pt25);
            firstdigit:=[((pt25[j]-base[j]) div p) mod p:j in [1..6]];
            repkey:=Sprint(ScalingRepresentative(firstdigit,sigma));
            if IsDefined(compatible_by_scaling_rep,repkey) then
                compatible_by_scaling_rep[repkey]+:=1;
            else
                compatible_by_scaling_rep[repkey]:=1;
            end if;
            lifts125 cat:=q125;
        end if;
    end for;

    lifts625:=[];
    compatible125:=[];
    for pt125 in lifts125 do
        q625,rhs125,ds625:=LiftOne(pt125,p^3,buckets);
        if #q625 gt 0 then
            Append(~compatible125,pt125);
            lifts625 cat:=q625;
        end if;
    end for;

    all_mod125 cat:=lifts125;
    all_mod625 cat:=lifts625;

    assert #compatible25 eq 5;
    assert #lifts125 eq 125;
    assert #compatible125 eq 25;
    assert #lifts625 eq 625;

    print "BASE_INDEX",base_index,"BASE",base;
    print " JAC_RANK",Rank(J),"SCALING_TANGENT",sigma;
    print " FIRST_CORRECTION_SCALING_QUOTIENT_REPS",reps;
    print " COUNT_MOD25",#lifts25;
    print " MOD25_RESIDUES_COMPATIBLE_WITH_MOD125",#compatible25;
    print " COMPATIBLE_MOD25_BY_SCALING_QUOTIENT",
          Sort([<key,compatible_by_scaling_rep[key]>
                :key in Keys(compatible_by_scaling_rep)]);
    if #compatible25 gt 0 then
        print " COMPATIBLE_MOD25_SAMPLES",compatible25[1..Min(10,#compatible25)];
    end if;
    print " COUNT_MOD125",#lifts125;
    print " MOD125_RESIDUES_COMPATIBLE_WITH_MOD625",#compatible125;
    if #compatible125 gt 0 then
        print " COMPATIBLE_MOD125_SAMPLES",compatible125[1..Min(10,#compatible125)];
    end if;
    print " COUNT_MOD625",#lifts625;
    if #lifts625 gt 0 then
        print " MOD625_SAMPLES",lifts625[1..Min(10,#lifts625)];
    end if;
end for;

grand_counts cat:=[#all_mod25,#all_mod125,#all_mod625];
assert grand_counts eq [4,100,500,2500];
print "TOTAL_COUNTS_MOD_5_25_125_625",grand_counts;
print "GROWTH_FACTORS",
      [grand_counts[i+1] div grand_counts[i]:i in [1..#grand_counts-1]];
print "CONTACT6_M612_WEIGHTED_E9_NEXT_LIFTS_DONE";
quit;
