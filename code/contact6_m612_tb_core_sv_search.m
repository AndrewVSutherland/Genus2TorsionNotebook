//////////////////////////////////////////////////////////////////////
// Height driver for the discovery-oriented (s,v)-slicing of the full
// T_B-halving + cubic-contact [3,3] fiber product.
//
// Each rational (s,v) slice is solved algebraically in (r,K,M,U), so r is
// not height-bounded.  Rational points must have K and M square.  Survivors
// are checked against the p=5 and p=7 boundary wall, reconstructed as an
// order-4 half plus independent order-3 class, simplicity-certified, and
// exact-torsion tested.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned height then height:=4;
elif Type(height) eq MonStgElt then height:=StringToInteger(height); end if;
if not assigned max_slices then max_slices:=0;
elif Type(max_slices) eq MonStgElt then max_slices:=StringToInteger(max_slices); end if;
if not assigned progress_interval then progress_interval:=25;
elif Type(progress_interval) eq MonStgElt then progress_interval:=StringToInteger(progress_interval); end if;
if not assigned max_hits then max_hits:=10;
elif Type(max_hits) eq MonStgElt then max_hits:=StringToInteger(max_hits); end if;
if not assigned use_projective_masks then use_projective_masks:=false;
elif Type(use_projective_masks) eq MonStgElt then
    use_projective_masks:=use_projective_masks in {"true","True","1","yes"};
end if;

Q:=Rationals(); Z:=Integers();
P<x>:=PolynomialRing(Q);
PF<Tfrob>:=PolynomialRing(Q);
load "code/contact6_m612_tb_core_tools.m";

function RationalParametersOfHeight(B)
    vals:=[]; seen:={};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num,den) ne 1 then continue; end if;
            q:=Q!num/den;
            if Sprint(q) notin seen then
                Include(~seen,Sprint(q)); Append(~vals,q);
            end if;
        end for;
    end for;
    return vals;
end function;

// These are the fiberwise-saturated projective masks computed by
// contact6_m612_tb_core_sv_projective_masks.m after K=m^2 and M=L^2.
// A failed residue is only a boundary-deferred slice: specialization need
// not commute with saturation.  Denominator/infinity residues are therefore
// retained automatically, and the final counters never claim that skipped
// slices have been rigorously excluded.
function PassesProjectivePriorityMask(s0,v0)
    masks:=AssociativeArray();
    masks[5]:={<0,0>,<0,2>,<0,3>};
    masks[7]:={<0,i>:i in [0..6]} join {<2,3>,<2,5>,<2,6>};
    masks[11]:={<0,i>:i in [0..10]} join
        {<1,3>,<1,8>,<3,10>,<6,0>,<6,5>,<6,9>,<6,10>,
         <7,2>,<7,5>,<7,8>,<8,10>,<9,4>,<9,5>,<9,7>};
    for p in [5,7,11] do
        if Denominator(s0) mod p eq 0 or Denominator(v0) mod p eq 0 then
            continue;
        end if;
        F:=GF(p);
        ss:=Z!(F!Numerator(s0)/(F!Denominator(s0)));
        vv:=Z!(F!Numerator(v0)/(F!Denominator(v0)));
        if <ss,vv> notin masks[p] then return false,p; end if;
    end for;
    return true,0;
end function;

function Primitive(poly)
    if poly eq 0 then return poly; end if;
    den:=LCM([Denominator(c):c in Coefficients(poly)]);
    g:=Parent(poly)!(den*poly);
    nums:=[Z!c:c in Coefficients(g)];
    cont:=GCD([Abs(n):n in nums|n ne 0]);
    return cont gt 1 select Parent(poly)!(g/cont) else g;
end function;

function ContactPolynomial(a,b)
    h:=1+a*x+b*x^2+x^3;
    return h^2-(x-1)^6,h;
end function;

function FactorDegrees(f)
    return Sort([Degree(ff[1]):ff in Factorization(f)]);
end function;

function IntegralModel(f)
    d:=LCM([Denominator(Coefficient(f,i)):i in [0..Degree(f)]]);
    return P!(d^2*f),d;
end function;

function FrobeniusPolynomial(C,p)
    ef:=EulerFactor(C,p); d:=Degree(ef);
    return &+[Q!Coefficient(ef,i)*Tfrob^(d-i):i in [0..d]];
end function;

function FullSimplicityCertificate(C)
    plist:=[p:p in PrimesUpTo(151)|p ge 3];
    for p in plist do
        try
            Phi:=FrobeniusPolynomial(C,p);
            fac:=Factorization(Phi);
            if #fac ne 1 or fac[1][2] ne 1 or Degree(fac[1][1]) ne 4 then
                continue;
            end if;
            Gal:=GaloisGroup(Phi);
            desc:="unknown";
            try desc:=TransitiveGroupDescription(Gal);
            catch e2 desc:="unknown"; end try;
            if Order(Gal) eq 8 and desc eq "D(4)" then
                return true,"D4",p,Phi;
            end if;
        catch e
            continue;
        end try;
    end for;
    for p in plist do
        try
            Phi:=FrobeniusPolynomial(C,p);
            fac:=Factorization(Phi);
            if #fac ne 1 or fac[1][2] ne 1 or Degree(fac[1][1]) ne 4 then
                continue;
            end if;
            NF<pi>:=NumberField(Phi);
            if &and[Degree(MinimalPolynomial(pi^n)) eq 4:n in [2..12]] then
                return true,"root_power",p,Phi;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false,"none",0,PF!0;
end function;

function OnBoundaryAtPrime(a,b,p)
    if Denominator(a) mod p eq 0 or Denominator(b) mod p eq 0 then
        return true;
    end if;
    kf:=GF(p);
    aa:=kf!Numerator(a)/kf!Denominator(a);
    bb:=kf!Numerator(b)/kf!Denominator(b);
    DB:=(aa-3)^2-8*(bb+3);
    DC:=(bb-3)^2-8*(aa+3);
    RR:=((bb+3)*(aa+3)-4)^2
        -(bb^2-2*aa-3)*(aa^2-2*bb-3);
    return (bb+3)*(aa+3)*DB*DC*RR eq 0;
end function;

function SolveSVSlice(s0,v0)
    A<a,r,K,M,U>:=PolynomialRing(Q,5,"grevlex");
    b0:=2*s0^2-3;
    A3:=a-3+4*s0^2*r-2*K;
    H1:=s0*((a-3)*r^2+4*r-K*(a+3))-r*A3;
    H2:=8*s0^2*(2*s0^2*r^2+2*(a-3)*r+2-K*(2*s0^2-6))
        -A3^2-32*s0^3*r;
    D:=Coefficient(H1,a,1);
    N:=-Evaluate(H1,[A!0,r,K,M,U]);
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

    R<rr,KK,MM,UU>:=PolynomialRing(Q,4,"grevlex");
    toR:=hom<A -> R | R!0,rr,KK,MM,UU>;
    G0:=Primitive(R!toR(HTB));
    G1:=Primitive(R!toR(SubA(F1)));
    G2:=Primitive(R!toR(SubA(F2)));
    G3:=Primitive(R!toR(SubA(F3)));
    D4:=R!toR(D); N4:=R!toR(N);
    I:=ideal<R|G0,G1,G2,G3>;
    boundary:=rr*KK*MM*D4*(UU^2-4*v0^2)*(N4+(b0+2)*D4);
    sat_ok:=true;
    try I:=Saturation(I,ideal<R|boundary>);
    catch e sat_ok:=false; end try;
    dim:=Dimension(I);
    if dim eq -1 then return [],dim,sat_ok,true; end if;
    if dim ne 0 then return [],dim,sat_ok,false; end if;
    pts:=[]; variety_ok:=false;
    try pts:=Variety(I); variety_ok:=true;
    catch e variety_ok:=false; end try;
    out:=[];
    if variety_ok then
        for pt in pts do
            r0:=Q!pt[1]; K0:=Q!pt[2]; M0:=Q!pt[3]; U0:=Q!pt[4];
            d0:=Q!Evaluate(D4,<r0,K0,M0,U0>);
            if d0 eq 0 then continue; end if;
            a0:=Q!Evaluate(N4,<r0,K0,M0,U0>)/d0;
            Append(~out,<r0,K0,M0,U0,a0,b0>);
        end for;
    end if;
    return out,dim,sat_ok,variety_ok;
end function;

function VerifyCandidate(s0,v0,rec)
    r0,K0,M0,U0,a0,b0:=Explode(rec);
    sqK,m0:=IsSquare(K0); sqM,L0:=IsSquare(M0);
    if not sqK or not sqM or K0 eq 0 or M0 eq 0 then
        return false,[],false,0,P!0;
    end if;
    f,h6:=ContactPolynomial(a0,b0);
    if Degree(f) ne 5 or Discriminant(f) eq 0 or
       FactorDegrees(f) ne [1,2,2] then
        return false,[],false,0,P!0;
    end if;
    if not OnBoundaryAtPrime(a0,b0,5) or not OnBoundaryAtPrime(a0,b0,7) then
        print "LOCAL_WALL_FAILURE","s",s0,"v",v0,"a",a0,"b",b0;
        return false,[],false,0,P!0;
    end if;
    okcore,L,ordD,ordE,E:=M612_VerifyCorePoint(a0,b0,M0,U0,v0);
    if not okcore then return false,[],false,0,P!0; end if;

    J:=Jacobian(HyperellipticCurve(f));
    D6:=J![x-1,Evaluate(h6,1)];
    qB:=x^2+((a0-3)/(b0+3))*x+2/(b0+3);
    A30:=a0-3+4*s0^2*r0-2*K0;
    p0:=A30/(4*s0^2);
    Ghalf:=x^2+p0*x+r0/s0;
    Mtool:=2*s0^2/m0;
    Ntool:=2*s0^2*r0/m0;
    ell:=qB*(Mtool*x+Ntool);
    TB:=J![qB,P!0];
    H4:=J![Ghalf,-ell mod Ghalf];
    if 3*D6 ne TB or 2*H4 ne TB or Order(H4) ne 4 then
        return false,[],false,0,P!0;
    end if;

    fI,d:=IntegralModel(f);
    CI:=HyperellipticCurve(fI);
    simple,method,pcert,Phi:=FullSimplicityCertificate(CI);
    if not simple then return true,[],false,0,fI; end if;
    TG,mp:=TorsionSubgroup(Jacobian(CI));
    return true,Invariants(TG),true,pcert,fI;
end function;

params:=RationalParametersOfHeight(height);
checked:=0; enumerated:=0; maskDeferred:=0; emptyDim:=0; dim0:=0; positiveDim:=0;
saturationFail:=0; varietyFail:=0;
rationalPoints:=0; KSquares:=0; MSquares:=0; bothSquares:=0;
verified:=0; simpleCount:=0; exactTests:=0; hits:=[]; seen:={};
start:=Cputime();

print "CONTACT6_M612_TB_CORE_SV_SEARCH";
print "height",height,"parameter_count",#params,"max_slices",max_slices,
      "use_projective_masks",use_projective_masks;
for s0 in params do
    if s0 eq 0 then continue; end if;
    for v0 in params do
        if v0 eq 0 then continue; end if;
        enumerated+:=1;
        if use_projective_masks then
            priority,badp:=PassesProjectivePriorityMask(s0,v0);
            if not priority then
                maskDeferred+:=1;
                continue;
            end if;
        end if;
        if max_slices gt 0 and checked ge max_slices then break s0; end if;
        if #hits ge max_hits then break s0; end if;
        checked+:=1;
        t0:=Cputime();
        recs,dim,sat_ok,variety_ok:=SolveSVSlice(s0,v0);
        elapsed:=Cputime(t0);
        if not sat_ok then saturationFail+:=1; end if;
        if dim eq -1 then emptyDim+:=1;
        elif dim eq 0 then dim0+:=1;
        else positiveDim+:=1; end if;
        if dim eq 0 and not variety_ok then varietyFail+:=1; end if;
        if dim gt 0 or not sat_ok or (dim eq 0 and not variety_ok) then
            print "EXCEPTIONAL_SLICE","s",s0,"v",v0,"dim",dim,
                  "sat",sat_ok,"variety",variety_ok,"seconds",elapsed;
        end if;
        rationalPoints+:=#recs;
        for rec in recs do
            r0,K0,M0,U0,a0,b0:=Explode(rec);
            sqK,m0:=IsSquare(K0); sqM,L0:=IsSquare(M0);
            if sqK then KSquares+:=1; end if;
            if sqM then MSquares+:=1; end if;
            if not sqK or not sqM then
                print "RATIONAL_NONSQUARE","s",s0,"v",v0,"r",r0,
                      "K",K0,"M",M0,"U",U0,"a",a0,"b",b0,
                      "Ksquare",sqK,"Msquare",sqM;
                continue;
            end if;
            bothSquares+:=1;
            key:=Sprint(<a0,b0>);
            if key in seen then continue; end if;
            Include(~seen,key);
            ok,inv,simple,pcert,fI:=VerifyCandidate(s0,v0,rec);
            if not ok then continue; end if;
            verified+:=1;
            if simple then simpleCount+:=1; exactTests+:=1; end if;
            print "VERIFIED_INTERSECTION","s",s0,"v",v0,"r",r0,
                  "K",K0,"M",M0,"U",U0,"a",a0,"b",b0,
                  "simple",simple,"pcert",pcert,"torsion",inv;
            if simple and inv eq [6,12] then
                Append(~hits,<s0,v0,r0,K0,M0,U0,a0,b0,pcert,fI>);
                print "HIT_6_12",hits[#hits];
            end if;
        end for;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress",checked,"empty",emptyDim,"dim0",dim0,
                  "positive",positiveDim,
                  "Qpoints",rationalPoints,"both_squares",bothSquares,
                  "verified",verified,"hits",#hits,"last_seconds",elapsed;
        end if;
    end for;
end for;

print "DONE","enumerated",enumerated,"checked",checked,
      "mask_boundary_deferred",maskDeferred,
      "empty_dim",emptyDim,"dim0",dim0,
      "positive_dim",positiveDim,
      "saturation_fail",saturationFail,"variety_fail",varietyFail,
      "rational_points",rationalPoints,"K_squares",KSquares,
      "M_squares",MSquares,"both_squares",bothSquares,
      "verified",verified,"simple",simpleCount,"exact_tests",exactTests,
      "hits",#hits,"seconds",Cputime(start);
for H in hits do print "H",H; end for;
quit;
