//////////////////////////////////////////////////////////////////////
// Search the cubic-contact [3,3] core for the cheap discriminant covers
// forced by halving the two simple distinguished Richelot-kernel classes.
//
// If R2=2*x^2-(a+3) is divisible by 2 on the distinguished dual, then
//
//   DC=(b-3)^2-8*(a+3)
//
// is a rational square.  Symmetrically, halving R3 forces
//
//   DB=(a-3)^2-8*(b+3)
//
// to be a rational square.  These are exactly the discriminants of the
// original contact factors C and B.  The mixed factor R1 can halve only if
// DB*DC is a square.  Earlier [6,6] searches rejected many of the first two
// covers because their source factor type is not [1,2,2]; here they are the
// intended input to the Richelot quotient.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned height then height:=10;
elif Type(height) eq MonStgElt then height:=StringToInteger(height); end if;
if not assigned progress_interval then progress_interval:=250;
elif Type(progress_interval) eq MonStgElt then
    progress_interval:=StringToInteger(progress_interval);
end if;
if not assigned max_hits then max_hits:=10;
elif Type(max_hits) eq MonStgElt then max_hits:=StringToInteger(max_hits); end if;

Q:=Rationals(); Z:=Integers();
P<x>:=PolynomialRing(Q); PF<Tfrob>:=PolynomialRing(Q);
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

function ContactPolynomial(a,b)
    h:=1+a*x+b*x^2+x^3;
    return h^2-(x-1)^6,h;
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
            Gal:=GaloisGroup(Phi); desc:="unknown";
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

function GenericSliceSolutions(b0,v0)
    R<M,U>:=PolynomialRing(Q,2,"grevlex");
    K:=FieldOfFractions(R); PA<a>:=PolynomialRing(K);
    b:=K!b0; v:=K!v0;
    c1:=2*a+6;
    c2:=a^2+2*b-15;
    c3:=2*a*b+22;
    c4:=2*a+b^2-15;
    c5:=2*b+6;
    B3:=c5*M+3*U;
    Delta3:=4*c4*M+12*(U^2+v^2)-B3^2;
    F3:=B3*Delta3+16*v^3-8*c3*M-8*U^3-48*U*v^2;
    F2:=Delta3^2+64*B3*v^3-64*c2*M-192*(U^2*v^2+v^4);
    F1:=Delta3*v^3-4*c1*M-12*U*v^4;
    D:=Coefficient(F1,1); N:=-Coefficient(F1,0);
    if D eq 0 then return [],false,-2; end if;

    function SubNum(poly)
        d:=Degree(poly); value:=K!0;
        for i in [0..d] do
            value+:=Coefficient(poly,i)*N^i*D^(d-i);
        end for;
        return M612_PrimitivePolynomial(R!Numerator(value));
    end function;

    G2:=SubNum(F2); G3:=SubNum(F3);
    boundary:=M*(U^2-4*v0^2)*(R!(N+(b0+2)*D));
    I:=ideal<R|G2,G3>; sat_ok:=true;
    try I:=Saturation(I,ideal<R|boundary>);
    catch e sat_ok:=false; end try;
    dim:=Dimension(I);
    if not sat_ok or dim ne 0 then return [],sat_ok,dim; end if;
    pts:=[];
    try pts:=Variety(I);
    catch e return [],sat_ok,-3; end try;
    out:=[];
    for pt in pts do
        M0:=Q!pt[1]; U0:=Q!pt[2];
        if M0 eq 0 then continue; end if;
        sq,L0:=IsSquare(M0);
        if not sq or L0 eq 0 then continue; end if;
        Nnum:=R!Numerator(N); Nden:=R!Denominator(N);
        Dnum:=R!Numerator(D); Dden:=R!Denominator(D);
        denN:=Evaluate(Nden,<M0,U0>); denD:=Evaluate(Dden,<M0,U0>);
        if denN eq 0 or denD eq 0 or Evaluate(Dnum,<M0,U0>) eq 0 then
            continue;
        end if;
        a0:=(Evaluate(Nnum,<M0,U0>)/denN) /
            (Evaluate(Dnum,<M0,U0>)/denD);
        Append(~out,<a0,b0,M0,U0,v0>);
    end for;
    return out,sat_ok,dim;
end function;

// On the rational special branch v=1, F1 loses its a-term.  Use F1 itself
// together with F2 after recovering a from the (generically linear) F3.
function SpecialSliceSolutions(b0)
    R<L,U>:=PolynomialRing(Q,2,"grevlex");
    K:=FieldOfFractions(R); PA<a>:=PolynomialRing(K);
    b:=K!b0; v:=K!1; M:=L^2;
    c1:=2*a+6;
    c2:=a^2+2*b-15;
    c3:=2*a*b+22;
    c4:=2*a+b^2-15;
    c5:=2*b+6;
    B3:=c5*M+3*U;
    Delta3:=4*c4*M+12*(U^2+v^2)-B3^2;
    F3:=B3*Delta3+16*v^3-8*c3*M-8*U^3-48*U*v^2;
    F2:=Delta3^2+64*B3*v^3-64*c2*M-192*(U^2*v^2+v^4);
    F1:=Delta3*v^3-4*c1*M-12*U*v^4;
    if Degree(F1) ne 0 then return [],false,-2; end if;
    H:=M612_PrimitivePolynomial(R!Coefficient(F1,0));
    D:=Coefficient(F3,1); N:=-Coefficient(F3,0);
    if D eq 0 then return [],false,-2; end if;
    d:=Degree(F2); value:=K!0;
    for i in [0..d] do value+:=Coefficient(F2,i)*N^i*D^(d-i); end for;
    J:=M612_PrimitivePolynomial(R!Numerator(value));
    DR:=R!Numerator(D); NR:=R!Numerator(N);
    I:=ideal<R|H,J>; sat_ok:=true;
    boundary:=L*(U^2-4)*(b0+3)*DR*(NR+(b0+2)*DR);
    if boundary eq 0 then return [],false,-4; end if;
    try I:=Saturation(I,ideal<R|boundary>);
    catch e sat_ok:=false; end try;
    dim:=Dimension(I);
    if not sat_ok or dim ne 0 then return [],sat_ok,dim; end if;
    pts:=[];
    try pts:=Variety(I);
    catch e return [],sat_ok,-3; end try;
    out:=[];
    for pt in pts do
        L0:=Q!pt[1]; U0:=Q!pt[2];
        if L0 eq 0 then continue; end if;
        d0:=Evaluate(DR,<L0,U0>);
        if d0 eq 0 then continue; end if;
        a0:=Evaluate(NR,<L0,U0>)/d0;
        Append(~out,<a0,b0,L0^2,U0,Q!1>);
    end for;
    return out,sat_ok,dim;
end function;

function TorsionData(C)
    f,h:=HyperellipticPolynomials(C);
    if h ne 0 then return false,[],P!0; end if;
    fI,d:=IntegralModel(P!f);
    try
        G,mp:=TorsionSubgroup(Jacobian(HyperellipticCurve(fI)));
        return true,Invariants(G),fI;
    catch e
        return false,[],fI;
    end try;
end function;

params:=RationalParametersOfHeight(height);
checked:=0; exceptional:=0; formalLifts:=0; formalSplitCover:=0;
verifiedCore:=0; splitCover:=0; simpleSources:=0;
richelotJacobians:=0; exactTests:=0;
hits:=[]; coverSources:={}; seenVerified:={}; seenTargets:={}; start:=Cputime();

print "CONTACT6_M612_DUAL_SPLIT_CORE_SEARCH";
print "height",height,"parameter_count",#params;
for b0 in params do
    // The x^5 coefficient of h^2-(x-1)^6 is 2*(b+3).  Thus b=-3 is a
    // degree-drop boundary, not a genus-2 source.  In particular it must
    // be removed explicitly: multiplying a saturation boundary by the
    // specialized scalar b0+3 does not remove this fiber.
    if b0 eq -3 then continue; end if;
    for v0 in params do
        if v0 eq 0 then continue; end if;
        if #hits ge max_hits then break b0; end if;
        checked+:=1;
        if v0 eq 1 then
            recs,sat_ok,dim:=SpecialSliceSolutions(b0);
        else
            recs,sat_ok,dim:=GenericSliceSolutions(b0,v0);
        end if;
        // Dimension -1 is the empty ideal.  The values -2 and -3 are the
        // solver's recovery/Variety failure sentinels.
        if not sat_ok or dim gt 0 or dim lt -1 then
            exceptional+:=1;
            print "EXCEPTIONAL_SLICE","b",b0,"v",v0,"sat",sat_ok,"dim",dim;
        end if;
        formalLifts+:=#recs;
        for rec in recs do
            a,b,M,U,v:=Explode(rec);
            DB:=(a-3)^2-8*(b+3);
            DC:=(b-3)^2-8*(a+3);
            sqB,wB:=IsSquare(DB); sqC,wC:=IsSquare(DC);
            sqMix,wMix:=IsSquare(DB*DC);
            if DB eq 0 then sqB:=false; end if;
            if DC eq 0 then sqC:=false; end if;
            if DB*DC eq 0 then sqMix:=false; end if;
            if sqB or sqC or sqMix then formalSplitCover+:=1; end if;

            // Verify smoothness, the cubic-contact identity, order, and
            // independence before interpreting a square-cover count as a
            // [3,3] source.  Earlier output applied this only after the
            // cover filter, so degree-drop and singular boundary points
            // dominated core_lifts/split_cover while verified_core stayed
            // zero.
            okcore,L,ordD,ordE,E:=M612_VerifyCorePoint(a,b,M,U,v);
            if not okcore then continue; end if;
            verifiedCore+:=1;
            if not sqB and not sqC and not sqMix then continue; end if;
            splitCover+:=1;
            sourceKey:=Sprint(<a,b>);
            Include(~coverSources,sourceKey);
            if sourceKey in seenVerified then continue; end if;
            Include(~seenVerified,sourceKey);
            f,h6:=ContactPolynomial(a,b); fI,d:=IntegralModel(f);
            CI:=HyperellipticCurve(fI);
            simple,method,pcert,Phi:=FullSimplicityCertificate(CI);
            print "DUAL_SPLIT_SOURCE","a",a,"b",b,"M",M,"U",U,"v",v,
                  "DBsquare",sqB,"DB",DB,"DCsquare",sqC,"DC",DC,
                  "DBDCsquare",sqMix,
                  "simple",simple,"method",method,"pcert",pcert;
            if not simple then continue; end if;
            simpleSources+:=1;
            J:=Jacobian(HyperellipticCurve(f));
            Rs:=RichelotIsogenousSurfaces(J);
            for i in [1..#Rs] do
                if Type(Rs[i]) ne JacHyp then continue; end if;
                richelotJacobians+:=1;
                good,inv,gI:=TorsionData(Curve(Rs[i]));
                if not good then continue; end if;
                exactTests+:=1;
                targetKey:=Sprint(gI);
                if targetKey in seenTargets then continue; end if;
                Include(~seenTargets,targetKey);
                print "RICHELOT_TARGET","source_a",a,"source_b",b,
                      "index",i,"torsion",inv,"curve",gI;
                if inv eq [6,12] then
                    Append(~hits,<a,b,M,U,v,sqB,DB,sqC,DC,sqMix,pcert,i,gI>);
                    print "HIT_6_12",hits[#hits];
                end if;
            end for;
        end for;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress",checked,"formal_lifts",formalLifts,
                  "formal_split_cover",formalSplitCover,
                  "verified_core",verifiedCore,"split_cover",splitCover,
                  "simple_sources",simpleSources,"exact_tests",exactTests,
                  "hits",#hits;
        end if;
    end for;
end for;

print "DONE","checked",checked,"exceptional",exceptional,
      "formal_lifts",formalLifts,"formal_split_cover",formalSplitCover,
      "verified_core",verifiedCore,"split_cover",splitCover,
      "distinct_sources",#coverSources,
      "verified_sources",#seenVerified,
      "simple_sources",simpleSources,"richelot_jacobians",richelotJacobians,
      "exact_tests",exactTests,"hits",#hits,"seconds",Cputime(start);
for H in hits do print "H",H; end for;
quit;
