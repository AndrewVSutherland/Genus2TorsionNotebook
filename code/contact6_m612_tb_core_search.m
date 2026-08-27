//////////////////////////////////////////////////////////////////////
// Direct fiber-product search: low-degree T_B-halving cover plus the
// independent cubic-contact [3,3] core.
//
// Best slicing: fix rational (s,r).  First solve the irreducible cubic
// H(s,r,K)=0, require K=m^2, and recover (a,b).  Then solve the three
// cubic-contact equations as a zero-dimensional ideal in (M,U,v), require
// M=L^2, and verify both marked classes exactly.  This triangular block
// solve is much smaller than fixing (s,v) and solving one four-variable
// high-degree ideal in (r,K,M,U).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned height then height:=10;
elif Type(height) eq MonStgElt then height:=StringToInteger(height); end if;
if not assigned prime_bound then prime_bound:=31;
elif Type(prime_bound) eq MonStgElt then prime_bound:=StringToInteger(prime_bound); end if;
if not assigned progress_interval then progress_interval:=10000;
elif Type(progress_interval) eq MonStgElt then progress_interval:=StringToInteger(progress_interval); end if;
if not assigned max_hits then max_hits:=10;
elif Type(max_hits) eq MonStgElt then max_hits:=StringToInteger(max_hits); end if;

Q:=Rationals(); Z:=Integers();
P<x>:=PolynomialRing(Q);
PF<Tfrob>:=PolynomialRing(Q);

load "code/contact6_m612_tb_core_tools.m";

function Primitive(poly)
    if poly eq 0 then return poly; end if;
    den:=LCM([Denominator(c):c in Coefficients(poly)]);
    vals:=[Z!(den*c):c in Coefficients(poly)];
    cont:=GCD(vals); if cont eq 0 then cont:=1; end if;
    return Parent(poly)!((Q!den/Q!Abs(cont))*poly);
end function;

// Exact cubic H(s,r,K) and the linear recovery a=anum/aden.
RA<aa,ss,rr,kk>:=PolynomialRing(Q,4,"grevlex");
A3A:=aa-3+4*ss^2*rr-2*kk;
H1A:=ss*((aa-3)*rr^2+4*rr-kk*(aa+3))-rr*A3A;
H2A:=8*ss^2*(2*ss^2*rr^2+2*(aa-3)*rr+2
              -kk*(2*ss^2-6))-A3A^2-32*ss^3*rr;
adenA:=Coefficient(H1A,aa,1);
anumA:=-Evaluate(H1A,[RA!0,ss,rr,kk]);
da:=Degree(H2A,aa);
HredA:=&+[Coefficient(H2A,aa,i)*anumA^i*adenA^(da-i):i in [0..da]];
assert IsDivisibleBy(HredA,ss^2*kk);
HcoreA:=ExactQuotient(HredA,ss^2*kk);

RH<s,r,K>:=PolynomialRing(Q,3,"grevlex");
toRH:=hom<RA -> RH | RH!0,s,r,K>;
Hcore:=Primitive(RH!toRH(HcoreA));
aden:=RH!toRH(adenA);
anum:=RH!toRH(anumA);
assert Degree(Hcore,K) eq 3 and IsIrreducible(Hcore);

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

function GoodPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function IntegralModel(f)
    d:=LCM([Denominator(Coefficient(f,i)):i in [0..Degree(f)]]);
    return P!(d^2*f),d;
end function;

function ResidueMask(p)
    kf:=GF(p);
    allowed:={}; boundary:={};
    squares:={z^2:z in kf|z ne 0};
    lcK:=Coefficient(Hcore,K,3);
    for sv in kf do
        for rv in kf do
            key:=<Z!sv,Z!rv>;
            if Evaluate(lcK,[sv,rv,kf!0]) eq 0 then Include(~boundary,key); end if;
            for kv in kf do
                if Evaluate(Hcore,[sv,rv,kv]) ne 0 then continue; end if;
                dv:=Evaluate(aden,[sv,rv,kv]);
                if sv*rv*kv*dv eq 0 then Include(~boundary,key);
                elif kv in squares then Include(~allowed,key); end if;
            end for;
        end for;
    end for;
    return allowed,boundary;
end function;

function ResidueOfRational(q,p)
    if Denominator(q) mod p eq 0 then return false,0; end if;
    kf:=GF(p);
    return true,Z!((kf!Numerator(q))/(kf!Denominator(q)));
end function;

function FrobeniusPolynomial(C,p)
    ef:=EulerFactor(C,p); d:=Degree(ef);
    return &+[Q!Coefficient(ef,i)*Tfrob^(d-i):i in [0..d]];
end function;

function FullSimplicityCertificate(C)
    for p in [q:q in PrimesUpTo(151)|q ge 3] do
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
    return false,"none",0,PF!0;
end function;

function VerifyFullTarget(s0,r0,K0,a0,b0,M0,U0,v0)
    okcore,L,ordD,ordE,E := M612_VerifyCorePoint(a0,b0,M0,U0,v0);
    if not okcore then return false,[],false,0,P!0; end if;
    isSq,m0:=IsSquare(K0);
    if not isSq or m0 eq 0 then return false,[],false,0,P!0; end if;
    f,h6:=ContactPolynomial(a0,b0);
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

maskPrimes:=[p:p in PrimesUpTo(prime_bound)|p notin {2,3}];
allowedByP:=AssociativeArray(); boundaryByP:=AssociativeArray();
print "CONTACT6_M612_TB_CORE_SEARCH";
print "height",height,"prime_bound",prime_bound,"mask_primes",maskPrimes;
for p in maskPrimes do
    al,bd:=ResidueMask(p);
    allowedByP[p]:=al; boundaryByP[p]:=bd;
    print "MASK",p,"allowed",#al,"boundary",#bd;
end for;

params:=RationalParametersOfHeight(height);
UK<t>:=PolynomialRing(Q);
checked:=0; maskSurvivors:=0; cubicRoots:=0; squareLifts:=0;
smooth:=0; coreSolves:=0; coreRationalPoints:=0; coreSquareM:=0;
coreVerified:=0; exactTests:=0; hits:=[]; seen:={};
totalStart:=Cputime();

for s0 in params do
    if s0 eq 0 then continue; end if;
    for r0 in params do
        if r0 eq 0 then continue; end if;
        checked+:=1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress",checked,"mask",maskSurvivors,"roots",cubicRoots,
                  "square",squareLifts,"smooth",smooth,"core_solves",coreSolves,
                  "core_Q",coreRationalPoints,"core_verified",coreVerified,
                  "hits",#hits;
        end if;
        killed:=false;
        for p in maskPrimes do
            oks,sv:=ResidueOfRational(s0,p);
            okr,rv:=ResidueOfRational(r0,p);
            if not oks or not okr then continue; end if;
            key:=<sv,rv>;
            if key notin allowedByP[p] and key notin boundaryByP[p] then
                killed:=true; break;
            end if;
        end for;
        if killed then continue; end if;
        maskSurvivors+:=1;
        sp:=hom<RH -> UK | s0,r0,t>;
        hk:=UK!sp(Hcore);
        if hk eq 0 then continue; end if;
        roots:=Roots(hk); cubicRoots+:=#roots;
        for rt in roots do
            K0:=Q!rt[1];
            if K0 eq 0 then continue; end if;
            sqK,m0:=IsSquare(K0);
            if not sqK or m0 eq 0 then continue; end if;
            squareLifts+:=1;
            den0:=Q!Evaluate(aden,[s0,r0,K0]);
            if den0 eq 0 then continue; end if;
            a0:=Q!Evaluate(anum,[s0,r0,K0])/den0;
            b0:=2*s0^2-3;
            keyab:=Sprint(<a0,b0>);
            if keyab in seen then continue; end if;
            Include(~seen,keyab);
            f,h6:=ContactPolynomial(a0,b0);
            if not GoodPolynomial(f) or Evaluate(h6,1) eq 0 then continue; end if;
            smooth+:=1;

            tcore:=Cputime();
            sols,dim,sat_ok,variety_ok:=M612_FixedCoreSolutions(a0,b0);
            coreSolves+:=1; coreRationalPoints+:=#sols;
            print "TB_CURVE","s",s0,"r",r0,"K",K0,"a",a0,"b",b0,
                  "core_dim",dim,"core_points",#sols,
                  "core_seconds",Cputime(tcore);
            for pt in sols do
                M0:=Q!pt[1]; U0:=Q!pt[2]; v0:=Q!pt[3];
                sqM,L0:=IsSquare(M0);
                if not sqM or M0 eq 0 then continue; end if;
                coreSquareM+:=1;
                okcore,L,ordD,ordE,E:=M612_VerifyCorePoint(a0,b0,M0,U0,v0);
                if not okcore then continue; end if;
                coreVerified+:=1;
                ok,inv,simple,pcert,fI:=
                    VerifyFullTarget(s0,r0,K0,a0,b0,M0,U0,v0);
                if not ok then continue; end if;
                if simple then exactTests+:=1; end if;
                print "CORE_INTERSECTION","s",s0,"r",r0,"K",K0,
                      "a",a0,"b",b0,"M",M0,"L",L,"U",U0,"v",v0,
                      "simple",simple,"pcert",pcert,"torsion",inv;
                if simple and inv eq [6,12] then
                    Append(~hits,<s0,r0,K0,a0,b0,M0,U0,v0,pcert,fI>);
                    print "HIT_6_12",hits[#hits];
                    if #hits ge max_hits then break s0; end if;
                end if;
            end for;
        end for;
    end for;
end for;

print "DONE","checked",checked,"mask_survivors",maskSurvivors,
      "cubic_roots",cubicRoots,"square_lifts",squareLifts,
      "smooth",smooth,"core_solves",coreSolves,
      "core_rational_points",coreRationalPoints,
      "core_square_M",coreSquareM,"core_verified",coreVerified,
      "exact_tests",exactTests,"hits",#hits,
      "seconds",Cputime(totalStart);
for H in hits do print "H",H; end for;
quit;
