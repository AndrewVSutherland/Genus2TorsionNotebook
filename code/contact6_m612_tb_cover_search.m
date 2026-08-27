//////////////////////////////////////////////////////////////////////
// Direct search on the low-degree T_B-halving cover of the contact-6
// family.  This parameterizes the 2-primary target first and lets exact
// torsion detect whether the second independent 3-direction is present.
//
// The generic chart is
//
//   b = 2*s^2-3,  K=m^2,
//
// with a recovered linearly from H1.  After eliminating a and removing the
// open factors s^2*K, the remaining equation H(s,r,K) is irreducible, has
// total degree 8, 31 terms, and degree 3 in K.  For fixed rational (s,r)
// we therefore solve only a cubic in K, require K to be a rational square,
// and reconstruct the rational half exactly.
//
// Run from torsion_jac, for example:
//
//   magma -b height:=10 prime_bound:=23 \
//     code/contact6_m612_tb_cover_search.m \
//     > data/contact6_m612_tb_cover_h10_p23.txt
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned height then height:=10;
elif Type(height) eq MonStgElt then height:=StringToInteger(height); end if;
if not assigned prime_bound then prime_bound:=23;
elif Type(prime_bound) eq MonStgElt then
    prime_bound:=StringToInteger(prime_bound);
end if;
if not assigned max_exact then max_exact:=100;
elif Type(max_exact) eq MonStgElt then
    max_exact:=StringToInteger(max_exact);
end if;
if not assigned max_hits then max_hits:=10;
elif Type(max_hits) eq MonStgElt then
    max_hits:=StringToInteger(max_hits);
end if;
if not assigned progress_interval then progress_interval:=10000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval:=StringToInteger(progress_interval);
end if;
if not assigned diagnostic_halves then diagnostic_halves:=false;
elif Type(diagnostic_halves) eq MonStgElt then
    diagnostic_halves:=diagnostic_halves in {"true","True","1","yes"};
end if;
if not assigned MemGB then MemGB:=8;
elif Type(MemGB) eq MonStgElt then MemGB:=StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

load "../../halving_mumford_tools.m";

Q:=Rationals(); Z:=Integers();
P<x>:=PolynomialRing(Q);
PF<Tfrob>:=PolynomialRing(Q);

function Primitive(poly)
    if poly eq 0 then return poly; end if;
    den:=LCM([Denominator(c):c in Coefficients(poly)]);
    vals:=[Z!(den*c):c in Coefficients(poly)];
    cont:=GCD(vals); if cont eq 0 then cont:=1; end if;
    return Parent(poly)!((Q!den/Q!Abs(cont))*poly);
end function;

// Build H(s,r,K) exactly from H1,H2.
RA<aa,ss,rr,kk>:=PolynomialRing(Q,4,"grevlex");
A3A:=aa-3+4*ss^2*rr-2*kk;
H1A:=ss*((aa-3)*rr^2+4*rr-kk*(aa+3))-rr*A3A;
H2A:=8*ss^2*(2*ss^2*rr^2+2*(aa-3)*rr+2
              -kk*(2*ss^2-6))-A3A^2-32*ss^3*rr;
adenA:=Coefficient(H1A,aa,1);
anumA:=-Evaluate(H1A,[RA!0,ss,rr,kk]);
assert H1A eq adenA*aa-anumA;
da:=Degree(H2A,aa);
HredA:=&+[Coefficient(H2A,aa,i)*anumA^i*adenA^(da-i)
           : i in [0..da]];
assert IsDivisibleBy(HredA,ss^2*kk);
HcoreA:=ExactQuotient(HredA,ss^2*kk);

RH<s,r,K>:=PolynomialRing(Q,3,"grevlex");
toRH:=hom<RA -> RH | RH!0,s,r,K>;
Hcore:=Primitive(RH!toRH(HcoreA));
aden:=RH!toRH(adenA);
anum:=RH!toRH(anumA);
assert Degree(Hcore,K) eq 3;
assert TotalDegree(Hcore) eq 8 and #Terms(Hcore) eq 31;
assert IsIrreducible(Hcore);

function RecoverAValues(s0,r0,K0)
    d0:=Q!Evaluate(aden,[s0,r0,K0]);
    n0:=Q!Evaluate(anum,[s0,r0,K0]);
    if d0 ne 0 then
        return [n0/d0],"generic";
    end if;
    if n0 ne 0 then
        return [],"inconsistent_H1";
    end if;
    UA<z>:=PolynomialRing(Q);
    sp:=hom<RA -> UA | z,s0,r0,K0>;
    h2:=UA!sp(H2A);
    if h2 eq 0 then
        return [],"positive_dimensional_a";
    end if;
    return [Q!rt[1]:rt in Roots(h2)],"exceptional_H1";
end function;

function RationalParametersOfHeight(B)
    vals:=[]; seen:={};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num,den) ne 1 then continue; end if;
            q:=Q!num/den; key:=Sprint(q);
            if key notin seen then
                Include(~seen,key); Append(~vals,q);
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
    L:=1;
    for i in [0..Degree(f)] do
        L:=LCM(L,Denominator(Coefficient(f,i)));
    end for;
    return P!(L^2*f),L;
end function;

function SquareQuarticDataLocal(S)
    if Degree(S) ne 4 then return false,Q!0,P!0; end if;
    s4:=Coefficient(S,4);
    if s4 eq 0 then return false,Q!0,P!0; end if;
    s3:=Coefficient(S,3); s2:=Coefficient(S,2);
    G:=x^2+(s3/(2*s4))*x+(4*s4*s2-s3^2)/(8*s4^2);
    return S eq s4*G^2,s4,G;
end function;

function Has612(invs)
    return #[n:n in invs | (Z!n) mod 6 eq 0] ge 2
       and #[n:n in invs | (Z!n) mod 12 eq 0] ge 1;
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
            if #fac ne 1 or fac[1][2] ne 1
               or Degree(fac[1][1]) ne 4 then
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
            if #fac ne 1 or fac[1][2] ne 1
               or Degree(fac[1][1]) ne 4 then
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

function RationalReduces(q,p)
    return Denominator(q) mod p ne 0;
end function;

function LocalTargetStatus(f,parameters,p)
    // Conservative: a parameter or coefficient denominator divisible by p,
    // or a singular displayed reduction, is a boundary and is allowed to
    // pass.  At a good displayed reduction, [6,12] must inject.
    if &or[not RationalReduces(q,p):q in parameters] then
        return true,"parameter_boundary",[];
    end if;
    if &or[Denominator(Coefficient(f,i)) mod p eq 0
           :i in [0..Degree(f)]] then
        return true,"coefficient_boundary",[];
    end if;
    try
        fp:=ChangeRing(f,GF(p));
        if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
            return true,"curve_boundary",[];
        end if;
        AG,mp:=AbelianGroup(Jacobian(HyperellipticCurve(fp)));
        inv:=Invariants(AG);
        if Has612(inv) then
            return true,"good_allowed",inv;
        end if;
        return false,"good_killed",inv;
    catch e
        return true,"unresolved_boundary",[];
    end try;
end function;

// Cheap polynomial residue masks on (s,r).  A pair is kept if H has a
// nonzero square K root, or if some H-root lies on a recovery boundary.
function ResidueMask(p)
    kf:=GF(p);
    allowed:={}; boundary:={};
    squares:={z^2:z in kf|z ne 0};
    lcK:=Coefficient(Hcore,K,3);
    for sv in kf do
        for rv in kf do
            key:=<Z!sv,Z!rv>;
            // K=infinity in the projective K-line.
            if Evaluate(lcK,[sv,rv,kf!0]) eq 0 then
                Include(~boundary,key);
            end if;
            for kv in kf do
                if Evaluate(Hcore,[sv,rv,kv]) ne 0 then continue; end if;
                dv:=Evaluate(aden,[sv,rv,kv]);
                if sv*rv*kv*dv eq 0 then
                    Include(~boundary,key);
                elif kv in squares then
                    Include(~allowed,key);
                end if;
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

maskPrimes:=[p:p in PrimesUpTo(prime_bound)|p notin {2,3}];
allowedByP:=AssociativeArray(); boundaryByP:=AssociativeArray();

print "CONTACT6_M612_TB_COVER_SEARCH";
print "H_SHAPE",<TotalDegree(Hcore),Degree(Hcore,s),Degree(Hcore,r),
                  Degree(Hcore,K),#Terms(Hcore)>;
print "height",height,"prime_bound",prime_bound,"mask_primes",maskPrimes,
      "max_exact",max_exact,"max_hits",max_hits,
      "diagnostic_halves",diagnostic_halves;
for p in maskPrimes do
    al,bd:=ResidueMask(p);
    allowedByP[p]:=al; boundaryByP[p]:=bd;
    print "MASK",p,"allowed",#al,"boundary",#bd;
end for;

params:=RationalParametersOfHeight(height);
UK<t>:=PolynomialRing(Q);
checked:=0; maskSurvivors:=0; cubicRoots:=0; squareLifts:=0;
smooth:=0; localSurvivors:=0; simpleSurvivors:=0; exactTests:=0;
hits:=[]; seen:={};
killCounts:=AssociativeArray();

for s0 in params do
    if s0 eq 0 then continue; end if;
    for r0 in params do
        if r0 eq 0 then continue; end if;
        checked+:=1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress",checked,"mask",maskSurvivors,
                  "roots",cubicRoots,"square",squareLifts,
                  "smooth",smooth,"local",localSurvivors,
                  "simple",simpleSurvivors,"exact",exactTests,
                  "hits",#hits;
        end if;

        killed:=false;
        for p in maskPrimes do
            oks,sv:=ResidueOfRational(s0,p);
            okr,rv:=ResidueOfRational(r0,p);
            if not oks or not okr then continue; end if;
            key:=<sv,rv>;
            if key notin allowedByP[p] and key notin boundaryByP[p] then
                killed:=true;
                if IsDefined(killCounts,p) then killCounts[p]+:=1;
                else killCounts[p]:=1; end if;
                break;
            end if;
        end for;
        if killed then continue; end if;
        maskSurvivors+:=1;

        sp:=hom<RH -> UK | s0,r0,t>;
        hk:=UK!sp(Hcore);
        if hk eq 0 then
            print "IDENTICALLY_ZERO_CUBIC","s",s0,"r",r0;
            continue;
        end if;
        roots:=Roots(hk);
        cubicRoots+:=#roots;
        for rt in roots do
            K0:=Q!rt[1];
            if K0 eq 0 then continue; end if;
            isSq,m0:=IsSquare(K0);
            if not isSq or m0 eq 0 then continue; end if;
            squareLifts+:=1;
            avals,recoveryStatus:=RecoverAValues(s0,r0,K0);
            if recoveryStatus ne "generic" then
                print "RECOVERY_BOUNDARY","s",s0,"r",r0,"K",K0,
                      "status",recoveryStatus,"rational_a",avals;
            end if;
            if #avals eq 0 then continue; end if;
            for a0 in avals do
            b0:=2*s0^2-3;
            keyab:=Sprint(<a0,b0>);
            if keyab in seen then continue; end if;
            Include(~seen,keyab);

            f,h6:=ContactPolynomial(a0,b0);
            if not GoodPolynomial(f) or Evaluate(h6,Q!1) eq 0 then
                continue;
            end if;
            smooth+:=1;

            // Polynomial-level half reconstruction is cheap and is checked
            // before any local gate.  Optional Jacobian diagnostics verify
            // the marked relation and the order-4 half even for candidates
            // which are immediately killed at a good prime.
            qB:=x^2+((a0-3)/(b0+3))*x+2/(b0+3);
            A30:=a0-3+4*s0^2*r0-2*K0;
            p0:=A30/(4*s0^2);
            Ghalf:=x^2+p0*x+r0/s0;
            Mtool:=2*s0^2/m0;
            Ntool:=2*s0^2*r0/m0;
            ell:=qB*(Mtool*x+Ntool);
            Sq:=ExactQuotient(ell^2-f,qB);
            sqOK,lambda,Gcheck:=SquareQuarticDataLocal(Sq);
            polyHalfOK:=sqOK and Gcheck eq Ghalf;
            assert polyHalfOK;

            jacHalfOK:="not_run";
            if diagnostic_halves then
                Cdiag:=HyperellipticCurve(f);
                Jdiag:=Jacobian(Cdiag);
                Ddiag:=Jdiag![x-1,Evaluate(h6,Q!1)];
                TBdiag:=Jdiag![qB,P!0];
                Hdiag:=Jdiag![Ghalf,-ell mod Ghalf];
                if Order(Ddiag) eq 6 and 3*Ddiag eq TBdiag
                   and 2*Hdiag eq TBdiag and Order(Hdiag) eq 4 then
                    jacHalfOK:="true";
                else
                    jacHalfOK:="false";
                end if;
            end if;

            // The invariant lists can have different lengths at different
            // primes, so keep the diagnostic tuples in an untyped list.
            statuses:=[* *];
            localOK:=true;
            firstKill:=0;
            localParameters:=[s0,r0,K0,a0,b0];
            for p in maskPrimes do
                ok,status,invp:=LocalTargetStatus(f,localParameters,p);
                Append(~statuses,<p,status,invp>);
                if not ok then
                    localOK:=false; firstKill:=p; break;
                end if;
            end for;
            print "SMOOTH_DIAGNOSTIC","s",s0,"r",r0,"K",K0,
                  "a",a0,"b",b0,"poly_half",polyHalfOK,
                  "jac_half",jacHalfOK,"first_kill",firstKill,
                  "statuses",statuses;
            if not localOK then continue; end if;
            localSurvivors+:=1;
            print "LOCAL_SURVIVOR","s",s0,"r",r0,"K",K0,
                  "a",a0,"b",b0,"statuses",statuses;

            // Exact marked order-6 and rational-half reconstruction.
            C:=HyperellipticCurve(f); J:=Jacobian(C);
            D6:=J![x-1,Evaluate(h6,Q!1)];
            if Order(D6) ne 6 then continue; end if;
            TB:=J![qB,P!0];
            if 3*D6 ne TB then continue; end if;
            if not sqOK or Gcheck ne Ghalf then
                print "HALF_RECONSTRUCTION_FAILURE",s0,r0,K0;
                continue;
            end if;
            H4:=J![Ghalf,-ell mod Ghalf];
            if 2*H4 ne TB or Order(H4) ne 4 then
                print "HALF_JACOBIAN_FAILURE",s0,r0,K0;
                continue;
            end if;
            D12:=H4+2*D6;
            assert 2*D12 eq D6 and Order(D12) eq 12;
            TA:=J![x,P!0];
            assert TA ne TB and 2*TA eq J!0;

            fI,scale:=IntegralModel(f);
            CI:=HyperellipticCurve(fI);
            simple,method,pcert,Phi:=FullSimplicityCertificate(CI);
            if not simple then
                print "NONSIMPLE_OR_UNCERTIFIED","a",a0,"b",b0;
                continue;
            end if;
            simpleSurvivors+:=1;
            print "SIMPLE_SURVIVOR","a",a0,"b",b0,
                  "method",method,"pcert",pcert,"Phi",Phi;
            if exactTests ge max_exact then continue; end if;
            exactTests+:=1;
            JI:=Jacobian(CI);
            TG,mp:=TorsionSubgroup(JI);
            inv:=Invariants(TG);
            print "EXACT","a",a0,"b",b0,"torsion",inv,"f",fI;
            if inv eq [6,12] then
                Append(~hits,<s0,r0,K0,a0,b0,method,pcert,fI>);
                print "HIT_6_12","s",s0,"r",r0,"K",K0,
                      "a",a0,"b",b0,"method",method,
                      "pcert",pcert,"f",fI;
                if #hits ge max_hits then break s0; end if;
            end if;
            end for;
        end for;
    end for;
end for;

print "DONE","checked",checked,"mask_survivors",maskSurvivors,
      "cubic_roots",cubicRoots,"square_lifts",squareLifts,
      "smooth",smooth,"local_survivors",localSurvivors,
      "simple_survivors",simpleSurvivors,"exact_tests",exactTests,
      "hits",#hits;
print "KILL_COUNTS",Sort([<p,killCounts[p]>:p in Keys(killCounts)]);
for hit in hits do print "HIT_RECORD",hit; end for;
quit;
