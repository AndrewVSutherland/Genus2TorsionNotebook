//////////////////////////////////////////////////////////////////////
// Rigorous finite-field parameter masks for an extra rational
// 3-direction on the rational endpoint-R3 P8 family.
//
// The P8 parameter is tau and the associated b=0 contact-6 source is
//   f=x*(3*x^2+(a-3)*x+2)*(2*x^2-3*x+(a+3)),  a=1/e(tau).
// It already has one rational order-3 line.  A second rational line over
// Q forces dim_F3 J(F_p)[3] >= 2 at every prime of good reduction p != 3.
// Bad/chart-boundary residues are retained, so every rejection is rigorous.
// Richelot is prime to 3, hence the same test applies on the P8 dual.
//
// Usage:
//   magma -b mode:=finite prime_bound:=43 code/..._sieve.m
//   magma -b mode:=search height:=200 prime_bound:=43 code/..._sieve.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode:="finite"; end if;
if not assigned height then height:=200;
elif Type(height) eq MonStgElt then height:=StringToInteger(height); end if;
if not assigned prime_bound then prime_bound:=71;
elif Type(prime_bound) eq MonStgElt then
    prime_bound:=StringToInteger(prime_bound);
end if;
if not assigned max_exact then max_exact:=20;
elif Type(max_exact) eq MonStgElt then max_exact:=StringToInteger(max_exact); end if;

Q:=Rationals(); Z:=Integers(); Qx<x>:=PolynomialRing(Q);
primes:=[p:p in PrimesUpTo(prime_bound)|p notin {2,3,5}];

function ThreeRank(invs)
    return #[n:n in invs|Z!n mod 3 eq 0];
end function;

function P8EFinite(k,u,is_infinity)
    if is_infinity then
        if k!409 eq 0 then return false,k!0; end if;
        return true,-k!200/k!409;
    end if;
    if u^2+k!6 eq 0 then return false,k!0; end if;
    t:=k!4*(u^2+u-k!6)/(u^2+k!6);
    de:=t^4-k!25*t^2+k!1250/k!3;
    if de eq 0 then return false,k!0; end if;
    e:=-(k!25/k!3)*t^2/de;
    if e eq 0 then return false,k!0; end if;
    return true,e;
end function;

function SourceFinite(k,e)
    P<X>:=PolynomialRing(k);
    if e eq 0 then return false,P!0; end if;
    a:=1/e;
    f:=X*(3*X^2+(a-3)*X+2)*(2*X^2-3*X+(a+3));
    if Degree(f) ne 5 or Discriminant(f) eq 0 then return false,f; end if;
    return true,f;
end function;

function ResidueTable(p)
    k:=GF(p); allowed:={Z|}; bad:={Z|}; counts:=AssociativeArray();
    samples:=[];
    // Keys 0,...,p-1 are affine tau; key p is tau=infinity.
    for key in [0..p] do
        inf:=key eq p; u:=inf select k!0 else k!key;
        oke,e:=P8EFinite(k,u,inf);
        if not oke then Include(~bad,key); continue; end if;
        okf,f:=SourceFinite(k,e);
        if not okf then Include(~bad,key); continue; end if;
        A,mp:=AbelianGroup(Jacobian(HyperellipticCurve(f)));
        invs:=Invariants(A); r3:=ThreeRank(invs);
        if IsDefined(counts,r3) then counts[r3]+:=1; else counts[r3]:=1; end if;
        if r3 ge 2 then
            Include(~allowed,key);
            if #samples lt 8 then Append(~samples,<key,invs>); end if;
        end if;
    end for;
    return allowed,bad,counts,samples;
end function;

function RationalParameters(B)
    vals:=[]; seen:={};
    for d in [1..B] do for n in [-B..B] do
        if GCD(n,d) ne 1 then continue; end if;
        q:=Q!n/d; s:=Sprint(q);
        if s in seen then continue; end if;
        Include(~seen,s); Append(~vals,q);
    end for; end for;
    return vals;
end function;

function ProjectiveResidue(q,p)
    k:=GF(p); n:=k!Numerator(q); d:=k!Denominator(q);
    return d eq 0 select p else Z!(n/d);
end function;

function SourceRational(tau)
    t:=4*(tau^2+tau-6)/(tau^2+6);
    de:=t^4-25*t^2+Q!1250/3;
    if de eq 0 then return false,Qx!0,Q!0; end if;
    e:=-(Q!25/3)*t^2/de;
    if e eq 0 then return false,Qx!0,Q!0; end if;
    a:=1/e;
    f:=x*(3*x^2+(a-3)*x+2)*(2*x^2-3*x+(a+3));
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false,f,a;
    end if;
    return true,f,a;
end function;

function IntegralSquareModel(f)
    d:=LCM([Denominator(Coefficient(f,i)):i in [0..Degree(f)]]);
    return Qx!(d^2*f);
end function;

allowed_by_p:=AssociativeArray(); bad_by_p:=AssociativeArray();
print "CONTACT6_M612_P8_EXTRA3_RESIDUE_SIEVE";
print "MODE",mode,"PRIMES",primes;
for p in primes do
    allowed,bad,counts,samples:=ResidueTable(p);
    allowed_by_p[p]:=allowed; bad_by_p[p]:=bad;
    print "PRIME",p,"P1_SIZE",p+1,"ALLOWED",#allowed,"BAD",#bad,
          "RANK_COUNTS",Sort([<r,counts[r]>:r in Keys(counts)]),
          "ALLOWED_KEYS",Sort(Setseq(allowed)),
          "BAD_KEYS",Sort(Setseq(bad)),"SAMPLES",samples;
end for;

if mode eq "search" then
    params:=RationalParameters(height); survivors:=[];
    kills:=AssociativeArray(); bad_uses:=AssociativeArray();
    for tau in params do
        killed:=false;
        for p in primes do
            key:=ProjectiveResidue(tau,p);
            if key in bad_by_p[p] then
                if IsDefined(bad_uses,p) then bad_uses[p]+:=1;
                else bad_uses[p]:=1; end if;
                continue;
            end if;
            if key notin allowed_by_p[p] then
                if IsDefined(kills,p) then kills[p]+:=1; else kills[p]:=1; end if;
                killed:=true; break;
            end if;
        end for;
        if not killed then Append(~survivors,tau); end if;
    end for;
    print "SEARCH_HEIGHT",height,"PARAMETERS",#params,
          "MASK_SURVIVORS",#survivors;
    print "KILL_COUNTS",Sort([<p,kills[p]>:p in Keys(kills)]);
    print "BAD_RESIDUE_USES",Sort([<p,bad_uses[p]>:p in Keys(bad_uses)]);
    print "SURVIVORS",survivors;

    exact:=0; hits:=[]; seen_a:={}; exact_skipped_duplicate:=0;
    for tau in survivors do
        if exact ge max_exact then break; end if;
        ok,f,a:=SourceRational(tau);
        if not ok then continue; end if;
        akey:=Sprint(a);
        if akey in seen_a then exact_skipped_duplicate+:=1; continue; end if;
        Include(~seen_a,akey);
        G,mp:=TorsionSubgroup(Jacobian(HyperellipticCurve(IntegralSquareModel(f))));
        invs:=Invariants(G); exact+:=1;
        print "EXACT",tau,"a",a,"TORSION",invs,"THREE_RANK",ThreeRank(invs);
        if ThreeRank(invs) ge 2 then Append(~hits,<tau,a,invs>); end if;
    end for;
    print "EXACT_UNIQUE_A",exact,"EXACT_DUPLICATE_A_SKIPPED",
          exact_skipped_duplicate,"HITS",hits;
elif mode ne "finite" then
    error "mode must be finite or search";
end if;

print "CONTACT6_M612_P8_EXTRA3_RESIDUE_SIEVE_DONE";
quit;
