//////////////////////////////////////////////////////////////////////
// Height-bounded rational a=1/e scan on the exact orthogonal support
// polynomial.  A modular mask retains a residue a only if
//
//       f12_a(v)=0,  M_a(v) is a square
//
// has a solution over F_p.  Coefficient-pole and projective-infinity
// residues are retained as bad, so mask rejection is rigorous.
//
// Normally invoked by contact6_m612_relative3_rational_a_scan.sh.  It can
// also be used as a continuation after the reconstruction variables f12,
// mrec, and urec have been defined.
//////////////////////////////////////////////////////////////////////

if not assigned height then height:=200;
elif Type(height) eq MonStgElt then height:=StringToInteger(height); end if;
if not assigned prime_bound then prime_bound:=43;
elif Type(prime_bound) eq MonStgElt then
    prime_bound:=StringToInteger(prime_bound);
end if;

scan_primes:=[p:p in PrimesUpTo(prime_bound)|p notin {2,3,5}];

function EvalQPolyMod(h,xx,k)
    ans:=k!0;
    for i in Reverse([0..Degree(h)]) do
        ans:=ans*xx+k!Coefficient(h,i);
    end for;
    return ans;
end function;

function EvalRatFunMod(c,xx,k)
    cn:=EvalQPolyMod(Numerator(c),xx,k);
    cd:=EvalQPolyMod(Denominator(c),xx,k);
    if cd eq 0 then return false,k!0; end if;
    return true,cn/cd;
end function;

function SpecializeMod(poly,ee,k,P)
    cs:=[];
    for i in [0..Degree(poly)] do
        ok,c:=EvalRatFunMod(Coefficient(poly,i),ee,k);
        if not ok then return false,P!0; end if;
        Append(~cs,c);
    end for;
    return true,P!cs;
end function;

function ResidueMask(p,f12,mrec)
    k:=GF(p); P<V>:=PolynomialRing(k);
    allowed:={Integers()|}; bad:={Integers()|p};
    zeroM:={Integers()|}; rootCounts:=AssociativeArray();
    // key p is a=infinity.  key 0 means e=infinity; both are retained.
    Include(~bad,0);
    for key in [1..p-1] do
        aa:=k!key; ee:=1/aa;
        okf,ff:=SpecializeMod(f12,ee,k,P);
        okm,mm:=SpecializeMod(mrec,ee,k,P);
        if not (okf and okm) then Include(~bad,key); continue; end if;
        roots:=Roots(ff); rootCounts[key]:=#roots;
        keep:=false;
        for rt in roots do
            mv:=Evaluate(mm,rt[1]);
            if mv eq 0 then
                keep:=true; Include(~zeroM,key);
            elif IsSquare(mv) then
                keep:=true;
            end if;
        end for;
        if keep then Include(~allowed,key); end if;
    end for;
    return allowed,bad,zeroM,rootCounts;
end function;

function RationalParameters(B,Q)
    vals:=[];
    for d in [1..B] do
        for n in [-B..B] do
            if GCD(Abs(n),d) eq 1 then Append(~vals,Q!n/d); end if;
        end for;
    end for;
    return vals;
end function;

function ProjectiveKey(a,p)
    k:=GF(p);
    dn:=k!Denominator(a);
    return dn eq 0 select p else Integers()!((k!Numerator(a))/dn);
end function;

function SpecializeExact(poly,ee,Q,Qz)
    cs:=[];
    for i in [0..Degree(poly)] do
        c:=Coefficient(poly,i);
        den:=Evaluate(Denominator(c),ee);
        if den eq 0 then return false,Qz!0; end if;
        Append(~cs,Evaluate(Numerator(c),ee)/den);
    end for;
    return true,Qz!cs;
end function;

print "CONTACT6_M612_RELATIVE3_RATIONAL_A_SCAN";
print "HEIGHT",height,"PRIMES",scan_primes;
allowedByP:=AssociativeArray(); badByP:=AssociativeArray();
for p in scan_primes do
    allowed,bad,zeroM,rootCounts:=ResidueMask(p,f12,mrec);
    allowedByP[p]:=allowed; badByP[p]:=bad;
    print "MASK_PRIME",p,"P1_SIZE",p+1,"ALLOWED",#allowed,
          "BAD",#bad,"ZERO_M_RESIDUES",#zeroM,
          "ALLOWED_KEYS",Sort(Setseq(allowed)),
          "BAD_KEYS",Sort(Setseq(bad));
end for;

params:=RationalParameters(height,Q); survivors:=[];
kills:=AssociativeArray(); badUses:=AssociativeArray();
for aa in params do
    killed:=false;
    for p in scan_primes do
        key:=ProjectiveKey(aa,p);
        if key in badByP[p] then
            if IsDefined(badUses,p) then badUses[p]+:=1;
            else badUses[p]:=1; end if;
            continue;
        end if;
        if key notin allowedByP[p] then
            if IsDefined(kills,p) then kills[p]+:=1; else kills[p]:=1; end if;
            killed:=true; break;
        end if;
    end for;
    if not killed then Append(~survivors,aa); end if;
end for;

print "SCAN_PARAMETERS",#params,"MASK_SURVIVORS",#survivors;
print "KILL_COUNTS",Sort([<p,kills[p]>:p in Keys(kills)]);
print "BAD_RESIDUE_USES",Sort([<p,badUses[p]>:p in Keys(badUses)]);
print "SURVIVORS",survivors;

Qxscan<Xscan>:=PolynomialRing(Q);
hits:=[]; boundaryHits:=[* *]; exactNoRoot:=0; exactNonsquare:=0;
for aa in survivors do
    if aa eq 0 then
        Append(~boundaryHits,<aa,"a=0/e=infinity">); continue;
    end if;
    ee:=1/aa;
    okf,ff:=SpecializeExact(f12,ee,Q,Qz);
    okm,mm:=SpecializeExact(mrec,ee,Q,Qz);
    oku,uu:=SpecializeExact(urec,ee,Q,Qz);
    if not (okf and okm and oku) then
        Append(~boundaryHits,<aa,"recovery pole">); continue;
    end if;
    rts:=Roots(ff);
    if #rts eq 0 then exactNoRoot+:=1; continue; end if;
    foundSquare:=false;
    for rt in rts do
        vv:=rt[1]; mv:=Evaluate(mm,vv); uv:=Evaluate(uu,vv);
        sq,Lv:=IsSquare(mv);
        if not sq then continue; end if;
        foundSquare:=true;
        source:=Xscan*(3*Xscan^2+(aa-3)*Xscan+2)
                      *(2*Xscan^2-3*Xscan+(aa+3));
        smooth:=Degree(source) eq 5 and Discriminant(source) ne 0;
        open:=mv ne 0 and vv ne 0 and uv^2-4*vv^2 ne 0;
        rec:=<aa,ee,vv,mv,Lv,uv,smooth,open>;
        if smooth and open then Append(~hits,rec);
        else Append(~boundaryHits,<aa,rec>); end if;
    end for;
    if not foundSquare then exactNonsquare+:=1; end if;
end for;

print "EXACT_NO_RATIONAL_SUPPORT_ROOT",exactNoRoot,
      "EXACT_ROOTS_BUT_M_NONSQUARE",exactNonsquare;
print "OPEN_SIGNED_HITS",hits;
print "BOUNDARY_HITS",boundaryHits;
print "CONTACT6_M612_RELATIVE3_RATIONAL_A_SCAN_DONE";
quit;
