//////////////////////////////////////////////////////////////////////
// Two-dimensional Mordell--Weil residue sieve for the first elliptic
// quotient of the signed shared-triple F1 fiber
//
//   (528,-726,-891,50*t^2).
//
// The quotient imposing the 528 and -726 square conditions has
//
//   E : y^2 = x^3 + 1520244*x + 202972176,
//   E(Q) = <(7810,-698896),(4741/9,-915607/27)> + <(-132,0)>.
//
// At each good prime we enumerate coefficient pairs (m,n) modulo the
// orders of the two reduced generators.  A class is retained only if the
// third square condition is locally soluble and either the genus-2 curve
// is on the reduction boundary or 3 divides #J(F_p).  We report the
// boundary/open classes at p=11,13, form their generalized-CRT classes,
// and exact-screen a bounded global lattice using every selected prime.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(8);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);
if not assigned N then N:=1000; elif Type(N) eq MonStgElt then N:=StringToInteger(N); end if;
if not assigned PrimeList then PrimeList:=[11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97];
elif Type(PrimeList) eq MonStgElt then PrimeList:=[StringToInteger(s):s in Split(PrimeList,",")]; end if;
if not assigned CRTPrimeList then CRTPrimeList:=[11,13,17,19,23,29,31];
elif Type(CRTPrimeList) eq MonStgElt then CRTPrimeList:=[StringToInteger(s):s in Split(CRTPrimeList,",")]; end if;
if not assigned HeightCeiling then HeightCeiling:=10^160;
elif Type(HeightCeiling) eq MonStgElt then HeightCeiling:=StringToInteger(HeightCeiling); end if;
if not assigned ExactN then ExactN:=N;
elif Type(ExactN) eq MonStgElt then ExactN:=StringToInteger(ExactN); end if;
if not assigned log_file then log_file:=Sprintf("results/target_22224_f1_rank2_mwsieve_N%o.log",N); end if;
if not assigned output_file then output_file:=Sprintf("results/target_22224_f1_rank2_mwsieve_N%o.tsv",N); end if;
if not assigned residue_file then residue_file:="results/target_22224_f1_rank2_residue_classes.tsv"; end if;
if not assigned crt_file then crt_file:="results/target_22224_f1_rank2_periodic_lattice.tsv"; end if;
if not assigned modular_file then modular_file:=Sprintf("results/target_22224_f1_rank2_modular_N%o.tsv",N); end if;
SetLogFile(log_file:Overwrite:=true);

fixed:=[Q!-1470,Q!-630,Q!336]; d0:=Q!25; t1:=Q!5;
// Use the (1,3) quotient.  Unlike (1,2), its t=1 normalization is
// integral at both 11 and 13, so the reduction map retains the projective
// t-coordinate at the two forced-boundary primes.
pair:=[1,3]; other:=2;
Rx:=(fixed[pair[1]]+d0*T^2)/(fixed[pair[1]]+d0);
Ry:=(fixed[pair[2]]+d0*T^2)/(fixed[pair[2]]+d0);
C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
Eraw,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi);
E,minmap:=MinimalModel(Eraw); minmapinv:=Inverse(minmap);
gens:=Generators(E); free:=[g:g in gens|Order(g) eq 0]; assert #free eq 2;
TG,tmap:=TorsionSubgroup(E); tors:=[tmap(g):g in TG];
assert #tors eq 2 and Invariants(TG) eq [2];
G1:=free[1]; G2:=free[2];
rawMap:=DefiningPolynomials(minmapinv); curveMap:=DefiningPolynomials(Einv);

function ReduceRat(q,F)
    q:=Q!q; ell:=Z!Characteristic(F);
    if Denominator(q) mod ell eq 0 then return false,F!0; end if;
    return true,F!Numerator(q)/F!Denominator(q);
end function;

function ReducePoint(PQ,EF)
    if PQ eq Parent(PQ)!0 then return true,EF!0; end if;
    cs:=[Q!PQ[i]:i in [1..3]]; den:=LCM([Denominator(z):z in cs]);
    ints:=[Z!(z*den):z in cs]; g:=GCD([Abs(z):z in ints]);
    if g ne 0 then ints:=[z div g:z in ints]; end if;
    try return true,EF![BaseRing(EF)!z:z in ints];
    catch e return false,EF!0; end try;
end function;

function ReduceMap(polys,F)
    out:=[];
    try
        for h in polys do Append(~out,ChangeRing(h,F)); end for;
        return true,out;
    catch e
        // A rational projective map may still reduce after clearing one
        // common denominator from all coordinates and removing its content.
        // This is essential at p=13 for F1's raw quartic model.
        allcoeffs:=&cat[Coefficients(h):h in polys];
        den:=LCM([Denominator(Q!c):c in allcoeffs]);
        scaled:=[den*h:h in polys];
        coeffs:=[Z!c:c in &cat[Coefficients(h):h in scaled]];
        cont:=GCD([Abs(c):c in coeffs]); if cont eq 0 then return false,[]; end if;
        scaled:=[h/cont:h in scaled]; out:=[];
        try for h in scaled do Append(~out,ChangeRing(h,F)); end for;
        catch e2 return false,[]; end try;
    end try;
    return true,out;
end function;

function EvalMap(polys,coords)
    vals:=[Evaluate(h,coords):h in polys];
    if &and[x eq 0:x in vals] then return false,vals; end if;
    return true,vals;
end function;

function TFromFinitePoint(PF,rawF,curveF,F)
    xyz:=[PF[i]:i in [1..3]];
    ok,raw:=EvalMap(rawF,xyz); if not ok then return true,true,F!0; end if;
    ok,cp:=EvalMap(curveF,raw); if not ok or cp[3] eq 0 then return true,true,F!0; end if;
    return true,false,cp[1]/cp[3];
end function;

function ExactTFromMW(PQ)
    try cp:=Einv(minmapinv(PQ)); catch e return false,true,Q!0; end try;
    if cp[3] eq 0 then return true,true,Q!0; end if;
    return true,false,Q!(cp[1]/cp[3]);
end function;

function ContactKey(vals,F)
    ell:=Z!Characteristic(F);
    if &and[F!z eq 0:z in vals] then return false,0; end if;
    best:=ell^4;
    for u in F do
        if u eq 0 then continue; end if;
        sq:=Sort([Z!((u*(F!z))^2):z in vals]);
        code:=sq[1]*ell^3+sq[2]*ell^2+sq[3]*ell+sq[4];
        if code lt best then best:=code; end if;
    end for;
    return true,best;
end function;

function DirectContact(vals,F)
    ell:=Z!Characteristic(F); ok,key:=ContactKey(vals,F);
    if not ok then return true; end if; // unresolved all-zero blow-up chart
    if ell eq 11 then
        return key in {1,12,133,136,137,159,1464,1488,1489,1490,1505,1512,1516};
    elif ell eq 13 then
        return key in {1,14,183,185,186,191,192,194,220,2380,2408,2414,2422,2428,2534};
    end if;
    return true;
end function;

// state: 0=third condition nonsquare / not a direct-contact boundary;
//        1=open but no 3-contact;
//        2=open and 3-contact; 3=boundary (conservatively retained).
function StateAtT(tt,F)
    ell:=Z!Characteristic(F); P<X>:=PolynomialRing(F);
    vals:=[F!fixed[1],F!fixed[2],F!fixed[3],F!d0*tt^2];
    if ell in {11,13} and not DirectContact(vals,F) then return 0,-1; end if;
    den:=F!(fixed[other]+d0);
    // At p=11 this denominator has valuation two.  The special fibre only
    // supplies the necessary contact class; square solubility is decided by
    // the exact lift after the coefficient sieve.
    if den eq 0 then return 3,-1; end if;
    r3:=(F!fixed[other]+F!d0*tt^2)/den;
    if not IsSquare(r3) then return 0,-1; end if;
    f:=X*&*[X+z^2:z in vals];
    if Degree(f) ne 5 or Discriminant(f) eq 0 or r3 eq 0 then return 3,-1; end if;
    nj:=Z!Evaluate(LPolynomial(HyperellipticCurve(f)),1);
    return (nj mod 3 eq 0) select 2 else 1,nj;
end function;

function CRT2(a,m,b,n)
    g,s,u:=XGCD(m,n); delta:=b-a;
    if delta mod g ne 0 then return false,0,0; end if;
    ng:=n div g; k:=(ng eq 1) select 0 else ((delta div g)*s) mod ng;
    M:=(m div g)*n; return true,(a+m*k) mod M,M;
end function;

function Primitive(vals)
    den:=LCM([Denominator(Q!z):z in vals]); v:=[Z!(den*z):z in vals];
    g:=GCD([Abs(z):z in v]); if g gt 1 then v:=[z div g:z in v]; end if;
    for z in v do if z ne 0 then if z lt 0 then v:=[-w:w in v]; end if; break; end if; end for;
    return v;
end function;

function ThreeBound(v)
    P<X>:=PolynomialRing(Q); f:=X*&*[X+(Q!z)^2:z in v]; g:=0; used:=[];
    for ell in PrimesUpTo(199) do
        try
            fp:=ChangeRing(f,GF(ell));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            nj:=Z!Evaluate(LPolynomial(HyperellipticCurve(fp)),1);
            g:=(g eq 0) select nj else GCD(g,nj); Append(~used,<ell,nj,g>);
            if g mod 3 ne 0 then return false,g,used; end if;
        catch e continue; end try;
    end for;
    return g ne 0 and g mod 3 eq 0,g,used;
end function;

print "F1_RANK2_MWSIEVE_START","N",N,"primes",PrimeList,
      "E",E,"rank_bounds",RankBounds(E),"generators",gens,
      "known_t",t1;
resout:=Open(residue_file,"w");
fprintf resout,"prime\ttorsion_coset\tm\tn\tord_g1\tord_g2\tstate\tt\n";
// Heterogeneous list: indexed-set universes can differ when (for example)
// a prime has no open-3 classes at all.
profiles:=[* *];
for ell in PrimeList do
    F:=GF(ell); ainv:=[]; goodE:=true;
    for q in aInvariants(E) do ok,v:=ReduceRat(q,F); if not ok then goodE:=false; break; end if; Append(~ainv,v); end for;
    if not goodE then print "F1_PRIME_SKIP",ell,"E_denominator"; continue; end if;
    try EF:=EllipticCurve(ainv); catch e print "F1_PRIME_SKIP",ell,"E_bad_reduction"; continue; end try;
    ok1,g1:=ReducePoint(G1,EF); ok2,g2:=ReducePoint(G2,EF);
    if not ok1 or not ok2 then print "F1_PRIME_SKIP",ell,"generator_reduction"; continue; end if;
    o1:=Order(g1); o2:=Order(g2);
    if ell eq 13 then
        // The chosen quartic map collapses in the special fibre.  There are
        // only 3*6*2 coefficient classes, so classify the rational function
        // t on exact small representatives.  A unit denominator makes the
        // t=0 exclusion rigorous on the whole residue disk; poles/basepoints
        // are conservatively retained as the t=infinity boundary chart.
        rawF:=[]; curveF:=[];
    else
        okr,rawF:=ReduceMap(rawMap,F); okc,curveF:=ReduceMap(curveMap,F);
        if not okr or not okc then print "F1_PRIME_SKIP",ell,"map_reduction"; continue; end if;
    end if;
    torsF:=[]; oktors:=true;
    for tq in tors do ok,tp:=ReducePoint(tq,EF); if not ok then oktors:=false; break; end if; Append(~torsF,tp); end for;
    if not oktors then print "F1_PRIME_SKIP",ell,"torsion_reduction"; continue; end if;
    stateCache:=[-1:i in [1..ell]]; orderCache:=[-2:i in [1..ell]];
    for k in [0..ell-1] do st,nj:=StateAtT(F!k,F); stateCache[k+1]:=st; orderCache[k+1]:=nj; end for;
    // Seed then remove a typed dummy so empty state sets retain the same
    // tuple universe (important when a boundary-only prime comes first).
    dummy:=<Z!0,Z!0,Z!0>;
    allowed:={ dummy }; boundary:={ dummy }; open3:={ dummy };
    Exclude(~allowed,dummy); Exclude(~boundary,dummy); Exclude(~open3,dummy);
    counts:=[0,0,0,0];
    for ti in [1..#torsF] do
        for m in [0..o1-1] do for n in [0..o2-1] do
            if ell eq 13 then
                okT,isinf,ttQ:=ExactTFromMW(m*G1+n*G2+tors[ti]);
                if not okT or isinf or Denominator(ttQ) mod ell eq 0 then
                    st:=3; ts:="inf";
                elif Numerator(ttQ) mod ell eq 0 then
                    st:=0; ts:="0";
                else
                    st:=3; ts:=Sprint(Z!(F!Numerator(ttQ)/F!Denominator(ttQ)));
                end if;
            else
                ep:=m*g1+n*g2+torsF[ti]; okT,isinf,tt:=TFromFinitePoint(ep,rawF,curveF,F);
                if not okT or isinf then st:=3; ts:="inf";
                else st:=stateCache[Z!tt+1]; ts:=Sprint(Z!tt); end if;
            end if;
            counts[st+1]+:=1; rec:=<ti,m,n>;
            if st in {2,3} then Include(~allowed,rec); end if;
            if st eq 3 then Include(~boundary,rec); end if;
            if st eq 2 then Include(~open3,rec); end if;
            if ell in CRTPrimeList and st in {2,3} then
                fprintf resout,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",ell,ti,m,n,o1,o2,st,ts;
            end if;
        end for; end for;
    end for;
    Append(~profiles,[*ell,o1,o2,allowed,boundary,open3*]);
    print "F1_PRIME_PROFILE","p",ell,"orders",o1,o2,
          "total",&+counts,"third_fail",counts[1],"open_no3",counts[2],
          "open_3",counts[3],"boundary",counts[4],"allowed",#allowed;
    if ell in CRTPrimeList then
        if #boundary + #open3 le 100 then
            print "F1_PRIORITY_CLASSES","p",ell,
                  "boundary",Setseq(boundary),"open_3",Setseq(open3);
        else
            print "F1_PRIORITY_CLASSES","p",ell,
                  "boundary_count",#boundary,"open_3_count",#open3,
                  "full_list_in",residue_file;
        end if;
    end if;
end for;
delete resout;

// Generalized CRT for the requested priority primes only.
crt:=[<ti,0,1,0,1>:ti in [1..#tors]];
for prof in profiles do
    ell:=prof[1]; o1:=prof[2]; o2:=prof[3];
    allowed:=prof[4]; boundary:=prof[5]; open3:=prof[6];
    if ell notin CRTPrimeList then continue; end if;
    nxt:=[]; seen:={};
    for c in crt do for loc in allowed do
        if c[1] ne loc[1] then continue; end if;
        okm,rm,Mm:=CRT2(c[2],c[3],loc[2],o1);
        okn,rn,Mn:=CRT2(c[4],c[5],loc[3],o2);
        if okm and okn then
            key:=Sprintf("%o:%o:%o:%o:%o",c[1],rm,Mm,rn,Mn);
            if key notin seen then Include(~seen,key); Append(~nxt,<c[1],rm,Mm,rn,Mn>); end if;
        end if;
    end for; end for;
    crt:=nxt; print "F1_CRT_STEP","p",ell,"classes",#crt;
end for;
print "F1_CRT_CLASSES","count",#crt;
if #crt le 500 then print crt; end if;
crtout:=Open(crt_file,"w");
fprintf crtout,"torsion_coset\tm_residue\tm_modulus\tn_residue\tn_modulus\n";
for c in crt do
    fprintf crtout,"%o\t%o\t%o\t%o\t%o\n",c[1],c[2],c[3],c[4],c[5];
end for;
delete crtout;

out:=Open(output_file,"w");
fprintf out,"m\tn\ttorsion_coset\tt_num\tt_den\theight\ta\tb\tc\td\tthree_survivor\tgcd_bound\n";
modout:=Open(modular_file,"w");
fprintf modout,"m\tn\ttorsion_coset\n";
boxTotal:=#tors*(2*N+1)^2; periodic:=0; modular:=0; deferred:=0;
exactFull:=0; heightSkip:=0; threeHits:=0;
// Every class is a rectangular periodic lattice.  Enumerating their points
// avoids scanning the full (2N+1)^2 box.  Conditions at non-CRT primes are
// then only cheap tuple lookups; no elliptic arithmetic occurs before all
// modular tests pass.
for c in crt do
    ti:=c[1]; rm:=c[2]; Mm:=c[3]; rn:=c[4]; Mn:=c[5];
    kmlo:=Ceiling(Q!(-N-rm)/Mm); kmhi:=Floor(Q!(N-rm)/Mm);
    knlo:=Ceiling(Q!(-N-rn)/Mn); knhi:=Floor(Q!(N-rn)/Mn);
    if kmlo gt kmhi or knlo gt knhi then continue; end if;
    for km in [kmlo..kmhi] do
      m:=rm+Mm*km;
      for kn in [knlo..knhi] do
        n:=rn+Mn*kn; periodic+:=1; pass:=true;
        for prof in profiles do
            ell:=prof[1]; if ell in CRTPrimeList then continue; end if;
            o1:=prof[2]; o2:=prof[3]; allowed:=prof[4];
            if <ti,m mod o1,n mod o2> notin allowed then pass:=false; break; end if;
        end for;
        if not pass then continue; end if; modular+:=1;
        fprintf modout,"%o\t%o\t%o\n",m,n,ti;
        // Exact coordinates of nP have height quadratic in n.  Large-index
        // survivors are recorded for a separate p-adic/reconstruction pass.
        if Max(Abs(m),Abs(n)) gt ExactN then deferred+:=1; continue; end if;
    ep:=m*G1+n*G2+tors[ti];
    try cp:=Einv(minmapinv(ep)); catch e continue; end try;
    if cp[3] eq 0 then continue; end if;
    tt:=Q!(cp[1]/cp[3]); dd:=d0*tt^2; good:=true;
    for z in fixed do ok,s:=IsSquare((z+dd)/(z+d0)); if not ok then good:=false; break; end if; end for;
    if not good then continue; end if; exactFull+:=1;
    v:=Primitive(fixed cat [dd]); h:=Maximum([Abs(z):z in v]);
    if h gt HeightCeiling then heightSkip+:=1; continue; end if;
    survives,g,used:=ThreeBound(v); if survives then threeHits+:=1; end if;
    fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
            m,n,ti,Numerator(tt),Denominator(tt),h,v[1],v[2],v[3],v[4],survives select 1 else 0,g;
    print "F1_EXACT_CLASS","m",m,"n",n,"ti",ti,"t",tt,
          "tuple",v,"three",survives,"gcd",g;
    if survives then
        Pq<X>:=PolynomialRing(Q); f:=X*&*[X+(Q!z)^2:z in v];
        G,jmap:=TorsionSubgroup(Jacobian(HyperellipticCurve(f)));
        print "F1_TARGET_EXACT_TORSION",m,n,ti,Invariants(G),f;
    end if;
      end for;
    end for;
end for;
delete modout; delete out;
print "F1_RANK2_MWSIEVE_DONE","box_total",boxTotal,
      "periodic_enumerated",periodic,"modular_survivors",modular,
      "deferred_large_index",deferred,"exact_full",exactFull,
      "height_skip",heightSkip,"three_hits",threeHits,
      "crt_file",crt_file,"modular_file",modular_file;
UnsetLogFile(); quit;
