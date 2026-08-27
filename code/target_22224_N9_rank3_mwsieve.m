//////////////////////////////////////////////////////////////////////
// Rank-3 Mordell--Weil residue sieve for the B=5000 N9 fiber
//
//   (-46,47,49,(-1081/50) T^2).
//
// We use quotient <1,2>, whose minimal elliptic model has rank 3 and
// torsion Z/2.  The quotient imposes the first two square conditions;
// the third is filtered locally together with corrected projective
// direct-contact incidence at p=11,13 and the necessary 3-divisibility
// of the genus-2 Jacobian order at every smooth good reduction.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(18);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);
if not assigned FamilyName then FamilyName:="N9"; end if;
if not assigned PairI then PairI:=1; elif Type(PairI) eq MonStgElt then PairI:=StringToInteger(PairI); end if;
if not assigned PairJ then PairJ:=2; elif Type(PairJ) eq MonStgElt then PairJ:=StringToInteger(PairJ); end if;
if not assigned N then N:=100; elif Type(N) eq MonStgElt then N:=StringToInteger(N); end if;
if not assigned ExactN then ExactN:=N; elif Type(ExactN) eq MonStgElt then ExactN:=StringToInteger(ExactN); end if;
if not assigned PrimeList then PrimeList:=[11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97];
elif Type(PrimeList) eq MonStgElt then PrimeList:=[StringToInteger(s):s in Split(PrimeList,",")]; end if;
if not assigned log_file then log_file:=Sprintf("results/target_22224_%o_rank3_mwsieve_%o%o_N%o.log",FamilyName,PairI,PairJ,N); end if;
if not assigned output_file then output_file:=Sprintf("results/target_22224_%o_rank3_mwsieve_%o%o_N%o.tsv",FamilyName,PairI,PairJ,N); end if;
if not assigned modular_file then modular_file:=Sprintf("results/target_22224_%o_rank3_modular_%o%o_N%o.tsv",FamilyName,PairI,PairJ,N); end if;
if not assigned residue_file then residue_file:=Sprintf("results/target_22224_%o_rank3_residue_classes_%o%o.tsv",FamilyName,PairI,PairJ); end if;
SetLogFile(log_file:Overwrite:=true);

if FamilyName eq "N9" then
    fixed:=[Q!-46,Q!47,Q!49]; d0:=Q!-1081/50; known_t:=Q!5/7;
elif FamilyName eq "S3" then
    fixed:=[Q!-9,Q!35,Q!65]; d0:=Q!-3025/91; known_t:=Q!5/11;
elif FamilyName eq "W1" then
    fixed:=[Q!-2254,Q!-2162,Q!2303]; d0:=Q!4900; known_t:=Q!7/5;
else
    error "unsupported rank-3 family";
end if;
assert PairI in [1..3] and PairJ in [1..3] and PairI lt PairJ;
qi:=PairI; qj:=PairJ;
rr:=[i:i in [1..3]|i notin {qi,qj}]; assert #rr eq 1; remaining:=rr[1];
Rx:=(fixed[qi]+d0*T^2)/(fixed[qi]+d0);
Ry:=(fixed[qj]+d0*T^2)/(fixed[qj]+d0);
C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
Eraw,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi);
E,minmap:=MinimalModel(Eraw); minmapinv:=Inverse(minmap);
gens:=Generators(E); free:=[g:g in gens|Order(g) eq 0]; assert #free eq 3;
TG,tmap:=TorsionSubgroup(E); tors:=[tmap(g):g in TG]; assert Invariants(TG) eq [2];
G1:=free[1]; G2:=free[2]; G3:=free[3];
rawMap:=DefiningPolynomials(minmapinv); curveMap:=DefiningPolynomials(Einv);
projDen:=LCM([Denominator(Q!z):z in fixed cat [d0]]);
projFixed:=[Z!(projDen*z):z in fixed]; projD0:=Z!(projDen*d0);

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
    try return true,EF![BaseRing(EF)!z:z in ints]; catch e return false,EF!0; end try;
end function;
function ReduceMap(polys,F)
    out:=[];
    try for h in polys do Append(~out,ChangeRing(h,F)); end for; return true,out;
    catch e
        allcoeffs:=&cat[Coefficients(h):h in polys]; den:=LCM([Denominator(Q!c):c in allcoeffs]);
        scaled:=[den*h:h in polys]; coeffs:=[Z!c:c in &cat[Coefficients(h):h in scaled]];
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
    if PF eq Parent(PF)!0 then return true,false,F!1; end if;
    ok,raw:=EvalMap(rawF,[PF[i]:i in [1..3]]); if not ok then return true,true,F!0; end if;
    ok,cp:=EvalMap(curveF,raw); if not ok or cp[3] eq 0 then return true,true,F!0; end if;
    return true,false,cp[1]/cp[3];
end function;
function ContactKey(vals,F)
    ell:=Z!Characteristic(F); if &and[F!z eq 0:z in vals] then return false,0; end if;
    best:=ell^4;
    for u in F do if u ne 0 then
        sq:=Sort([Z!((u*(F!z))^2):z in vals]);
        code:=sq[1]*ell^3+sq[2]*ell^2+sq[3]*ell+sq[4]; if code lt best then best:=code; end if;
    end if; end for;
    return true,best;
end function;
function DirectContact(vals,F)
    ell:=Z!Characteristic(F); ok,key:=ContactKey(vals,F); if not ok then return true; end if;
    if ell eq 11 then return key in {1,12,133,136,137,159,1464,1488,1489,1490,1505,1512,1516};
    elif ell eq 13 then return key in {1,14,183,185,186,191,192,194,220,2380,2408,2414,2422,2428,2534};
    end if;
    return true;
end function;
// states: 0 local failure; 1 open/no 3; 2 open/3; 3 direct boundary.
function StateAtT(tt,F)
    ell:=Z!Characteristic(F); P<X>:=PolynomialRing(F);
    denrat:=Q!(fixed[remaining]+d0);
    oka,alpha:=ReduceRat(Q!fixed[remaining]/denrat,F);
    okb,beta:=ReduceRat(Q!d0/denrat,F);
    if not oka or not okb then return 3,-1; end if;
    r3:=alpha+beta*tt^2; if not IsSquare(r3) then return 0,-1; end if;
    vals:=[F!projFixed[k]:k in [1..3]] cat [F!projD0*tt^2];
    if ell in {11,13} and not DirectContact(vals,F) then return 0,-1; end if;
    f:=X*&*[X+z^2:z in vals];
    if Degree(f) ne 5 or Discriminant(f) eq 0 or r3 eq 0 then return 3,-1; end if;
    nj:=Z!Evaluate(LPolynomial(HyperellipticCurve(f)),1);
    return (nj mod 3 eq 0) select 2 else 1,nj;
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
        try fp:=ChangeRing(f,GF(ell)); if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            nj:=Z!Evaluate(LPolynomial(HyperellipticCurve(fp)),1);
            g:=(g eq 0) select nj else GCD(g,nj); Append(~used,<ell,nj,g>);
            if g mod 3 ne 0 then return false,g,used; end if;
        catch e continue; end try;
    end for;
    return g ne 0 and g mod 3 eq 0,g,used;
end function;

print "RANK3_MWSIEVE_START","family",FamilyName,"pair",<qi,qj>,
      "remaining",remaining,"N",N,"E",E,"rank_bounds",RankBounds(E),
      "generators",gens,"known_t",known_t,"primes",PrimeList;
profiles:=[* *]; resout:=Open(residue_file,"w");
fprintf resout,"prime\ttorsion_coset\tm\tn\tk\tord_g1\tord_g2\tord_g3\tstate\tt\n";
for ell in PrimeList do
    F:=GF(ell); ainv:=[]; goodE:=true;
    for q in aInvariants(E) do ok,v:=ReduceRat(q,F); if not ok then goodE:=false; break; end if; Append(~ainv,v); end for;
    if not goodE then print "PRIME_SKIP",ell,"E_denominator"; continue; end if;
    try EF:=EllipticCurve(ainv); catch e print "PRIME_SKIP",ell,"E_bad_reduction"; continue; end try;
    ok1,g1:=ReducePoint(G1,EF); ok2,g2:=ReducePoint(G2,EF); ok3,g3:=ReducePoint(G3,EF);
    if not ok1 or not ok2 or not ok3 then print "PRIME_SKIP",ell,"generator_reduction"; continue; end if;
    o1:=Order(g1); o2:=Order(g2); o3:=Order(g3);
    okr,rawF:=ReduceMap(rawMap,F); okc,curveF:=ReduceMap(curveMap,F);
    if not okr or not okc then print "PRIME_SKIP",ell,"map_reduction"; continue; end if;
    torsF:=[]; okt:=true;
    for tq in tors do ok,tp:=ReducePoint(tq,EF); if not ok then okt:=false; break; end if; Append(~torsF,tp); end for;
    if not okt then print "PRIME_SKIP",ell,"torsion_reduction"; continue; end if;
    stateCache:=[0:i in [0..ell-1]]; for r in [0..ell-1] do st,nj:=StateAtT(F!r,F); stateCache[r+1]:=st; end for;
    dummy:=<Z!0,Z!0,Z!0,Z!0>; allowed:={dummy}; Exclude(~allowed,dummy); counts:=[0,0,0,0];
    mg:=[m*g1:m in [0..o1-1]]; ng:=[n*g2:n in [0..o2-1]]; kg:=[k*g3:k in [0..o3-1]];
    for ti in [1..#torsF] do for m in [0..o1-1] do for n in [0..o2-1] do for k in [0..o3-1] do
        ep:=mg[m+1]+ng[n+1]+kg[k+1]+torsF[ti]; okT,isinf,tt:=TFromFinitePoint(ep,rawF,curveF,F);
        if not okT or isinf then st:=3; ts:="inf"; else st:=stateCache[Z!tt+1]; ts:=Sprint(Z!tt); end if;
        counts[st+1]+:=1;
        if st in {2,3} then Include(~allowed,<ti,m,n,k>);
            fprintf resout,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",ell,ti,m,n,k,o1,o2,o3,st,ts;
        end if;
    end for; end for; end for; end for;
    Append(~profiles,[*ell,o1,o2,o3,allowed*]);
    print "PRIME_PROFILE","p",ell,"orders",o1,o2,o3,"total",&+counts,
          "local_fail",counts[1],"open_no3",counts[2],"open3",counts[3],
          "boundary",counts[4],"allowed",#allowed;
end for;
delete resout;

out:=Open(output_file,"w"); fprintf out,"m\tn\tk\ttorsion_coset\tt\ta\tb\tc\td\tthree_survivor\tgcd_bound\n";
modout:=Open(modular_file,"w"); fprintf modout,"m\tn\tk\ttorsion_coset\n";
tested:=0; modular:=0; exactFull:=0; threeHits:=0; seenT:={Q|};
for m in [-N..N] do for n in [-N..N] do for k in [-N..N] do for ti in [1..#tors] do
    tested+:=1; pass:=true;
    for prof in profiles do
        if <ti,m mod prof[2],n mod prof[3],k mod prof[4]> notin prof[5] then pass:=false; break; end if;
    end for;
    if not pass then continue; end if; modular+:=1; fprintf modout,"%o\t%o\t%o\t%o\n",m,n,k,ti;
    if Max([Abs(m),Abs(n),Abs(k)]) gt ExactN then continue; end if;
    ep:=m*G1+n*G2+k*G3+tors[ti];
    try cp:=Einv(minmapinv(ep)); catch e continue; end try;
    if cp[3] eq 0 then continue; end if;
    tt:=Q!(cp[1]/cp[3]); if tt in seenT then continue; end if; Include(~seenT,tt);
    dd:=d0*tt^2; good:=true;
    for z in fixed do ok,s:=IsSquare((z+dd)/(z+d0)); if not ok then good:=false; break; end if; end for;
    if not good then continue; end if; exactFull+:=1; v:=Primitive(fixed cat [dd]);
    survives,g,used:=ThreeBound(v); if survives then threeHits+:=1; end if;
    fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",m,n,k,ti,tt,v[1],v[2],v[3],v[4],survives select 1 else 0,g;
    print "EXACT_FULL","coeffs",m,n,k,ti,"t",tt,"tuple",v,"three",survives,"gcd",g;
    if survives then
        Pq<X>:=PolynomialRing(Q); f:=X*&*[X+(Q!z)^2:z in v]; G,jmap:=TorsionSubgroup(Jacobian(HyperellipticCurve(f)));
        print "TARGET_EXACT_TORSION",m,n,k,ti,Invariants(G),f;
    end if;
end for; end for; end for; end for;
delete modout; delete out;
print "RANK3_MWSIEVE_DONE","family",FamilyName,"pair",<qi,qj>,
      "tested",tested,"modular",modular,
      "distinct_exact_t",#seenT,"exact_full",exactFull,"three_hits",threeHits,
      "modular_file",modular_file;
UnsetLogFile(); quit;
