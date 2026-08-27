//////////////////////////////////////////////////////////////////////
// Coefficient-lattice MW sieve for the rank-2 pair quotient of one of
// the repeated-triple full-cover fibers.  The quotient enforces two of
// the three square conditions.  At each good prime we retain precisely
// the lattice classes for which the third condition is a square and the
// resulting genus-2 Jacobian can have rational 3-torsion (boundary, or
// 3 divides its finite-field order).  Exact arithmetic is used only on
// the final bounded coefficient classes.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(18);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);

if not assigned Fiber then Fiber:=4; elif Type(Fiber) eq MonStgElt then Fiber:=StringToInteger(Fiber); end if;
if not assigned N then N:=500; elif Type(N) eq MonStgElt then N:=StringToInteger(N); end if;
if not assigned ExactN then ExactN:=N; elif Type(ExactN) eq MonStgElt then ExactN:=StringToInteger(ExactN); end if;
if not assigned PairI then PairI:=1; elif Type(PairI) eq MonStgElt then PairI:=StringToInteger(PairI); end if;
if not assigned PairJ then PairJ:=2; elif Type(PairJ) eq MonStgElt then PairJ:=StringToInteger(PairJ); end if;
if not assigned PrimeList then PrimeList:=[11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97];
elif Type(PrimeList) eq MonStgElt then PrimeList:=[StringToInteger(s):s in Split(PrimeList,",")]; end if;
if not assigned HeightCeiling then HeightCeiling:=10^180;
elif Type(HeightCeiling) eq MonStgElt then HeightCeiling:=StringToInteger(HeightCeiling); end if;

fibers := [
 <"F1",[Q!-1470,-630,336],Q!25,Q!5>,
 <"F2",[Q!-720,20,300],Q!-363,Q!17/11>,
 <"F3",[Q!-612,34,289],Q!-338,Q!25/13>,
 <"F4",[Q!-126,28,49],Q!-50,Q!17/5>,
 <"F5",[Q!-112,14,49],Q!-50,Q!13/5>,
 <"F6",[Q!-50,30,45],Q!-48,Q!2>,
 <"F7",[Q!-18,1,16],Q!-50,Q!29/5>,
 // Repeated fibres first seen in the live B=5000 stream.
 <"N1",[Q!-3,10,18],Q!-48/5,Q!2/5>,
 <"N2",[Q!-25,-1,40],Q!4205/98,Q!49/29>,
 <"N3",[Q!-23,-4,50],Q!578/23,Q!37/17>,
 <"N4",[Q!-36,-3,52],Q!768/13,Q!5/2>,
 <"N5",[Q!-31,-2,64],Q!1250/31,Q!41/25>,
 <"N6",[Q!-3,-2,18],Q!588/169,Q!13/7>,
 <"N7",[Q!-15,16,24],Q!-72/5,Q!3/5>,
 <"N8",[Q!-14,-11,35],Q!1375/98,Q!7/5>,
 <"N9",[Q!-46,47,49],Q!-1081/50,Q!5/7>,
 <"N10",[Q!-45,-12,50],Q!3630/49,Q!49/11>,
 <"N11",[Q!-64,120,169],Q!-1352/15,Q!1/13>,
 <"N12",[Q!-28,50,55],Q!-3610/77,Q!1/19>,
 <"O1",[Q!-1,4,6],Q!-98/27,Q!3/7>,
 <"O2",[Q!-1,6,9],Q!-507/98,Q!7/13>,
 <"O3",[Q!-3,-2,12],Q!162/49,Q!7/3>,
 <"O4",[Q!-9,10,15],Q!-75/8,Q!1/2>,
 <"O5",[Q!-15,-10,16],Q!50/3,Q!5/3>,
 <"O6",[Q!-15,-1,25],Q!605/27,Q!21/11>,
 <"O7",[Q!-9,-5,30],Q!75/8,Q!5/2>,
 <"O8",[Q!-22,55,70],Q!-1372/25,Q!5/7>,
 <"O9",[Q!1,55,99],Q!5/9,Q!15>,
 <"M1",[Q!-16,-5,30],Q!50/3,Q!11/3>,
 <"M2",[Q!-1,25,40],Q!-605/32,Q!17/22>,
 <"M3",[Q!-25,-1,65],Q!405/13,Q!37/27>,
 <"M4",[Q!-49,-2,100],Q!3362/49,Q!61/41>,
 <"M5",[Q!99,244,4026],Q!6,Q!29>,
 // First rank-one quotient from the live B=10000 a~779 snapshot.
 <"R1",[Q!-1,2,25],Q!-529/338,Q!13/23>,
 <"R2",[Q!-44,-16,49],Q!676/11,Q!43/13>,
 <"R3",[Q!-17,-16,50],Q!338/17,Q!53/13>,
 <"R4",[Q!-7,16,56],Q!-392/25,Q!5/13>,
 <"R5",[Q!-23,-18,64],Q!578/23,Q!65/17>,
 <"R6",[Q!-3,5,75],Q!-405/121,Q!11/21>,
 <"R7",[Q!-47,-4,98],Q!2738/47,Q!65/37>,
 // Final three fibers in the completed B=10000, a<=1000 stream.
 <"S1",[Q!-4,9,30],Q!-2166/245,Q!7/19>,
 <"S2",[Q!-64,120,169],Q!-1352/15,Q!1/13>,
 <"S3",[Q!-9,35,65],Q!-3025/91,Q!5/11>,
 <"U1",[Q!-119,124,126],Q!-1054/9,Q!3/7>,
 <"U2",[Q!-1071,-1054,1116],Q!1134,Q!7/3>,
 <"V1",[Q!-18,20,75],Q!-1470/121,Q!11/49>,
 <"V2",[Q!-9,14,49],Q!-2023/162,Q!3/17>,
 <"W1",[Q!-2254,-2162,2303],Q!4900,Q!7/5>,
 <"W2",[Q!-8,9,25],Q!-2209/338,Q!13/47>,
 <"X1",[Q!-49,-21,81],Q!6069/121,Q!143/17>,
 <"Y1",[Q!-15,28,50],Q!-5290/189,Q!3/23>,
 <"Z1",[Q!-169,-120,225],Q!5415/32,Q!26/19>,
 <"AA1",[Q!-64,-40,65],Q!9800/117,Q!9/7>
];
assert Fiber ge 1 and Fiber le #fibers and PairI in [1..3] and PairJ in [1..3] and PairI lt PairJ;
fam:=fibers[Fiber]; name:=fam[1]; fixed:=fam[2]; d0:=fam[3]; known_t:=fam[4];
if not assigned log_file then log_file:=Sprintf("results/target_22224_%o_rank2_mwsieve_%o%o_N%o.log",name,PairI,PairJ,N); end if;
if not assigned output_file then output_file:=Sprintf("results/target_22224_%o_rank2_mwsieve_%o%o_N%o.tsv",name,PairI,PairJ,N); end if;
if not assigned residue_file then residue_file:=Sprintf("results/target_22224_%o_rank2_mwsieve_%o%o_residues.tsv",name,PairI,PairJ); end if;
if not assigned modular_file then modular_file:=Sprintf("results/target_22224_%o_rank2_mwsieve_%o%o_N%o_modular.tsv",name,PairI,PairJ,N); end if;
SetLogFile(log_file:Overwrite:=true);

remaining:=[k:k in [1..3]|k notin {PairI,PairJ}]; assert #remaining eq 1; rem:=remaining[1];
Rx:=(fixed[PairI]+d0*T^2)/(fixed[PairI]+d0);
Ry:=(fixed[PairJ]+d0*T^2)/(fixed[PairJ]+d0);
C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
Eraw,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi);
E,minmap:=MinimalModel(Eraw); minmapinv:=Inverse(minmap);
gens:=Generators(E); free:=[g:g in gens|Order(g) eq 0]; assert #free in {1,2};
TG,tmap:=TorsionSubgroup(E); tors:=[tmap(g):g in TG];
rankone:=#free eq 1; G1:=free[1]; G2:=rankone select E!0 else free[2];
rawMap:=DefiningPolynomials(minmapinv); curveMap:=DefiningPolynomials(Einv);

// Canonical projective root multiset {a^2,b^2,c^2,d^2}.  This removes
// permutations, independent signs, and common scaling of the branches.
function ProjectiveRootKey(vals,ell)
    F:=GF(ell); ss:=[(F!z)^2:z in vals]; choices:=[];
    for i in [1..4] do if ss[i] ne 0 then
        zz:=Sort([Z!(q/ss[i]):q in ss]);
        Append(~choices,<zz[1],zz[2],zz[3],zz[4]>);
    end if; end for;
    if #choices eq 0 then return Sprint(<0,0,0,0>); end if;
    return Sprint(Min(choices));
end function;
function ReadDirectBoundaryMask(ell)
    path:=Sprintf("results/target_22224_direct_contact_deep13_boundary_p%o.tsv",ell);
    lines:=Split(Read(path),"\n"); out:={"__dummy__"}; Exclude(~out,"__dummy__");
    for i in [2..#lines] do if #lines[i] eq 0 then continue; end if;
        z:=Split(lines[i],"\t"); if #z lt 8 then continue; end if;
        vals:=[StringToInteger(z[j]):j in [5..8]];
        Include(~out,ProjectiveRootKey(vals,ell));
    end for;
    return out;
end function;
directMasks:=AssociativeArray();
for ell in [11,13] do directMasks[ell]:=ReadDirectBoundaryMask(ell); end for;

// One integral projective representative of the four branch parameters.
// This matters when d0 has a p-power denominator: coercing d0 and the fixed
// entries separately to GF(p) is invalid, while clearing the one common
// projective denominator gives the correct boundary specialization.
projDen:=LCM([Denominator(Q!z):z in fixed cat [d0]]);
projFixed:=[Z!(projDen*z):z in fixed];
projD0:=Z!(projDen*d0);

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
    out:=[]; try for h in polys do Append(~out,ChangeRing(h,F)); end for;
    catch e return false,[]; end try; return true,out;
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

// 0 = third square fails; 1 = open but no 3-contact; 2 = open 3-contact;
// 3 = boundary or map infinity (retained conservatively).
function StateAtT(tt,F)
    ell:=Z!Characteristic(F); P<X>:=PolynomialRing(F);
    denrat:=Q!(fixed[rem]+d0);
    oka,alpha:=ReduceRat(Q!fixed[rem]/denrat,F);
    okb,beta:=ReduceRat(Q!d0/denrat,F);
    if not oka or not okb then return 3,-1; end if;
    r3:=alpha+beta*tt^2;
    if not IsSquare(r3) then return 0,-1; end if;
    vals:=[F!projFixed[k]:k in [1..3]] cat [F!projD0*tt^2];
    f:=X*&*[X+z^2:z in vals];
    if Degree(f) ne 5 or Discriminant(f) eq 0 or r3 eq 0 then
        // A singular reduction is useful only when its projective root
        // multiset actually occurs on the direct-contact boundary.  Mere
        // root collision is not enough at the priority primes.
        if ell in {11,13} then
            return ProjectiveRootKey(vals,ell) in directMasks[ell] select 3 else 1,-1;
        end if;
        return 3,-1;
    end if;
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
        try fp:=ChangeRing(f,GF(ell));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            nj:=Z!Evaluate(LPolynomial(HyperellipticCurve(fp)),1);
            g:=(g eq 0) select nj else GCD(g,nj); Append(~used,<ell,nj,g>);
            if g mod 3 ne 0 then return false,g,used; end if;
        catch e continue; end try;
    end for;
    return g ne 0 and g mod 3 eq 0,g,used;
end function;

print "REPEATED_RANK2_MWSIEVE_START",name,"pair",<PairI,PairJ>,"remaining",rem,
      "N",N,"fixed",fixed,"d0",d0,"known_t",known_t,
      "E",E,"rank_bounds",RankBounds(E),"generators",gens,"torsion",Invariants(TG),
      "direct_boundary_masks",#directMasks[11],#directMasks[13];
resout:=Open(residue_file,"w");
fprintf resout,"prime\ttorsion_coset\tm\tn\tord_g1\tord_g2\tstate\tt\n";
candidates:=[]; firstProfile:=true; usedPrimes:=[];
for ell in PrimeList do
    F:=GF(ell); ainv:=[]; goodE:=true;
    for q in aInvariants(E) do ok,v:=ReduceRat(q,F); if not ok then goodE:=false; break; end if; Append(~ainv,v); end for;
    if not goodE then print "PRIME_SKIP",ell,"E_denominator"; continue; end if;
    try EF:=EllipticCurve(ainv); catch e print "PRIME_SKIP",ell,"E_bad_reduction"; continue; end try;
    ok1,g1:=ReducePoint(G1,EF); ok2,g2:=ReducePoint(G2,EF);
    if not ok1 or not ok2 then print "PRIME_SKIP",ell,"generator_reduction"; continue; end if;
    o1:=Order(g1); o2:=Order(g2);
    okr,rawF:=ReduceMap(rawMap,F); okc,curveF:=ReduceMap(curveMap,F);
    if not okr or not okc then print "PRIME_SKIP",ell,"map_reduction"; continue; end if;
    torsF:=[]; oktors:=true;
    for tq in tors do ok,tp:=ReducePoint(tq,EF); if not ok then oktors:=false; break; end if; Append(~torsF,tp); end for;
    if not oktors then print "PRIME_SKIP",ell,"torsion_reduction"; continue; end if;
    dummy:=<Z!0,Z!0,Z!0>; allowed:={dummy}; boundary:={dummy}; open3:={dummy};
    Exclude(~allowed,dummy); Exclude(~boundary,dummy); Exclude(~open3,dummy);
    // State depends only on t, not on the (usually many) MW coefficient
    // presentations of the same finite-field point.
    stateCache:=[0:k in [0..ell-1]];
    for k in [0..ell-1] do st0,nj0:=StateAtT(F!k,F); stateCache[k+1]:=st0; end for;
    counts:=[0,0,0,0]; infcount:=0; finiteBoundary:=0;
    for ti in [1..#torsF] do for m in [0..o1-1] do for n in [0..o2-1] do
        ep:=m*g1+n*g2+torsF[ti]; okT,isinf,tt:=TFromFinitePoint(ep,rawF,curveF,F);
        if not okT or isinf then st:=3; ts:="inf"; infcount+:=1;
        else st:=stateCache[Z!tt+1]; ts:=Sprint(Z!tt); if st eq 3 then finiteBoundary+:=1; end if; end if;
        counts[st+1]+:=1; key:=<ti,m,n>;
        if st in {2,3} then Include(~allowed,key); fprintf resout,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",ell,ti,m,n,o1,o2,st,ts; end if;
        if st eq 3 then Include(~boundary,key); elif st eq 2 then Include(~open3,key); end if;
    end for; end for; end for;
    oldCount:=firstProfile select #tors*(2*N+1)*(rankone select 1 else 2*N+1) else #candidates;
    if firstProfile then
        nrange:=rankone select [0] else [-N..N];
        for m in [-N..N] do for n in nrange do for ti in [1..#tors] do
            if <ti,m mod o1,n mod o2> in allowed then Append(~candidates,<m,n,ti>); end if;
        end for; end for; end for;
        firstProfile:=false;
    else
        candidates:=[c:c in candidates|<c[3],c[1] mod o1,c[2] mod o2> in allowed];
    end if;
    Append(~usedPrimes,ell);
    print "PRIME_PROFILE",name,"p",ell,"orders",o1,o2,"total",&+counts,
          "third_fail",counts[1],"open_no3",counts[2],"open_3",counts[3],
          "boundary",counts[4],"finite_direct_boundary",finiteBoundary,"infinity",infcount,
          "allowed",#allowed,"lattice_before",oldCount,
          "lattice_after",#candidates;
    if ell in [11,13] then print "PRIORITY_CLASSES",name,"p",ell,
        "boundary_count",#boundary,"open3_count",#open3,
        "boundary",(#boundary le 100) select Setseq(boundary) else [];
    end if;
    if #candidates eq 0 then break; end if;
end for;
delete resout;

out:=Open(output_file,"w");
fprintf out,"m\tn\ttorsion_coset\tt_num\tt_den\theight\ta\tb\tc\td\tthree_survivor\tgcd_bound\n";
modout:=Open(modular_file,"w"); fprintf modout,"m\tn\ttorsion_coset\n";
exactFull:=0; threeHits:=0; deferred:=0; seenT:={Q|};
for c in candidates do m,n,ti:=Explode(c);
    fprintf modout,"%o\t%o\t%o\n",m,n,ti;
    if Max(Abs(m),Abs(n)) gt ExactN then deferred+:=1; continue; end if;
    ep:=m*G1+n*G2+tors[ti];
    try cp:=Einv(minmapinv(ep)); catch e continue; end try;
    if cp[3] eq 0 then continue; end if;
    tt:=Q!(cp[1]/cp[3]); if tt in seenT then continue; end if; Include(~seenT,tt);
    dd:=d0*tt^2; good:=true;
    for z in fixed do ok,s:=IsSquare((z+dd)/(z+d0)); if not ok then good:=false; break; end if; end for;
    if not good then continue; end if; exactFull+:=1;
    v:=Primitive(fixed cat [dd]); h:=Maximum([Abs(z):z in v]); if h gt HeightCeiling then continue; end if;
    survives,g,used:=ThreeBound(v); if survives then threeHits+:=1; end if;
    fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
            m,n,ti,Numerator(tt),Denominator(tt),h,v[1],v[2],v[3],v[4],survives select 1 else 0,g;
    print "EXACT_FULL",name,"m",m,"n",n,"ti",ti,"t",tt,"tuple",v,
          "three",survives,"gcd",g;
    if survives then
        Pq<X>:=PolynomialRing(Q); f:=X*&*[X+(Q!z)^2:z in v];
        G,jmap:=TorsionSubgroup(Jacobian(HyperellipticCurve(f)));
        print "TARGET_EXACT_TORSION",name,m,n,ti,Invariants(G),f;
    end if;
end for;
delete modout; delete out;
print "REPEATED_RANK2_MWSIEVE_DONE",name,"N",N,"used_primes",usedPrimes,
      "lattice_survivors",#candidates,"distinct_exact_t",#seenT,
      "deferred",deferred,"exact_full",exactFull,"three_hits",threeHits,
      "output",output_file,"modular_file",modular_file;
UnsetLogFile(); quit;
