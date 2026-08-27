//////////////////////////////////////////////////////////////////////
// Q_p evaluation and p^5 direct-contact tube filter for N9 rank-3
// coefficient survivors.  Run first at p=13, then feed its output to
// p=11.  Scalar multiplication is p-adic, so large MW indices are cheap.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(18);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);
if not assigned FamilyName then FamilyName:="N9"; end if;
if not assigned PairI then PairI:=1; elif Type(PairI) eq MonStgElt then PairI:=StringToInteger(PairI); end if;
if not assigned PairJ then PairJ:=2; elif Type(PairJ) eq MonStgElt then PairJ:=StringToInteger(PairJ); end if;
if not assigned Prime then Prime:=13; elif Type(Prime) eq MonStgElt then Prime:=StringToInteger(Prime); end if;
if not assigned Precision then Precision:=60; elif Type(Precision) eq MonStgElt then Precision:=StringToInteger(Precision); end if;
if not assigned input_file then input_file:="results/target_22224_N9_rank3_periodic_N1000.tsv"; end if;
if not assigned output_file then output_file:=Sprintf("results/target_22224_N9_rank3_periodic_N1000_p%odeep.tsv",Prime); end if;
if not assigned log_file then log_file:=Sprintf("results/target_22224_N9_rank3_periodic_N1000_p%odeep.log",Prime); end if;
assert Prime in {11,13}; SetLogFile(log_file:Overwrite:=true);

if FamilyName eq "N9" then
    fixed:=[Q!-46,Q!47,Q!49]; d0:=Q!-1081/50;
elif FamilyName eq "S3" then
    fixed:=[Q!-9,Q!35,Q!65]; d0:=Q!-3025/91;
else
    error "unsupported rank-3 family";
end if;
assert PairI in [1..3] and PairJ in [1..3] and PairI lt PairJ;
Rx:=(fixed[PairI]+d0*T^2)/(fixed[PairI]+d0);
Ry:=(fixed[PairJ]+d0*T^2)/(fixed[PairJ]+d0);
C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
Eraw,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi);
E,minmap:=MinimalModel(Eraw); minmapinv:=Inverse(minmap);
if FamilyName eq "S3" and PairI eq 1 and PairJ eq 3 then
    // Avoid recomputing the rank-3 Mordell--Weil basis on every p-adic
    // filtering pass.  These are exactly the saturated generators printed
    // by Generators(E) in the profile-generation run, in the same order.
    free:=[
      E![Q!4247,Q!1242150,Q!1],
      E![Q!314353/16,Q!-203211679/64,Q!1],
      E![Q!-59537/9,Q!18091450/27,Q!1]
    ];
    tors:=[E!0,E![Q!-8493,Q!0,Q!1]];
else
    gens:=Generators(E); free:=[g:g in gens|Order(g) eq 0];
    TG,tmap:=TorsionSubgroup(E); tors:=[tmap(g):g in TG];
end if;
assert #free eq 3;
rawMap:=DefiningPolynomials(minmapinv); curveMap:=DefiningPolynomials(Einv);

p:=Prime; mod5:=p^5; K:=pAdicField(p,Precision); EK:=BaseChange(E,K);
function KP(P)
    if P eq Parent(P)!0 then return EK!0; end if;
    return EK![K!P[1],K!P[2],K!P[3]];
end function;
G:=[KP(P):P in free]; torsK:=[KP(P):P in tors];
PK<X,Y,Zc>:=PolynomialRing(K,3); rawK:=[PK!h:h in rawMap]; curveK:=[PK!h:h in curveMap];
function AbsRes(a,m) a:=a mod m; return Min(a,(m-a) mod m); end function;
function DeepKey(vals,p,m)
    best:=<m,m,m,m>;
    for i in [1..4] do x:=vals[i] mod m; if x mod p eq 0 then continue; end if;
        q:=InverseMod(x,m); zz:=Sort([AbsRes((z mod m)*q,m):z in vals]);
        key:=<zz[1],zz[2],zz[3],zz[4]>; if key lt best then best:=key; end if;
    end for;
    return best;
end function;
function ReadBank(path,p,m)
    lines:=Split(Read(path),"\n"); dummy:=<Z!m,Z!m,Z!m,Z!m>; out:={dummy}; Exclude(~out,dummy);
    for i in [2..#lines] do if #lines[i] eq 0 then continue; end if;
        z:=Split(lines[i],"\t"); if #z lt 8 then continue; end if;
        vals:=[StringToInteger(z[j]):j in [5..8]]; Include(~out,DeepKey(vals,p,m));
    end for;
    return out;
end function;
bank:=ReadBank(Sprintf("results/target_22224_direct_contact_deep13_padic_tangent_p%o.tsv",p),p,mod5);
function EvalMap(polys,coords)
    vals:=[Evaluate(h,coords):h in polys]; if &and[z eq 0:z in vals] then return false,vals; end if;
    return true,vals;
end function;
function PadicT(P)
    if P eq EK!0 then return true,false,K!1; end if;
    ok,raw:=EvalMap(rawK,[P[i]:i in [1..3]]); if not ok then return true,true,K!0; end if;
    ok,cp:=EvalMap(curveK,raw); if not ok or cp[3] eq 0 then return true,true,K!0; end if;
    try return true,false,cp[1]/cp[3]; catch err return true,true,K!0; end try;
end function;
function BranchKey(tt,isinf)
    if isinf then return DeepKey([0,0,0,1],p,mod5),-1,-1; end if;
    vals:=[K!z:z in fixed] cat [K!d0*tt^2];
    vv:=[Valuation(z):z in vals|z ne 0]; minv:=Minimum(vv); scale:=K!(p^(-minv));
    ints:=[z eq 0 select Z!0 else (Z!(scale*z)) mod mod5:z in vals];
    return DeepKey(ints,p,mod5),Valuation(tt),(Valuation(tt) ge 0) select (Z!tt) mod mod5 else -1;
end function;

lines:=Split(Read(input_file),"\n"); out:=Open(output_file,"w");
fprintf out,"m\tn\tk\ttorsion_coset\tt_valuation\tt_mod_p5\n";
tested:=0; finite:=0; infinity:=0; deep:=0; unresolved:=0;
print "RANK3_DEEP_PADIC_FILTER_START","family",FamilyName,
      "pair",<PairI,PairJ>,"p",p,"precision",Precision,"bank",#bank,"input",input_file;
for i in [2..#lines] do if #lines[i] eq 0 then continue; end if;
    z:=Split(lines[i],"\t"); if #z lt 4 then continue; end if;
    m:=StringToInteger(z[1]); n:=StringToInteger(z[2]); k:=StringToInteger(z[3]); ti:=StringToInteger(z[4]); tested+:=1;
    P:=m*G[1]+n*G[2]+k*G[3]+torsK[ti]; ok,isinf,tt:=PadicT(P);
    if not ok then unresolved+:=1; continue; end if;
    if isinf or Valuation(tt) lt 0 then infinity+:=1; else finite+:=1; end if;
    key,v,tr:=BranchKey(tt,isinf); if key notin bank then continue; end if;
    deep+:=1; fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\n",m,n,k,ti,v,tr;
end for;
delete out;
print "RANK3_DEEP_PADIC_FILTER_DONE","family",FamilyName,
      "pair",<PairI,PairJ>,"p",p,"tested",tested,"finite",finite,
      "infinity",infinity,"unresolved",unresolved,"deep_p5",deep,"output",output_file;
UnsetLogFile(); quit;
