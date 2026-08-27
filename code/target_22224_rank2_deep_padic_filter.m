//////////////////////////////////////////////////////////////////////
// Generic Q_p / p^5 tube filter for rank-2 repeated-fiber survivors.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(18);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);
if not assigned FamilyName then FamilyName:="M5"; end if;
if not assigned Prime then Prime:=13; elif Type(Prime) eq MonStgElt then Prime:=StringToInteger(Prime); end if;
if not assigned Precision then Precision:=60; elif Type(Precision) eq MonStgElt then Precision:=StringToInteger(Precision); end if;
if FamilyName eq "M5" then
    fixed:=[Q!99,Q!244,Q!4026]; d0:=Q!6; qi:=1; qj:=2;
elif FamilyName eq "P8" then
    fixed:=[Q!528,Q!-726,Q!-891]; d0:=Q!50; qi:=1; qj:=2;
elif FamilyName eq "N2" then
    fixed:=[Q!-25,Q!-1,Q!40]; d0:=Q!4205/98; qi:=1; qj:=3;
elif FamilyName eq "O5" then
    fixed:=[Q!-15,Q!-10,Q!16]; d0:=Q!50/3; qi:=1; qj:=2;
elif FamilyName eq "O6" then
    fixed:=[Q!-15,Q!-1,Q!25]; d0:=Q!605/27; qi:=1; qj:=2;
elif FamilyName eq "R1" then
    fixed:=[Q!-1,Q!2,Q!25]; d0:=Q!-529/338; qi:=1; qj:=2;
elif FamilyName eq "R6" then
    // B=10000 repeated fiber (-3,5,75,-405/121*T^2).
    // Use the <1,3> quotient, of rank 2 with torsion Z/2.
    fixed:=[Q!-3,Q!5,Q!75]; d0:=Q!-405/121; qi:=1; qj:=3;
elif FamilyName eq "S1" then
    // Final B=10000 repeated fiber with a rank-one <1,2> quotient.
    fixed:=[Q!-4,Q!9,Q!30]; d0:=Q!-2166/245; qi:=1; qj:=2;
elif FamilyName eq "U2" then
    fixed:=[Q!-1071,Q!-1054,Q!1116]; d0:=Q!1134; qi:=1; qj:=3;
else error "unsupported family"; end if;
if not assigned input_file then input_file:=Sprintf("results/target_22224_%o_rank2_modular_N10000.tsv",FamilyName); end if;
if not assigned output_file then output_file:=Sprintf("results/target_22224_%o_rank2_modular_N10000_p%odeep.tsv",FamilyName,Prime); end if;
if not assigned log_file then log_file:=Sprintf("results/target_22224_%o_rank2_modular_N10000_p%odeep.log",FamilyName,Prime); end if;
assert Prime in {11,13}; SetLogFile(log_file:Overwrite:=true);

Rx:=(fixed[qi]+d0*T^2)/(fixed[qi]+d0); Ry:=(fixed[qj]+d0*T^2)/(fixed[qj]+d0);
C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
Eraw,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi); E,minmap:=MinimalModel(Eraw); minmapinv:=Inverse(minmap);
if FamilyName eq "U2" then
    // Match the saturated basis used by target_22224_p8_rank2_mwsieve.m.
    // Magma can return a different unimodular basis in this model when the
    // surrounding call sequence differs, which would invalidate imported
    // coefficient triples even though both bases generate E(Q).
    free:=[E![Q!1543/121,Q!8819100/1331,Q!1],
           E![Q!2008,Q!-105300,Q!1]];
    tors:=[E!0,E![Q!-17,Q!0,Q!1]];
    gens:=tors[2..2] cat free;
else
    gens:=Generators(E); free:=[g:g in gens|Order(g) eq 0];
    TG,tmap:=TorsionSubgroup(E); tors:=[tmap(g):g in TG];
end if;
assert #free in {1,2};
rawMap:=DefiningPolynomials(minmapinv); curveMap:=DefiningPolynomials(Einv);
p:=Prime; mod5:=p^5; K:=pAdicField(p,Precision); EK:=BaseChange(E,K);
function KP(P)
    if P eq Parent(P)!0 then return EK!0; end if; return EK![K!P[1],K!P[2],K!P[3]];
end function;
G:=[KP(P):P in free]; torsK:=[KP(P):P in tors];
PK<X,Y,Zc>:=PolynomialRing(K,3); rawK:=[PK!h:h in rawMap]; curveK:=[PK!h:h in curveMap];
function AbsRes(a,m) a:=a mod m; return Min(a,(m-a) mod m); end function;
function DeepKey(vals,p,m)
    best:=<m,m,m,m>;
    for i in [1..4] do x:=vals[i] mod m; if x mod p eq 0 then continue; end if;
        q:=InverseMod(x,m); zz:=Sort([AbsRes((z mod m)*q,m):z in vals]);
        key:=<zz[1],zz[2],zz[3],zz[4]>; if key lt best then best:=key; end if;
    end for; return best;
end function;
function ReadBank(path,p,m)
    lines:=Split(Read(path),"\n"); dummy:=<Z!m,Z!m,Z!m,Z!m>; out:={dummy}; Exclude(~out,dummy);
    for i in [2..#lines] do if #lines[i] eq 0 then continue; end if;
        z:=Split(lines[i],"\t"); if #z lt 8 then continue; end if;
        Include(~out,DeepKey([StringToInteger(z[j]):j in [5..8]],p,m));
    end for; return out;
end function;
bank:=ReadBank(Sprintf("results/target_22224_direct_contact_deep13_padic_tangent_p%o.tsv",p),p,mod5);
function EvalMap(polys,coords)
    vals:=[Evaluate(h,coords):h in polys]; if &and[z eq 0:z in vals] then return false,vals; end if; return true,vals;
end function;
function PadicT(P)
    if P eq EK!0 then return true,false,K!1; end if;
    ok,raw:=EvalMap(rawK,[P[i]:i in [1..3]]); if not ok then return true,true,K!0; end if;
    ok,cp:=EvalMap(curveK,raw); if not ok or cp[3] eq 0 then return true,true,K!0; end if;
    try return true,false,cp[1]/cp[3]; catch err return true,true,K!0; end try;
end function;
function BranchKey(tt,isinf)
    if isinf then return DeepKey([0,0,0,1],p,mod5),-1,-1; end if;
    vals:=[K!z:z in fixed] cat [K!d0*tt^2]; vv:=[Valuation(z):z in vals|z ne 0];
    minv:=Minimum(vv); scale:=K!(p^(-minv)); ints:=[z eq 0 select Z!0 else (Z!(scale*z)) mod mod5:z in vals];
    return DeepKey(ints,p,mod5),Valuation(tt),(Valuation(tt) ge 0) select (Z!tt) mod mod5 else -1;
end function;

lines:=Split(Read(input_file),"\n"); out:=Open(output_file,"w");
fprintf out,"m\tn\ttorsion_coset\tt_valuation\tt_mod_p5\n";
tested:=0; finite:=0; infinity:=0; deep:=0; unresolved:=0;
print "RANK2_DEEP_PADIC_START",FamilyName,"p",p,"precision",Precision,
      "bank",#bank,"generators",gens,"input",input_file;
for i in [2..#lines] do if #lines[i] eq 0 then continue; end if;
    z:=Split(lines[i],"\t"); if #z lt 3 then continue; end if;
    m:=StringToInteger(z[1]); n:=StringToInteger(z[2]); ti:=StringToInteger(z[3]); tested+:=1;
    P:=m*G[1]+torsK[ti]; if #G eq 2 then P+:=n*G[2]; end if;
    ok,isinf,tt:=PadicT(P); if not ok then unresolved+:=1; continue; end if;
    if isinf or Valuation(tt) lt 0 then infinity+:=1; else finite+:=1; end if;
    key,v,tr:=BranchKey(tt,isinf); if key notin bank then continue; end if;
    deep+:=1; fprintf out,"%o\t%o\t%o\t%o\t%o\n",m,n,ti,v,tr;
end for;
delete out;
print "RANK2_DEEP_PADIC_DONE",FamilyName,"p",p,"tested",tested,"finite",finite,
      "infinity",infinity,"unresolved",unresolved,"deep_p5",deep,"output",output_file;
UnsetLogFile(); quit;
