//////////////////////////////////////////////////////////////////////
// Intersect the matching-basis V1 N=10000 MW survivors with the intrinsic
// p=13 full-cover/contact components.  The basis is fixed explicitly to
// match target_22224_V1_rank2_modular_N10000.tsv.
//////////////////////////////////////////////////////////////////////
// Match the map/basis construction in target_22224_p8_rank2_mwsieve.m.
// The elliptic curve model alone is not enough: a different random seed can
// compose the inverse quartic map with a nontrivial elliptic automorphism.
SetColumns(0); SetSeed(8);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);
p:=13; depth:=5; modulus:=p^depth; DeepDepth:=12; deepmod:=p^DeepDepth;
Precision:=80;
K:=pAdicField(p,Precision);
input_file:="results/target_22224_V1_rank2_modular_N10000.tsv";
output_file:="results/target_22224_V1_MW_p13_components.tsv";
log_file:="results/target_22224_V1_MW_p13_components.log";
SetLogFile(log_file:Overwrite:=true);

fixed:=[Q!-18,Q!20,Q!75]; d0:=Q!-1470/121;
Rx:=(fixed[1]+d0*T^2)/(fixed[1]+d0);
Ry:=(fixed[2]+d0*T^2)/(fixed[2]+d0);
C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
Eraw,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi);
E,minmap:=MinimalModel(Eraw); minmapinv:=Inverse(minmap);
assert aInvariants(E) eq [Q!0,Q!1,Q!0,Q!35967,Q!252063];
free:=[E![Q!3603,Q!216600,Q!1],E![Q!93,Q!-2100,Q!1]];
tors:=[E!0,E![Q!-7,Q!0,Q!1]];

EK:=BaseChange(E,K);
function KP(P)
    if P eq Parent(P)!0 then return EK!0; end if;
    return EK![K!P[1],K!P[2],K!P[3]];
end function;
G:=[KP(P):P in free]; torsK:=[KP(P):P in tors];
PK<X,Y,Zc>:=PolynomialRing(K,3);
rawK:=[PK!h:h in DefiningPolynomials(minmapinv)];
curveK:=[PK!h:h in DefiningPolynomials(Einv)];

function EvalMap(polys,coords)
    vals:=[Evaluate(h,coords):h in polys];
    if &and[z eq 0:z in vals] then return false,vals; end if;
    return true,vals;
end function;

function PadicT(P)
    if P eq EK!0 then return true,false,K!1,"identity"; end if;
    ok,raw:=EvalMap(rawK,[P[i]:i in [1..3]]);
    if not ok then return false,true,K!0,"raw_basepoint"; end if;
    ok,cp:=EvalMap(curveK,raw);
    if not ok then return false,true,K!0,"curve_basepoint"; end if;
    if cp[3] eq 0 then return true,true,K!0,"projective_pole"; end if;
    try return true,false,cp[1]/cp[3],"finite";
    catch err return false,true,K!0,"division_error"; end try;
end function;

function CoverPattern(tt)
    // Integral projective scaling of [-18,20,75,(-1470/121)T^2].
    a:=K!-2178; b:=K!2420; c:=K!9075; d:=K!-1470*tt^2;
    rr:=[a*b*c*d,
         a*(a+b)*(a+c)*(a+d),
         b*(b+a)*(b+c)*(b+d),
         c*(c+a)*(c+b)*(c+d)];
    return [IsSquare(q):q in rr],[Valuation(q):q in rr];
end function;

// The complete viable depth-five mask derived independently by
// target_22224_V1_depth5_Tkeys.py.  The exceptional strict transform is a
// formal graph over all four residue disks.  The +/-1 key restriction is
// precisely the ramified full-cover square condition, tested by CoverPattern.
function ComponentAtKey(tkey,full)
    if not full then return false,"none"; end if;
    r:=tkey mod p;
    if r eq 5 then return true,"exceptional_plus5"; end if;
    if r eq 8 then return true,"exceptional_minus5"; end if;
    if r eq 1 then return true,"exceptional_plus1"; end if;
    if r eq 12 then return true,"exceptional_minus1"; end if;
    return false,"none";
end function;

lines:=Split(Read(input_file),"\n"); out:=Open(output_file,"w");
fprintf out,"m\tn\ttorsion_coset\tmap_status\tv13_T\tT_mod13\tT_mod13pow5\tT_mod13pow12\tcomponent\tcover_pattern\tcover_valuations\n";
tested:=0; finite:=0; nonintegral:=0; maploss:=0;
coverfull:=0; componenthits:=0; component_nonfull:=0;
print "V1_MW_P13_COMPONENTS_START","input",input_file,
      "precision",Precision,"basis",free,"torsion",tors;
for i in [2..#lines] do
    if #lines[i] eq 0 then continue; end if;
    w:=Split(lines[i],"\t"); if #w lt 3 then continue; end if;
    m:=StringToInteger(w[1]); n:=StringToInteger(w[2]); ti:=StringToInteger(w[3]);
    tested+:=1; P:=m*G[1]+n*G[2]+torsK[ti];
    ok,isinf,tt,status:=PadicT(P);
    if not ok then
        // Removable p-adic map presentations are recovered exactly when
        // the coefficient height is small enough (the current file has one).
        recovered:=false;
        if Max(Abs(m),Abs(n)) le 20 then
            try
                cpQ:=Einv(minmapinv(m*free[1]+n*free[2]+tors[ti]));
                if cpQ[3] ne 0 then
                    tq:=Q!(cpQ[1]/cpQ[3]);
                    tt:=K!tq; ok:=true; isinf:=false;
                    status:="exact_recovery"; recovered:=true;
                    print "V1_MW_EXACT_RECOVERY",m,n,ti,"T",tq;
                end if;
            catch err
                recovered:=false;
            end try;
        end if;
        if not recovered then
            maploss+:=1;
            fprintf out,"%o\t%o\t%o\t%o\tunresolved\t-1\t-1\t-1\tunresolved\tunresolved\tunresolved\n",
                    m,n,ti,status;
            print "V1_MW_MAPLOSS",m,n,ti,status;
            continue;
        end if;
    end if;
    if isinf or Valuation(tt) lt 0 then
        nonintegral+:=1;
        fprintf out,"%o\t%o\t%o\t%o\t%o\t-1\t-1\t-1\tinfinity_killed\tfalse\tna\n",
                m,n,ti,status,isinf select -1000000 else Valuation(tt);
        continue;
    end if;
    finite+:=1;
    // Precision is far beyond five digits for every finite row retained.
    assert AbsolutePrecision(tt) ge depth;
    tkey:=(Z!tt) mod modulus; deepkey:=(Z!tt) mod deepmod;
    tmod:=tkey mod p;
    pat,vals:=CoverPattern(tt); full:=&and pat;
    if full then coverfull+:=1; end if;
    hit,component:=ComponentAtKey(tkey,full);
    if hit and full then componenthits+:=1;
    elif hit then component_nonfull+:=1; end if;
    fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
            m,n,ti,status,Valuation(tt),tmod,tkey,deepkey,component,pat,vals;
    print "V1_MW_ROW","m",m,"n",n,"ti",ti,"tmod13",tmod,
          "tkey",tkey,"component",component,"cover",pat,
          "valuations",vals;
end for;
delete out;
print "V1_MW_P13_COMPONENTS_DONE","tested",tested,"finite",finite,
      "nonintegral",nonintegral,"maploss",maploss,"coverfull",coverfull,
      "component_hits",componenthits,"component_but_nonfull",component_nonfull,
      "output",output_file;
UnsetLogFile(); quit;
