//////////////////////////////////////////////////////////////////////
// p=13 bad-reduction filter for the deep N2 coefficient-lattice
// survivors.  Scalar multiplication is performed over Q_13, so even
// coefficient indices around 10^4 remain cheap.  The projective t-line
// contact set is {0,+-1,+-5,infinity} modulo 13.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(18);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);
if not assigned Precision then Precision:=20; elif Type(Precision) eq MonStgElt then Precision:=StringToInteger(Precision); end if;
if not assigned input_file then input_file:="results/target_22224_N2_rank2_modular_N10000.tsv"; end if;
if not assigned output_file then output_file:="results/target_22224_N2_rank2_modular_N10000_p13.tsv"; end if;
if not assigned log_file then log_file:="results/target_22224_N2_rank2_modular_N10000_p13.log"; end if;
SetLogFile(log_file:Overwrite:=true);

fixed:=[Q!-25,Q!-1,Q!40]; d0:=Q!4205/98;
Rx:=(fixed[1]+d0*T^2)/(fixed[1]+d0);
Ry:=(fixed[3]+d0*T^2)/(fixed[3]+d0);
C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
Eraw,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi);
E,minmap:=MinimalModel(Eraw); minmapinv:=Inverse(minmap);
gens:=Generators(E); free:=[g:g in gens|Order(g) eq 0]; assert #free eq 2;
TG,tmap:=TorsionSubgroup(E); tors:=[tmap(g):g in TG];
rawMap:=DefiningPolynomials(minmapinv); curveMap:=DefiningPolynomials(Einv);

K:=pAdicField(13,Precision); EK:=BaseChange(E,K);
mod5:=13^5;
function AbsRes(a,m)
    a:=a mod m; return Min(a,(m-a) mod m);
end function;
function DeepKey(vals,p,m)
    best:=<m,m,m,m>;
    for i in [1..4] do
        x:=vals[i] mod m; if x mod p eq 0 then continue; end if;
        q:=InverseMod(x,m); zz:=Sort([AbsRes((z mod m)*q,m):z in vals]);
        key:=<zz[1],zz[2],zz[3],zz[4]>; if key lt best then best:=key; end if;
    end for;
    return best;
end function;
function ReadDeepBank(path,p,m)
    lines:=Split(Read(path),"\n"); dummy:=<Z!m,Z!m,Z!m,Z!m>;
    out:={dummy}; Exclude(~out,dummy);
    for i in [2..#lines] do if #lines[i] eq 0 then continue; end if;
        z:=Split(lines[i],"\t"); if #z lt 8 then continue; end if;
        vals:=[StringToInteger(z[j]):j in [5..8]];
        Include(~out,DeepKey(vals,p,m));
    end for;
    return out;
end function;
deepBank:=ReadDeepBank("results/target_22224_direct_contact_deep13_padic_tangent_p13.tsv",13,mod5);
function KP(P)
    if P eq Parent(P)!0 then return EK!0; end if;
    return EK![K!P[1],K!P[2],K!P[3]];
end function;
G1:=KP(free[1]); G2:=KP(free[2]); torsK:=[KP(P):P in tors];
PK<X,Y,Zc>:=PolynomialRing(K,3);
rawK:=[PK!h:h in rawMap]; curveK:=[PK!h:h in curveMap];

function EvalMap(polys,coords)
    vals:=[Evaluate(h,coords):h in polys];
    if &and[z eq 0:z in vals] then return false,vals; end if;
    return true,vals;
end function;
function PadicT(P)
    if P eq EK!0 then return true,false,K!1; end if;
    ok,raw:=EvalMap(rawK,[P[i]:i in [1..3]]); if not ok then return false,true,K!0; end if;
    ok,cp:=EvalMap(curveK,raw); if not ok or cp[3] eq 0 then return true,true,K!0; end if;
    try return true,false,cp[1]/cp[3];
    catch err
        // A projective denominator known only as O(13^e) is the pole chart
        // at the available precision; retain it conservatively as infinity.
        return true,true,K!0;
    end try;
end function;

lines:=Split(Read(input_file),"\n");
out:=Open(output_file,"w");
fprintf out,"m\tn\ttorsion_coset\tt_valuation\tt_mod13\tt_mod13_5\tchart\tdeep_p5\n";
tested:=0; finite:=0; infinity:=0; allowed:=0; unresolved:=0; deep:=0;
allowedFinite:={Z|0,1,5,8,12};
print "N2_P13_PADIC_FILTER_START","precision",Precision,"input",input_file,
      "E",E,"generators",gens,"projective_finite",Sort(Setseq(allowedFinite));
for i in [2..#lines] do if #lines[i] eq 0 then continue; end if;
    z:=Split(lines[i],"\t"); if #z lt 3 then continue; end if;
    m:=StringToInteger(z[1]); n:=StringToInteger(z[2]); ti:=StringToInteger(z[3]); tested+:=1;
    P:=m*G1+n*G2+torsK[ti]; ok,isinf,tt:=PadicT(P);
    if not ok then unresolved+:=1;
        fprintf out,"%o\t%o\t%o\t999\t-1\t-1\tunresolved\t0\n",m,n,ti;
        continue;
    end if;
    v:=Valuation(tt);
    if isinf or v lt 0 then
        infinity+:=1; allowed+:=1;
        fprintf out,"%o\t%o\t%o\t%o\t-1\t-1\tinfinity\t0\n",m,n,ti,v;
    else
        finite+:=1; tr:=(Z!tt) mod mod5; r:=tr mod 13;
        if r in allowedFinite then
            allowed+:=1;
            valsK:=[K!z:z in fixed] cat [K!d0*tt^2];
            valsZ:=[(Z!z) mod mod5:z in valsK];
            dp:=DeepKey(valsZ,13,mod5) in deepBank; if dp then deep+:=1; end if;
            fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\tfinite\t%o\n",m,n,ti,v,r,tr,dp select 1 else 0;
        end if;
    end if;
end for;
delete out;
print "N2_P13_PADIC_FILTER_DONE","tested",tested,"finite",finite,
      "infinity",infinity,"unresolved",unresolved,"allowed",allowed,
      "deep_p5",deep,"deep_bank",#deepBank,"output",output_file;
UnsetLogFile(); quit;
