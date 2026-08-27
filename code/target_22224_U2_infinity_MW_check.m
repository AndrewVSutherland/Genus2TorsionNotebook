//////////////////////////////////////////////////////////////////////
// Locate the p=13 pole/infinity class among the 281 matching-basis U2
// Mordell--Weil survivors.  The complete local obstruction itself is in
// target_22224_U2_infinity_contact_lift.py; this file supplies the requested
// intersection with the imported coefficient list.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(181313);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);
p:=13; Precision:=100; K:=pAdicField(p,Precision);
input_file:="results/target_22224_U2_rank2_modular_N10000_fast.tsv";
output_file:="results/target_22224_U2_infinity_MW_check.tsv";
log_file:="results/target_22224_U2_infinity_MW_check.log";
SetLogFile(log_file:Overwrite:=true);

fixed:=[Q!-1071,Q!-1054,Q!1116]; d0:=Q!1134;
Rx:=(fixed[1]+d0*T^2)/(fixed[1]+d0);
Ry:=(fixed[3]+d0*T^2)/(fixed[3]+d0);
C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
Eraw,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi);
E,minmap:=MinimalModel(Eraw); minmapinv:=Inverse(minmap);
free:=[E![Q!1543/121,Q!8819100/1331,Q!1],
       E![Q!2008,Q!-105300,Q!1]];
tors:=[E!0,E![Q!-17,Q!0,Q!1]];
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
    if not ok then return true,true,K!0,"raw_basepoint"; end if;
    ok,cp:=EvalMap(curveK,raw);
    if not ok then return true,true,K!0,"curve_basepoint"; end if;
    if cp[3] eq 0 then return true,true,K!0,"projective_pole"; end if;
    try return true,false,cp[1]/cp[3],"finite";
    catch err return false,true,K!0,"division_error"; end try;
end function;

function CoverRadicands(base)
    a,b,c,d:=Explode(base);
    return [a*b*c*d,
            a*(a+b)*(a+c)*(a+d),
            b*(b+a)*(b+c)*(b+d),
            c*(c+a)*(c+b)*(c+d)];
end function;

lines:=Split(Read(input_file),"\n"); out:=Open(output_file,"w");
fprintf out,"m\tn\ttorsion_coset\tmap_status\tv13_T\tv13_z\tz_unit_mod13\tcover_square_pattern\n";
tested:=0; finite:=0; nonintegral:=0; exactpole:=0; unresolved:=0; fullcover:=0;
print "U2_INFINITY_MW_START","input",input_file,"precision",Precision,
      "generators",free,"torsion",tors;
for i in [2..#lines] do
    if #lines[i] eq 0 then continue; end if;
    w:=Split(lines[i],"\t"); if #w lt 3 then continue; end if;
    m:=StringToInteger(w[1]); n:=StringToInteger(w[2]); ti:=StringToInteger(w[3]);
    P:=m*G[1]+n*G[2]+torsK[ti]; tested+:=1;
    ok,isinf,tt,status:=PadicT(P);
    if not ok then
        PQ:=m*free[1]+n*free[2]+tors[ti]; exactstatus:="exact_map_error";
        try
            cpQ:=Einv(minmapinv(PQ));
            if cpQ[3] eq 0 then
                exactstatus:="exact_projective_pole"; exactpole+:=1;
                fprintf out,"%o\t%o\t%o\t%o\tinf\tinf\t-1\tdegenerate_z0\n",m,n,ti,exactstatus;
            else
                exactstatus:="exact_finite_after_padic_loss";
                tq:=Q!(cpQ[1]/cpQ[3]);
                print "U2_INFINITY_MW_EXACT_RECOVERY","t",tq;
                if Valuation(tq,13) ge 0 then
                    finite+:=1; tm:=(Numerator(tq) mod p)*InverseMod(Denominator(tq) mod p,p) mod p;
                    fprintf out,"%o\t%o\t%o\t%o\t%o\tna\t%o\tfinite_t_mod13_%o\n",
                            m,n,ti,exactstatus,Valuation(tq,13),tm,tm;
                else nonintegral+:=1; end if;
            end if;
        catch err
            unresolved+:=1;
            fprintf out,"%o\t%o\t%o\t%o\tunresolved\tunresolved\t-1\tunresolved\n",m,n,ti,exactstatus;
        end try;
        print "U2_INFINITY_MW_SPECIAL","m",m,"n",n,"ti",ti,
              "padic_status",status,"exact_status",exactstatus,"PQ",PQ;
        continue;
    end if;
    if isinf then
        exactpole+:=1;
        fprintf out,"%o\t%o\t%o\t%o\tinf\tinf\t-1\tdegenerate_z0\n",m,n,ti,status;
        print "U2_INFINITY_MW_CLASS","m",m,"n",n,"ti",ti,"status",status;
        continue;
    end if;
    vt:=Valuation(tt);
    if vt ge 0 then finite+:=1; continue; end if;
    nonintegral+:=1; z:=1/tt; vz:=Valuation(z);
    zu:=Z!(z/(K!p^vz)) mod p;
    base:=[K!fixed[j]*z^2:j in [1..3]] cat [K!d0];
    rads:=CoverRadicands(base); pattern:=[IsSquare(q):q in rads];
    if &and pattern then fullcover+:=1; end if;
    fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
            m,n,ti,status,vt,vz,zu,pattern;
    print "U2_INFINITY_MW_CLASS","m",m,"n",n,"ti",ti,
          "status",status,"vT",vt,"vz",vz,"zunit",zu,
          "cover_square_pattern",pattern;
end for;
delete out;
print "U2_INFINITY_MW_DONE","tested",tested,"finite",finite,
      "nonintegral",nonintegral,"exact_pole",exactpole,
      "unresolved",unresolved,"infinity_fullcover",fullcover,
      "output",output_file;
UnsetLogFile(); quit;
