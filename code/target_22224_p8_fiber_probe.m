//////////////////////////////////////////////////////////////////////
// Lightweight rank/maps/MW walk on the P8 shared-triple genus-5 fiber.
// Avoids the expensive full MordellWeilGroup saturation.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(8);
if not assigned MultipleBound then MultipleBound:=1000; end if;
if Type(MultipleBound) eq MonStgElt then MultipleBound:=StringToInteger(MultipleBound); end if;
if not assigned SearchBound then SearchBound:=10000; end if;
if Type(SearchBound) eq MonStgElt then SearchBound:=StringToInteger(SearchBound); end if;
if not assigned LatticeBound then LatticeBound:=8; end if;
if Type(LatticeBound) eq MonStgElt then LatticeBound:=StringToInteger(LatticeBound); end if;
if not assigned log_file then log_file:="results/target_22224_p8_fiber_probe.log"; end if;
if not assigned output_file then output_file:="results/target_22224_p8_fiber_probe.tsv"; end if;
SetLogFile(log_file : Overwrite:=true);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);
fixed:=[Q!528,-726,-891]; d0:=Q!50; t1:=Q!19/5;
out:=Open(output_file,"w");
fprintf out,"pair\tsource\tt\ta\tb\tc\td\tthree_survivor\tgcd_bound\n";

function ExactSquare(x)
    x:=Q!x; if x lt 0 then return false,Q!0; end if;
    a,ra:=IsSquare(Z!Numerator(x)); b,rb:=IsSquare(Z!Denominator(x));
    return a and b, (a and b) select Q!ra/rb else Q!0;
end function;

function FullAt(t)
    vals:=fixed cat [d0*t^2];
    a,b,c,d:=Explode(vals);
    rr:=[a*b*c*d,a*(a+b)*(a+c)*(a+d),
        b*(b+a)*(b+c)*(b+d),c*(c+a)*(c+b)*(c+d)];
    return &and[ExactSquare(x):x in rr],vals;
end function;

prime_list:=[11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97];
function ThreeBound(vals)
    f:=T*&*[T+z^2:z in vals]; g:=0;
    for p in prime_list do
        try
            fp:=ChangeRing(f,GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
            n:=Z!Evaluate(LPolynomial(HyperellipticCurve(fp)),1);
            g:=(g eq 0) select n else GCD(g,n);
            if g mod 3 ne 0 then return false,g; end if;
        catch e
            continue;
        end try;
    end for;
    return g ne 0 and g mod 3 eq 0,g;
end function;

print "P8_LIGHT_START","fixed",fixed,"d0",d0,"known_t",t1,
      "explicit_conics",[25*T^2+264,726-50*T^2,891-50*T^2];
seen:={Q|}; fullcount:=0; newcount:=0; threecount:=0;
for ij in [<1,2>,<1,3>,<2,3>] do
    i,j:=Explode(ij); x:=fixed[i]; y:=fixed[j];
    Rx:=(x+d0*T^2)/(x+d0); Ry:=(y+d0*T^2)/(y+d0);
    quartic:=Rx*Ry; C:=HyperellipticCurve(quartic);
    P0:=C![Q!1,Q!1,Q!1];
    okx,sx:=ExactSquare(Evaluate(Rx,t1)); oky,sy:=ExactSquare(Evaluate(Ry,t1));
    assert okx and oky; P1:=C![t1,sx*sy,Q!1];
    E,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi); Q1:=phi(P1);
    Emin:=MinimalModel(E); lo,hi:=RankBounds(Emin);
    TG,tmap:=TorsionSubgroup(E); tors:={tmap(g):g in TG};
    print "P8_LIGHT_QUOTIENT",ij,"quartic",quartic,"rank_bounds",lo,hi,
          "torsion",Invariants(TG),"minimal",Emin,"known_point",Q1,
          "known_order",Order(Q1),"inverse_map",DefiningPolynomials(Einv);
    gens:={Q1}; indep:=E!0; detmax:=RealField()!0;
    try
        pts:=PointSearch(E,SearchBound);
        for ep in pts do
            if Order(ep) ne 0 then continue; end if;
            Include(~gens,ep);
            try
                H:=HeightPairingMatrix([Q1,ep]); dd:=Abs(Determinant(H));
                if dd gt detmax then detmax:=dd;indep:=ep; end if;
            catch e
                dummy:=0;
            end try;
        end for;
        print "P8_POINTSEARCH",ij,"points",#pts,"non_torsion_generators",#gens,
              "best_pair_det",detmax,"second_point",indep;
    catch e
        print "P8_POINTSEARCH_ERROR",ij,e`Object;
    end try;
    epoints:={E!0};
    for gen in gens do for n in [-MultipleBound..MultipleBound] do Include(~epoints,n*gen); end for; end for;
    if indep ne E!0 then
        for m in [-LatticeBound..LatticeBound] do
            for n in [-LatticeBound..LatticeBound] do Include(~epoints,m*Q1+n*indep); end for;
        end for;
    end if;
    withtors:={ep+tor:ep in epoints,tor in tors};
    print "P8_WALK",ij,"epoints",#withtors;
    for ep in withtors do
        try cp:=Einv(ep); catch e continue; end try;
        if cp[3] eq 0 then continue; end if;
        tt:=Q!(cp[1]/cp[3]); if tt in seen then continue; end if;
        Include(~seen,tt); good,vals:=FullAt(tt); if not good then continue; end if;
        fullcount+:=1; known:=tt in {Q!1,Q!-1,t1,-t1}; if not known then newcount+:=1; end if;
        survives,g:=ThreeBound(vals); if survives then threecount+:=1; end if;
        fprintf out,"%o-%o\tMW\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
                i,j,tt,vals[1],vals[2],vals[3],vals[4],survives select 1 else 0,g;
        print "P8_FULL_POINT","pair",ij,"t",tt,"known",known,
              "three_survivor",survives,"gcd",g;
        if survives then
            f:=T*&*[T+z^2:z in vals]; G,mp:=TorsionSubgroup(Jacobian(HyperellipticCurve(f)));
            print "P8_EXACT_TORSION",tt,Invariants(G);
        end if;
    end for;
end for;
delete out;
print "P8_LIGHT_DONE","t_tested",#seen,"full",fullcount,"new",newcount,
      "three_survivors",threecount,"output",output_file;
UnsetLogFile(); quit;
