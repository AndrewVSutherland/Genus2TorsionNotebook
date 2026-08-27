// Probe the one-dimensional p-adic closure of the two V1 MW generators.
SetColumns(0); SetSeed(8);
Q:=Rationals(); p:=13; K:=pAdicField(p,80);
E:=EllipticCurve([Q!0,Q!1,Q!0,Q!35967,Q!252063]);
EK:=BaseChange(E,K);
G1:=EK![K!3603,K!216600,K!1];
G2:=EK![K!93,K!-2100,K!1];
H1:=4*G1; H2:=4*G2;
function FVal(P)
    if P eq EK!0 then return 1000000; end if;
    return Valuation(-P[1]/P[2]);
end function;
print "V1_P13_PERIOD_PROBE","vH1",FVal(H1),"vH2",FVal(H2),
      "H1",H1,"H2",H2;

// Lift c so that H2-c*H1 lies successively deeper in the formal group.
c:=0; modulus:=1;
for level in [1..12] do
    bestc:=c; bestv:=-1000000;
    for digit in [0..p-1] do
        cc:=c+modulus*digit; vv:=FVal(H2-cc*H1);
        if vv gt bestv then bestv:=vv; bestc:=cc; end if;
    end for;
    c:=bestc; modulus*:=p;
    print "RELATION_LEVEL",level,"c",c,"modulus",modulus,
          "difference_formal_v",bestv;
end for;

// Match the periodic sieve's inverse map and measure its formal depth.
R<T>:=PolynomialRing(Q); fixed:=[Q!-18,Q!20,Q!75]; d0:=Q!-1470/121;
Rx:=(fixed[1]+d0*T^2)/(fixed[1]+d0);
Ry:=(fixed[2]+d0*T^2)/(fixed[2]+d0);
C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
Eraw,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi);
Emin,minmap:=MinimalModel(Eraw); minmapinv:=Inverse(minmap);
PK<X,Y,Zc>:=PolynomialRing(K,3);
rawK:=[PK!h:h in DefiningPolynomials(minmapinv)];
curveK:=[PK!h:h in DefiningPolynomials(Einv)];
function PT(P)
    if P eq EK!0 then return K!1; end if;
    raw:=[Evaluate(h,[P[i]:i in [1..3]]):h in rawK];
    cp:=[Evaluate(h,raw):h in curveK];
    return cp[1]/cp[3];
end function;
for q in [0..5] do
    tq:=PT(q*H1); print "T_SAMPLE",q,"v(T-1)",Valuation(tq-1),"T",tq;
end for;
for k in [1..7] do
    step:=p^k; t0:=PT(H1); t1q:=PT(H1+step*H1);
    print "T_PERIOD_STEP",k,"coefficient_step",step,
          "v_delta_T",Valuation(t1q-t0);
end for;
quit;
