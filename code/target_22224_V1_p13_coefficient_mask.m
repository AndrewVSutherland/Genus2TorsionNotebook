//////////////////////////////////////////////////////////////////////
// Convert the complete V1 p=13 depth-five T mask into matching-basis
// coefficient congruences.
//
// Put H1=4G1 and H2=4G2.  In the one-dimensional formal group,
//
//     H2 = c H1  mod depth relevant to T mod 13^5,
//     c = 10123 mod 13^4.
//
// Hence for m=r+4a, n=s+4b, T mod 13^5 depends only on
// q=a+c*b mod 13^4 and the base pair (r,s) mod 4.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(8);
Q:=Rationals(); Z:=Integers(); p:=13; depth:=5;
period:=p^4; cformal:=10123; K:=pAdicField(p,40);
R<T>:=PolynomialRing(Q);
fixed:=[Q!-18,Q!20,Q!75]; d0:=Q!-1470/121;
Rx:=(fixed[1]+d0*T^2)/(fixed[1]+d0);
Ry:=(fixed[2]+d0*T^2)/(fixed[2]+d0);
C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
Eraw,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi);
E,minmap:=MinimalModel(Eraw); minmapinv:=Inverse(minmap);
EK:=BaseChange(E,K);
G1:=EK![K!3603,K!216600,K!1];
G2:=EK![K!93,K!-2100,K!1]; H1:=4*G1; H2:=4*G2;

function FVal(P)
    if P eq EK!0 then return 1000000; end if;
    return Valuation(-P[1]/P[2]);
end function;
assert FVal(H2-cformal*H1) ge 5;

// Compose the two projective maps symbolically and cancel their common
// factor before p-adic evaluation.  Evaluating the two stages separately
// creates removable basepoints on several residue disks.
S<Xq,Yq,Zq>:=PolynomialRing(Q,3);
rawQ:=[S!h:h in DefiningPolynomials(minmapinv)];
curveQ:=[S!h:h in DefiningPolynomials(Einv)];
compQ:=[Evaluate(h,rawQ):h in curveQ];
common:=GCD(GCD(compQ[1],compQ[2]),compQ[3]);
redQ:=[ExactQuotient(h,common):h in compQ];
PK<X,Y,Zc>:=PolynomialRing(K,3);
tnum:=PK!redQ[1]; tden:=PK!redQ[3];
function PadicT(P)
    coords:=[P[i]:i in [1..3]];
    num:=Evaluate(tnum,coords); den:=Evaluate(tden,coords);
    if den eq 0 then return true,true,K!0; end if;
    try return true,false,num/den;
    catch err return false,true,K!0; end try;
end function;
function SquareCongruence(q,modulus)
    q:=Z!q mod modulus;
    if q eq 0 then return true; end if; // ramified boundary key: retain
    v:=0;
    while q mod p eq 0 do q div:=p; v+:=1; end while;
    return v mod 2 eq 0 and IsSquare(GF(p)!q);
end function;
function FullCoverKey(tkey)
    modulus:=p^depth; a:=Z!-2178; b:=Z!2420; cc:=Z!9075;
    d:=Z!-1470*tkey^2;
    rr:=[a*b*cc*d,
         a*(a+b)*(a+cc)*(a+d),
         b*(b+a)*(b+cc)*(b+d),
         cc*(cc+a)*(cc+b)*(cc+d)];
    return &and[SquareCongruence(z,modulus):z in rr];
end function;

output_file:="results/target_22224_V1_p13_coefficient_mask.tsv";
log_file:="results/target_22224_V1_p13_coefficient_mask.log";
SetLogFile(log_file:Overwrite:=true);
out:=Open(output_file,"w");
fprintf out,"m_mod4\tn_mod4\tq_residue\tq_modulus\tT_mod13pow5\tT_mod13\n";
total:=0; finite:=0; infinity:=0; unresolved:=0; allowed:=0;
counts:=AssociativeArray();
for r in [0..3] do for s in [0..3] do
    key:=Sprintf("%o,%o",r,s); localcount:=0; Pbase:=r*G1+s*G2;
    for q in [0..period-1] do
        // Fresh scalar multiplication avoids cumulative p-adic precision
        // loss from 28,560 successive additions.
        P:=Pbase+q*H1;
        total+:=1; ok,isinf,tt:=PadicT(P);
        if not ok then unresolved+:=1;
        elif isinf or Valuation(tt) lt 0 then infinity+:=1;
        else
            finite+:=1;
            tkey:=(Z!tt) mod p^depth;
            if FullCoverKey(tkey) then
                fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\n",
                        r,s,q,period,tkey,tkey mod p;
                allowed+:=1; localcount+:=1;
            end if;
        end if;
    end for;
    counts[key]:=localcount;
    print "COEFFICIENT_BASE",<r,s>,"allowed",localcount;
end for; end for;
delete out;
print "V1_P13_COEFFICIENT_MASK_DONE","period",period,"cformal",cformal,
      "total",total,"finite",finite,"infinity",infinity,
      "unresolved",unresolved,"allowed",allowed,"counts",counts,
      "output",output_file;
UnsetLogFile(); quit;
