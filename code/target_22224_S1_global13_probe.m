//////////////////////////////////////////////////////////////////////
// Preliminary 13-adic structure probe for the S1 rank-one quotient.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(1813);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);
SetLogFile("results/target_22224_S1_global13_probe.log":Overwrite:=true);
fixed:=[Q!-4,Q!9,Q!30]; d0:=Q!-2166/245;
Rx:=(fixed[1]+d0*T^2)/(fixed[1]+d0);
Ry:=(fixed[2]+d0*T^2)/(fixed[2]+d0);
C:=HyperellipticCurve(Rx*Ry); P0:=C![Q!1,Q!1,Q!1];
Eraw,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi);
E,minmap:=MinimalModel(Eraw); minmapinv:=Inverse(minmap);
gens:=Generators(E); free:=[g:g in gens|Order(g) eq 0];
TG,tmap:=TorsionSubgroup(E); tors:=[tmap(g):g in TG];
print "S1_GLOBAL13_PROBE_START";
print "E",E;
print "discriminant",Discriminant(E),"v13",Valuation(Z!Discriminant(E),13);
print "gens",gens,"free",free,"tors",tors;
print "raw_map",DefiningPolynomials(minmapinv);
print "curve_map",DefiningPolynomials(Einv);
raw:=DefiningPolynomials(minmapinv);
// On the affine Weierstrass chart, the quartic parameter is
// t = raw[2]/(2*raw[1]+raw[2]).  Clear a common denominator from the
// two linear forms and inspect their unique base intersection.
PW<Xw,Yw,Zw>:=PolynomialRing(Q,3);
nq:=PW!raw[2]; dq:=PW!(2*raw[1]+raw[2]);
denND:=LCM([Denominator(Q!c):c in Coefficients(nq) cat Coefficients(dq)]);
nz:=denND*nq; dz:=denND*dq;
contND:=GCD([Abs(Z!c):c in Coefficients(nz) cat Coefficients(dz)]);
nz/:=contND; dz/:=contND;
print "t_num_integral",nz;
print "t_den_integral",dz;
// Coefficient order in a linear polynomial is constant,X,Y,Z internally;
// extract by evaluation to avoid relying on sparse Coefficients ordering.
function LinVec(h)
    return [Z!Evaluate(h,[Q!1,Q!0,Q!0])-Z!Evaluate(h,[Q!0,Q!0,Q!0]),
            Z!Evaluate(h,[Q!0,Q!1,Q!0])-Z!Evaluate(h,[Q!0,Q!0,Q!0]),
            Z!Evaluate(h,[Q!0,Q!0,Q!1])-Z!Evaluate(h,[Q!0,Q!0,Q!0])];
end function;
nv:=LinVec(nz); dv:=LinVec(dz);
B0:=[nv[2]*dv[3]-nv[3]*dv[2],nv[3]*dv[1]-nv[1]*dv[3],nv[1]*dv[2]-nv[2]*dv[1]];
gb:=GCD([Abs(z):z in B0]); B0:=[z div gb:z in B0];
print "nvec",nv,"dvec",dv,"base_intersection",B0;
a:=aInvariants(E);
Fhom:=Yw^2*Zw+a[1]*Xw*Yw*Zw+a[3]*Yw*Zw^2-Xw^3-a[2]*Xw^2*Zw-a[4]*Xw*Zw^2-a[5]*Zw^3;
fb:=Z!Evaluate(Fhom,[Q!z:z in B0]);
print "base_curve_value",fb,"valuation13",Valuation(fb,13);
Bpt:=E![Q!B0[1],Q!B0[2],Q!B0[3]];
print "base_point_on_E",Bpt;
for ti in [1..#tors] do for m in [-30..30] do
    if m*free[1]+tors[ti] eq Bpt then print "base_point_MW",m,ti; end if;
end for; end for;

print "cInvariants",cInvariants(E);
try print "local_information",LocalInformation(E,13); catch err print "local_information_error",err; end try;
F:=GF(13); RF<X>:=PolynomialRing(F);
print "reduced_cubic_factorization",Factorization(X^3+X^2+(F!3392)*X+F!62288);
for n in [1..20] do
    PP:=n*free[1];
    if PP eq E!0 then print "multiple",n,"zero";
    else print "multiple",n,"vx",Valuation(Q!PP[1],13),"vy",Valuation(Q!PP[2],13); end if;
end for;

K:=pAdicField(13,100); EK:=BaseChange(E,K);
GK:=EK![K!free[1][1],K!free[1][2],K!1];
o:=14;
H:=o*GK;
if H ne EK!0 then
    print "oG_vx",Valuation(H[1]),"oG_vy",Valuation(H[2]);
    // At infinity, z=-x/y is a standard formal parameter.
    print "oG_formal_z",-H[1]/H[2],"valuation",Valuation(-H[1]/H[2]);
end if;
for k in [1..7] do
    HH:=(o*13^(k-1))*GK;
    print "formal_depth",k,"mult",o*13^(k-1),
          "zval",HH eq EK!0 select 100 else Valuation(-HH[1]/HH[2]);
end for;
for k in [3..6] do
    mm:=3+14*13^k; WQ:=tors[2]; PP:=mm*GK+EK![K!WQ[1],K!WQ[2],K!1];
    nn:=K!-89425*PP[1]+K!1463*PP[2]+K!4228700*PP[3];
    dd:=K!-86279*PP[1]+K!1463*PP[2]-K!7215668*PP[3];
    print "near_base_depth",k+2,"m",mm,"vN",Valuation(nn),"vD",Valuation(dd);
end for;
print "S1_GLOBAL13_PROBE_DONE";
UnsetLogFile(); quit;
