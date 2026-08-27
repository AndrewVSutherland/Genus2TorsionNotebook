SetColumns(0);
SetMemoryLimit(12*10^9);
Q := Rationals();
Px<x> := PolynomialRing(Q);
P2<R,S> := PolynomialRing(Q, 2);
F13 := R^3 - R^2*S^4 + 5*R^2*S^3 - 9*R^2*S^2 + 4*R^2*S - 2*R^2 - R*S^3 + 6*R*S^2 - 3*R*S + R - S^3;
// s0 = 1 fiber: from the lane: K = Q(sqrt 193)
A2q<uu,vv> := AffineSpace(Q,2);
hcr := hom< P2 -> CoordinateRing(A2q) | [uu,vv] >;
CX := ProjectiveClosure(Curve(A2q, hcr(F13)));
okh, XH, mph := IsHyperelliptic(CX);
DPm := DefiningPolynomials(mph);
RA := CoordinateRing(Ambient(CX));
hz1 := hom< RA -> P2 | [P2.1, P2.2, 1] >;
XPa := hz1(DPm[1]); ZPa := hz1(DPm[3]);
K := QuadraticField(193);
w := K.1;
R2K := PolynomialRing(K, 2);
h2 := hom< P2 -> R2K | [R2K.1, R2K.2] >;
IK := ideal< R2K | h2(F13), h2(XPa) - 1*h2(ZPa) >;
V := Variety(IK);
printf "fiber points: %o\n", #V;
E := 0; b0 := 0; c0 := 0;
for pv in V do
    r0 := pv[1]; sr := pv[2];
    if r0 in Q and sr in Q then continue; end if;
    c0 := sr*(r0-1); b0 := r0*c0;
    if b0 eq 0 then continue; end if;
    E0 := EllipticCurve([1-c0, -b0, -b0, 0, 0]);
    if Order(E0![0,0]) eq 13 then E := E0; break; end if;
end for;
printf "E(K) tors: %o\n", Invariants(TorsionSubgroup(E));
sig := hom< K -> K | -w >;
Es := EllipticCurve([sig(a) : a in aInvariants(E)]);
printf "j equal: %o\n", jInvariant(E) eq jInvariant(Es);
istw, dtw := IsQuadraticTwist(E, Es);
printf "IsQuadraticTwist: %o d=%o\n", istw, istw select dtw else 0;
f1 := HyperellipticPolynomials(WeierstrassModel(E));
g1 := HyperellipticPolynomials(WeierstrassModel(Es));
L := SplittingField(f1);
printf "L degree %o\n", Degree(L);
RL := PolynomialRing(L);
al := [ rr[1] : rr in Roots(RL!f1) ];
gm := [ rr[1] : rr in Roots(RL!g1) ];
printf "roots %o %o\n", #al, #gm;
// which scalings map al to gm?
for i in [1..3] do
    sc := gm[i]/al[1];
    prm := [ Index(gm, L!sc*al[k]) : k in [1..3] ];
    printf "scaling by gm[%o]/al[1]: prm = %o; sc in K: %o\n", i, prm, sc in K;
end for;
quit;
// Automorphisms path
GA := Automorphisms(L);
GK := [ gA : gA in GA | gA(L!w) eq L!w ];
printf "autos: %o total, %o over K\n", #GA, #GK;
Gaction := func< rs0 | [[ Index(rs0, gA(rr)) : rr in rs0 ] : gA in GK ] >;
nequi := 0;
for sgm0 in SymmetricGroup(3) do
    beta0 := [ gm[i^sgm0] : i in [1..3] ];
    if Gaction(al) eq Gaction(beta0) then nequi +:= 1; printf "equivariant matching: %o\n", sgm0; end if;
end for;
printf "total equivariant: %o\n", nequi;
quit;
