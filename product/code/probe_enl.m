// probe_enl.m — verify the u=2 (si5) enlargement end-to-end: original
// greedy cover gens vs Saturation(...,100) basis; regulator ratio; and the
// integral expression of the old gens in the new basis (the index).
SetColumns(0);
SetClassGroupBounds("GRH");
RQ := Rationals();
Qx<x> := PolynomialRing(RQ);
function AV(a, b, xv)
    if a lt b then
        if a eq 1 then
            return b eq 2 select (xv+3)*(xv-5) else 2*(xv-3);
        else
            return -(xv-1)*(xv-9);
        end if;
    else
        return -AV(b, a, xv);
    end if;
end function;
sg := [2,3,1];   // si5
u0 := RQ! 2;
d := 2*AV(sg[1], sg[3], u0); cA := AV(sg[1], sg[2], u0);
quart := cA*(x^2+6*d)*(x^2-2*d);
C := HyperellipticCurve(quart);
pts := Points(C : Bound := 1000);
aff := [ P : P in pts | P[3] ne 0 ];
E2, mE := EllipticCurve(C, aff[1]);
Em, phi := MinimalModel(E2);
printf "Em: %o\n", aInvariants(Em);
rlo, rhi := RankBounds(Em);
printf "bounds [%o,%o]\n", rlo, rhi;
S, mps := TwoDescent(Em);
cands := [];
for k in [1..#S] do
    cpts := Points(S[k] : Bound := 3000);
    for cp in cpts do
        EP := mps[k](cp);
        if Order(EP) ne 0 then continue; end if;
        if exists{ Q : Q in cands | Q eq EP or Q eq -EP } then continue; end if;
        Append(~cands, EP);
        if #cands ge 12 then break; end if;
    end for;
    if #cands ge 12 then break; end if;
end for;
gens := [];
for P in cands do
    if #gens ge 2 then break; end if;
    trial := gens cat [P];
    if Determinant(HeightPairingMatrix(trial)) gt 10.0^-6 then Append(~gens, P); end if;
end for;
printf "gens heights: %o\n", [ Height(P) : P in gens ];
Rsub := Determinant(HeightPairingMatrix(gens));
printf "Rsub = %o\n", Rsub;
gsat := Saturation(gens, 100);
printf "Saturation returned %o points, heights %o, orders %o\n",
    #gsat, [ Height(P) : P in gsat ], [ Order(P) : P in gsat ];
sel := [];
for P in gsat do
    if Order(P) ne 0 then continue; end if;
    if #sel ge 2 then break; end if;
    trial := sel cat [P];
    if Determinant(HeightPairingMatrix(trial)) gt 10.0^-8 then Append(~sel, P); end if;
end for;
Rafter := Determinant(HeightPairingMatrix(sel));
printf "Rafter = %o, ratio Rsub/Rafter = %o, index ~ %o\n", Rafter, Rsub/Rafter, Sqrt(Rsub/Rafter);
// express old gens in the new basis
HM := HeightPairingMatrix(sel);
HMi := HM^-1;
T0, mT := TorsionSubgroup(Em);
tors := [ mT(g) : g in T0 ];
for g in gens do
    vec := Matrix(RealField(30), 2, 1, [ HeightPairing(g, s) : s in sel ]);
    cf := HMi * vec;
    c1 := cf[1][1]; c2 := cf[2][1];
    Q := Round(c1)*sel[1] + Round(c2)*sel[2];
    exact := exists{ t : t in tors | Q + t eq g };
    printf "old gen = %o*new1 + %o*new2 (+tors): rounded (%o,%o) exact=%o\n", c1, c2, Round(c1), Round(c2), exact;
end for;
printf "PROBE_ENL_DONE\n";
quit;
