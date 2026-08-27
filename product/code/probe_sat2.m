// probe_sat2.m — Saturation cost scaling on a heavy fiber's generators
// (si2 u=-16): time Saturation(gens, B) for B = 100, 500, 1000, 5000.
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
sg := [2,1,3]; u0 := RQ! -16;
d := 2*AV(sg[1], sg[3], u0); cA := AV(sg[1], sg[2], u0);
quart := cA*(x^2+6*d)*(x^2-2*d);
C := HyperellipticCurve(quart);
pts := Points(C : Bound := 200);
aff := [ P : P in pts | P[3] ne 0 ];
E2, mE := EllipticCurve(C, aff[1]);
Em, phi := MinimalModel(E2);
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
printf "gens %o heights %o\n", #gens, [ Height(P) : P in gens ];
for B in [100, 500, 1000, 5000] do
    t0 := Realtime();
    gs := Saturation(gens, B);
    printf "Saturation(gens, %o): %o pts %o s\n", B, #gs, Realtime()-t0;
end for;
printf "PROBE_SAT2_DONE\n";
quit;
