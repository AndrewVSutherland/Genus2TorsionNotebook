// probe_stuck.m — stage-by-stage timing of the stuck si3 fiber u=35/8.
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
sg := [3,2,1];   // si3
u0 := RQ! 35/8;
d := 2*AV(sg[1], sg[3], u0); cA := AV(sg[1], sg[2], u0);
quart := cA*(x^2+6*d)*(x^2-2*d);
C := HyperellipticCurve(quart);
tv := RQ! 13/8;   // join partner
oks, sv := IsSquare((tv-3)*d);
oky, Yv := IsSquare(Evaluate(quart, sv));
error if not (oks and oky), "base build failed";
P0 := C![sv, Yv];
printf "base from join partner t=13/8\n";
E2, mE := EllipticCurve(C, P0);
Em, phi := MinimalModel(E2);
t0 := Realtime();
rlo, rhi := RankBounds(Em);
printf "RankBounds [%o,%o]: %o s\n", rlo, rhi, Realtime()-t0;
t0 := Realtime();
S, mps := TwoDescent(Em);
printf "TwoDescent: %o covers, %o s\n", #S, Realtime()-t0;
t0 := Realtime();
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
printf "cover points: %o cands, %o s\n", #cands, Realtime()-t0;
gens := [];
for P in cands do
    if #gens ge Min(rlo eq rhi select rlo else rhi, 4) then break; end if;
    trial := gens cat [P];
    if Determinant(HeightPairingMatrix(trial)) gt 10.0^-6 then Append(~gens, P); end if;
end for;
printf "gens: %o, heights %o\n", #gens, [ Height(P) : P in gens ];
t0 := Realtime();
gs := Saturation(gens, 20);
printf "Saturation(20): %o s\n", Realtime()-t0;
t0 := Realtime();
gs := Saturation(gens, 100);
printf "Saturation(100): %o s\n", Realtime()-t0;
t0 := Realtime();
okC := true; lb0 := 0.0; ub0 := 0.0;
try lb0, ub0 := CPSHeightBounds(Em); catch e okC := false; end try;
printf "CPSHeightBounds: ok=%o lb=%o ub=%o, %o s\n", okC, lb0, ub0, Realtime()-t0;
t0 := Realtime();
sb := SilvermanBound(Em);
printf "SilvermanBound: %o, %o s\n", sb, Realtime()-t0;
ub := Min(Max(ub0, 0.0), Max(sb, 0.0));
Bs := Ceiling(Exp(0.05 + ub));
printf "search bound would be %o\n", Bs;
t0 := Realtime();
sp := Points(Em : Bound := Min(Bs, 2000000));
printf "Points(%o): %o pts, %o s\n", Min(Bs, 2000000), #sp, Realtime()-t0;
printf "PROBE_STUCK_DONE\n";
quit;
