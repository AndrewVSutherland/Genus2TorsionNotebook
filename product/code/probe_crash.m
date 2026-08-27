// probe_crash.m — isolate the si2 u=-49/2 segfault: run the fib2 pipeline
// stages with prints.
SetColumns(0);
SetMemoryLimit(5*10^9);
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
sg := [2,1,3];   // si2
u0 := RQ! -49/2;
d  := 2*AV(sg[1], sg[3], u0);
cA := AV(sg[1], sg[2], u0);
quart := cA*(x^2+6*d)*(x^2-2*d);
printf "quart = %o\n", quart;
C := HyperellipticCurve(quart);
// base point from the join partner t = 19/51 (as in the lane)
tv := RQ! 19/51;
oks, sv := IsSquare((tv-3)*d);
error if not oks, "no sv";
oky, Yv := IsSquare(Evaluate(quart, sv));
error if not oky, "no Yv";
P0 := C![sv, Yv];
printf "base from join: s=%o\n", sv;
E, mE := EllipticCurve(C, P0);
Em, phi := MinimalModel(E);
printf "Em = %o\n", aInvariants(Em);
rlo, rhi := RankBounds(Em);
printf "RankBounds [%o,%o]\n", rlo, rhi;
T, mT := TorsionSubgroup(Em);
printf "torsion %o\n", Invariants(T);
S, mps := TwoDescent(Em);
printf "covers: %o\n", #S;
cands := [];
for k in [1..#S] do
    cpts := Points(S[k] : Bound := 3000);
    printf "cover %o: %o pts\n", k, #cpts;
    for cp in cpts do
        EP := mps[k](cp);
        if Order(EP) ne 0 then continue; end if;
        if exists{ Q : Q in cands | Q eq EP or Q eq -EP } then continue; end if;
        Append(~cands, EP);
        if #cands ge 12 then break; end if;
    end for;
    if #cands ge 12 then break; end if;
end for;
printf "cands: %o\n", #cands;
gens := [];
for P in cands do
    if #gens ge 2 then break; end if;
    trial := gens cat [P];
    HM := HeightPairingMatrix(trial);
    if Determinant(HM) gt 10.0^-6 then Append(~gens, P); end if;
end for;
printf "indep gens: %o, heights %o\n", #gens, [ Height(P) : P in gens ];
printf "calling Saturation...\n";
gsat := Saturation(gens, 100);
printf "Saturation OK: %o pts\n", #gsat;
printf "PROBE_CRASH_DONE\n";
quit;
