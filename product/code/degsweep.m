SetColumns(0);
SetMemoryLimit(6*10^9);
// run from product/code (house convention); set TORSION_JAC_ROOT to
// run from elsewhere
root := GetEnv("TORSION_JAC_ROOT");
if root ne "" then ChangeDirectory(root cat "/product/code"); end if;
load "split_lab.m";
D := CremonaDatabase();
// golden windows: (5,6), (5,24), (7,10): isogeny-induced psi = c*phi|E[N] is
// Galois-equivariant AND Kani-nondegenerate (p-adic descents) => the glued
// (N,N)-quotient is a Jacobian over Q with T1 x T2 injecting.
function WindowOK(m, n12)
    ok5 := GCD(5, n12) eq 1 and (m mod 25 ne 0) and (m in {6, 24});
    ok7 := GCD(7, n12) eq 1 and m eq 10;
    return ok5 or ok7, ok5 select (m eq 6 select "5-6" else "5-24") else "7-10";
end function;
nfound := 0;
for N in [11..300000] do
    nc := NumberOfIsogenyClasses(D, N);
    for i in [1..nc] do
        cls := EllipticCurves(D, N, i);
        if #cls lt 2 then continue; end if;
        Ts := [ Invariants(TorsionSubgroup(E)) : E in cls ];
        ords := [ IsEmpty(T) select 1 else &*T : T in Ts ];
        if Max(ords) lt 8 then continue; end if;
        for k1 in [1..#cls] do
            if ords[k1] lt 8 then continue; end if;
            for k2 in [k1+1..#cls] do
                n12 := ords[k1]*ords[k2];
                if n12 lt 96 then continue; end if;
                prod := Invariants(AbelianGroup(Ts[k1] cat Ts[k2]));
                if prod in KNOWN then continue; end if;
                okI, mp := IsIsogenous(cls[k1], cls[k2]);
                if not okI then continue; end if;
                m := Degree(mp);
                okW, tag := WindowOK(m, n12);
                printf "DEGPAIR cond=%o cls=%o deg=%o T=%ox%o prod=%o %o\n",
                    N, i, m, Ts[k1], Ts[k2], prod, okW select "WINDOW-" cat tag else "off-window";
                if okW then nfound +:= 1; end if;
            end for;
        end for;
    end for;
    if N mod 50000 eq 0 then printf "PROGRESSD %o (found %o)\n", N, nfound; end if;
end for;
printf "DEGSWEEP_DONE found=%o\n", nfound;
quit;
