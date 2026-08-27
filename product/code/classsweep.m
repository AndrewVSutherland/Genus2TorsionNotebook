SetColumns(0);
SetMemoryLimit(4*10^9);
D := CremonaDatabase();
// sweep isogeny classes, flag those whose torsion multiset contains a
// high-torsion pair with coprime-to-5-or-7 product not in KNOWN
// run from product/code (house convention); set TORSION_JAC_ROOT to
// run from elsewhere
root := GetEnv("TORSION_JAC_ROOT");
if root ne "" then ChangeDirectory(root cat "/product/code"); end if;
load "split_lab.m";
best := AssociativeArray();
for N in [20001..250000] do
    nc := NumberOfIsogenyClasses(D, N);
    for i in [1..nc] do
        cls := EllipticCurves(D, N, i);
        if #cls lt 2 then continue; end if;
        Ts := [ Invariants(TorsionSubgroup(E)) : E in cls ];
        for k1 in [1..#cls] do
            for k2 in [k1+1..#cls] do
                T1 := Ts[k1]; T2 := Ts[k2];
                n12 := (IsEmpty(T1) select 1 else &*T1) * (IsEmpty(T2) select 1 else &*T2);
                if n12 lt 96 then continue; end if;
                prod := Invariants(AbelianGroup(T1 cat T2));
                if prod in KNOWN then continue; end if;
                if GCD(35, n12) eq 35 then continue; end if;  // need N=5 or 7 coprime
                key := Sprintf("%o", prod);
                if IsDefined(best, key) then continue; end if;
                best[key] := <N, i, T1, T2>;
                printf "CLASSPAIR cond=%o class=%o T1=%o T2=%o prod=%o\n", N, i, T1, T2, prod;
            end for;
        end for;
    end for;
end for;
printf "CLASSSWEEP_DONE\n";
quit;
