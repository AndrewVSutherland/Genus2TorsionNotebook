// claude_g2rec_run.m -- reconstruct a genus-2 curve /Q for a dim-2 weight-2
// newform given hecke cutters (Costa et al. ModularAbelianSurfaces pipeline).
// HC format: "p1:c0,c1,c2;p2:c0,c1,c2;..." (minpoly coefficients of a_p, ascending).
// Usage: nice -n 15 magma -b Lab:=16043.2.a.a N:=16043 HC:=... Prec:=300 claude_g2rec_run.m
// Dependency clones: pass ChimpDir:=... MASDir:=... on the command line, else
// default to ~/projects/ (the aws-spot-11 layout this pipeline ran on).
if not assigned ChimpDir then ChimpDir := GetEnv("HOME") cat "/projects/CHIMP"; end if;
if not assigned MASDir then MASDir := GetEnv("HOME") cat "/projects/ModularAbelianSurfaces"; end if;
AttachSpec(ChimpDir cat "/CHIMP.spec");
AttachSpec(MASDir cat "/spec");
SetColumns(0);
SetNthreads(1);
SetMemoryLimit(150*10^9);
SetVerbose("ModAbVarRec", 2);
SetVerbose("CurveRec", 1);
SetVerbose("HomologyModularSymbols", 1);
N := StringToInteger(N);
if not assigned Prec then Prec := 300; elif Type(Prec) eq MonStgElt then Prec := StringToInteger(Prec); end if;
hc := [];
for part in Split(HC, ";") do
    s := Split(part, ":");
    Append(~hc, <StringToInteger(s[1]), [StringToInteger(u) : u in Split(s[2], ",")]>);
end for;
printf "G2REC %o N=%o prec=%o cutters=%o\n", Lab, N, Prec, hc;
t0 := Cputime();
f := MakeNewformModSym(N, hc);
printf "MODSYM built: dim %o (%o s)\n", Dimension(f), Cputime(t0);
assert Dimension(f) eq 4;
for MaxEnd in [false, true] do
    printf "=== MaximalEnd = %o ===\n", MaxEnd;
    t1 := Cputime();
    found := false;
    try
        res := RationalGenus2Curves(f : Precision:=Prec, Quotient:=false, MaximalEnd:=MaxEnd, OnlyOne:=true);
        printf "#results = %o\n", #res;
        for i in [1..#res] do
            b, F, C, e := Explode(res[i]);
            printf "TRY %o found=%o\n", i, b;
            if b then
                found := true;
                printf "CURVE %o MaxEnd=%o : %o\n", Lab, MaxEnd, C;
                fC, hC := HyperellipticPolynomials(C);
                printf "CURVEEQ %o fcoeffs=%o hcoeffs=%o\n", Lab, Coefficients(fC), Coefficients(hC);
            else
                printf "NOCURVE try %o reason: %o ; data %o\n", i, e, C;
            end if;
        end for;
    catch er
        printf "ERROR MaxEnd=%o : %o\n", MaxEnd, er`Object;
    end try;
    printf "MaxEnd=%o done (%o s)\n", MaxEnd, Cputime(t1);
    if found then break; end if;
end for;
printf "G2RECDONE %o\n", Lab;
quit;
