//////////////////////////////////////////////////////////////////////
//  Emit a mod-343 (R,w) residue gate for the smooth p=7 strata of the
//  M_1(8,4) [4,16] second-halving cover.
//
//  This is a reusable version of the level-2 lifting code in
//  agent_m18_416_p7_blowup_level2.m.  It tests whether each (R,w) mod
//  343 over the selected mod-7 strata admits an aux lift
//  (a,b,c,d,e) mod 343 satisfying the exact cleared E416 system.
//
//  Output file format:
//      # comments...
//      R W R0 W0
//
//  Usage:
//    magma -b output_file:=data/agent_m18_416_smooth_strata_mod343_gate.txt \
//        code/agent_m18_416_smooth_strata_mod343_gate.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
Z := Integers();
p := 7;

if not assigned kercap then KerDimCap := 3;
elif Type(kercap) eq MonStgElt then KerDimCap := StringToInteger(kercap);
else KerDimCap := kercap; end if;

if not assigned output_file then
    output_file := "data/agent_m18_416_smooth_strata_mod343_gate.txt";
end if;

strata := [[3,3],[3,4],[4,0],[5,0],[5,2],[5,5]];
if assigned only then
    strata := [];
    for pr in Split(only, ";") do
        xy := Split(pr, ",");
        Append(~strata, [StringToInteger(xy[1]), StringToInteger(xy[2])]);
    end for;
end if;

// Optional partitioning over strata, useful if a future stratum is expensive.
if not assigned NParts then NParts := 1;
elif Type(NParts) eq MonStgElt then NParts := StringToInteger(NParts); end if;
if not assigned Part then Part := 0;
elif Type(Part) eq MonStgElt then Part := StringToInteger(Part); end if;

// ---- build cleared integer E416 ----
Rng<R,w,a,b,c,d,e> := PolynomialRing(Q, 7, "grevlex");
KF := FieldOfFractions(Rng);
PX<x> := PolynomialRing(KF);
t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
A := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
B := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x + (2*R^2 + R + 4*t);
f := x*A*B;
c4 := R + 2 + 4*t;
q416 := x^2 + a*x + b;
ell416 := c*x^2 + d*x + e;
F416 := f - ell416^2 - c4*(x+R)*q416^2;
E416 := [];
for i in [0..4] do
    ci := Numerator(Coefficient(F416, i));
    den := LCM([Denominator(cc) : cc in Coefficients(ci)]);
    ci := den*ci;
    g := GCD([Z!cc : cc in Coefficients(ci)]);
    Append(~E416, Rng!(ci/g));
end for;

auxidx := [3,4,5,6,7];
DerAux := [[Derivative(E416[i], auxidx[j]) : j in [1..5]] : i in [1..5]];
Fp := GF(p);

function Ev(pol, v)
    return Z!Evaluate(pol, v);
end function;

function AuxMod7(R0, w0)
    sols := [];
    for aa in [0..6] do for bb in [0..6] do for cc in [0..6] do
      for dd in [0..6] do for ee in [0..6] do
        v := [R0,w0,aa,bb,cc,dd,ee];
        ok := true;
        for i in [1..5] do
            if Ev(E416[i], v) mod p ne 0 then ok := false; break; end if;
        end for;
        if ok then Append(~sols, [aa,bb,cc,dd,ee]); end if;
    end for; end for; end for; end for; end for;
    return sols;
end function;

function LiftStep(Rw, aux, k)
    pk := p^k;
    base := [Rw[1], Rw[2]] cat aux;
    carry := [Fp!((Ev(E416[i], base)) div pk) : i in [1..5]];
    J := Matrix(Fp, 5, 5,
        [Fp!(Ev(DerAux[i][j], base)) : i in [1..5], j in [1..5]]);
    rhs := Vector(Fp, [-carry[i] : i in [1..5]]);
    consistent, part := IsConsistent(Transpose(J), rhs);
    if not consistent then return []; end if;
    ker := KernelMatrix(Transpose(J));
    dk := Nrows(ker);
    if dk gt KerDimCap then return "CAP"; end if;
    out := [];
    for cnt in [0..p^dk - 1] do
        coeffs := [];
        x0 := cnt;
        for _ in [1..dk] do
            Append(~coeffs, x0 mod p);
            x0 := x0 div p;
        end for;
        sig := part;
        for r in [1..dk] do
            sig := sig + (Fp!coeffs[r])*ker[r];
        end for;
        Append(~out, [aux[j] + pk*(Z!sig[j]) : j in [1..5]]);
    end for;
    return out;
end function;

function Survives343(Rw, aux7list)
    capped := false;
    for aux0 in aux7list do
        cur := [aux0];
        okchain := true;
        for k in [1,2] do
            nxt := [];
            for aux in cur do
                r := LiftStep(Rw, aux, k);
                if Type(r) eq MonStgElt then
                    capped := true;
                    okchain := false;
                    break;
                end if;
                nxt cat:= r;
            end for;
            if not okchain then break; end if;
            cur := nxt;
            if #cur eq 0 then break; end if;
        end for;
        if okchain and #cur gt 0 then return true; end if;
    end for;
    if capped then return "CAP"; end if;
    return false;
end function;

out := Open(output_file, "w");
fprintf out, "# M_1(8,4) [4,16] smooth-strata mod-343 gate\n";
fprintf out, "# kercap %o NParts %o Part %o\n", KerDimCap, NParts, Part;
fprintf out, "# columns: R W R0 W0\n";

printf "M18_416_SMOOTH_STRATA_MOD343_GATE kercap=%o output=%o\n",
    KerDimCap, output_file;
printf "strata=%o NParts=%o Part=%o\n", strata, NParts, Part;

total_survive := 0;
total_cap := 0;
total_killed := 0;

for si in [1..#strata] do
    if ((si-1) mod NParts) ne Part then continue; end if;
    rw0 := strata[si];
    R0 := rw0[1];
    w0 := rw0[2];
    aux7 := AuxMod7(R0, w0);
    nsurv := 0;
    ncap := 0;
    nkill := 0;
    for ri in [0..48] do for wi in [0..48] do
        Rv := R0 + p*ri;
        Wv := w0 + p*wi;
        s := Survives343([Rv,Wv], aux7);
        if Type(s) eq MonStgElt then
            ncap +:= 1;
        elif s then
            nsurv +:= 1;
            fprintf out, "%o %o %o %o\n", Rv, Wv, R0, w0;
        else
            nkill +:= 1;
        end if;
    end for; end for;
    total_survive +:= nsurv;
    total_cap +:= ncap;
    total_killed +:= nkill;
    printf "stratum <%o,%o>: aux7=%o survive=%o capped=%o killed=%o total=%o\n",
        R0, w0, #aux7, nsurv, ncap, nkill, nsurv+ncap+nkill;
end for;

printf "DONE survive=%o capped=%o killed=%o\n",
    total_survive, total_cap, total_killed;
print "Wrote", output_file;
quit;

