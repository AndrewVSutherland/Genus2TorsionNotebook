// claude_g2rec_mdeep.m -- deep verification + eigenvalue/cutter extraction for
// one dim-2 weight-2 newform, via sign +1 modular symbols (deepverify.m pattern).
// Emits machine-readable lines:
//   APDATA p tp np      exact trace/norm of a_p for every good p <= 199
//   CUTTER p [c0,..]    minimal polynomial of a_p (hecke cutter material)
//   RESULT <label> M_odd=<m> over <k> primes
// Usage: nice -n 15 magma -b Lab:=16043.2.a.a N:=16043 TR:=2,-1,... claude_g2rec_mdeep.m
SetColumns(0);
SetMemoryLimit(80*10^9);
N := StringToInteger(N);
tr := [StringToInteger(s) : s in Split(TR, ",")];
printf "MDEEP %o N=%o start %o\n", Lab, N, Cputime();
M := ModularSymbols(N, 2, +1);
S := CuspidalSubspace(M);
NS := NewSubspace(S);
printf "new subspace dim %o (%o s)\n", Dimension(NS), Cputime();
D := NewformDecomposition(NS);
printf "decomposed: %o factors (%o s)\n", #D, Cputime();
target := 0;
for i in [1..#D] do
    if Dimension(D[i]) ne 2 then continue; end if;
    f := qEigenform(D[i], 41);
    K := BaseRing(Parent(f));
    if Degree(K) ne 2 then continue; end if;
    ok := true;
    for n in [2..40] do
        if Trace(K!Coefficient(f, n)) ne tr[n] then ok := false; break; end if;
    end for;
    if ok then target := i; break; end if;
end for;
if target eq 0 then printf "NOMATCH %o\n", Lab; quit; end if;
printf "matched factor %o (%o s)\n", target, Cputime();
A := D[target];
f := qEigenform(A, 200);
K := BaseRing(Parent(f));
printf "HECKEFIELD defpoly %o disc %o\n", DefiningPolynomial(K), Discriminant(MaximalOrder(K));
Z := Integers();
Mgcd := 0; used := 0; ncut := 0;
for p in PrimesInInterval(2, 199) do
    if N mod p eq 0 then continue; end if;
    ap := K!Coefficient(f, p);
    tp := Z!Trace(ap); np := Z!Norm(ap);
    printf "APDATA %o %o %o\n", p, tp, np;
    if ncut lt 6 then
        mp := MinimalPolynomial(ap);
        printf "CUTTER %o %o\n", p, [Z!c : c in Coefficients(mp)];
        ncut +:= 1;
    end if;
    if p gt 2 then
        loc := Z!Norm((1+p) - ap);
        Mgcd := Gcd(Mgcd, loc); used +:= 1;
    end if;
end for;
Modd := Mgcd; while Modd ne 0 and Modd mod 2 eq 0 do Modd div:= 2; end while;
printf "RESULT %o M_odd=%o over %o primes; fac=%o (%o s)\n", Lab, Modd, used, Factorization(Modd), Cputime();
printf "MDEEPDONE %o\n", Lab;
quit;
