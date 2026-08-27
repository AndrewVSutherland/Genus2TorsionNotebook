// claude_ov_b4_validate.m -- Lane 4 (route B4): validate the gp Galois-stable
// Richelot implementation (code/claude_ov_b4_richelot.gp) against Magma's
// RichelotIsogenousSurfaces, on (i) the recorded 2-rank-1 control and (ii) the
// smallest Flynn (6,2) member.  Also pins the seeds' exact torsion.
//
// Run: magma -b claude_ov_b4_validate.m

SetColumns(0);
if not assigned MemGB then MemGB := 8; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q);

function TwoRankOf(J) return #Invariants(TwoTorsionSubgroup(J)); end function;

// TorsionSubgroup needs an INTEGRAL model: f = c*g, c = a/b lowest terms,
// then (b*y)^2 = a*b*g, so a*b*g is an equivalent integral model.
function IntModel(f)
    cs := Coefficients(f);
    d := LCM([Denominator(c) : c in cs]);
    g := f*d;
    n := GCD([Numerator(c) : c in Coefficients(g)]);
    g := g/n;                       // integral primitive
    c := d/n;                       // f = (1/c)*g  =>  y^2=f  ~  (c*y)^2 = c*g
    return (Numerator(1/c)*Denominator(1/c))*g;
end function;

procedure Report(name, f)
    f := IntModel(f);
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    printf "%o: f = %o\n", name, f;
    printf "%o: factor degs = %o\n", name, Sort([Degree(t[1]) : t in Factorization(f)]);
    printf "%o: 2-rank = %o\n", name, TwoRankOf(J);
    printf "%o: G2Invariants = %o\n", name, G2Invariants(C);
    T, mp := TorsionSubgroup(J);
    printf "%o: TORSION = %o\n", name, Invariants(T);
end procedure;

// ---------------- control ------------------------------------------------
f0 := (x^2+x+1)*(x^4-4*x^2+2);
printf "=== CONTROL SEED ===\n";
C0 := HyperellipticCurve(f0); J0 := Jacobian(C0);
printf "control 2-rank = %o\n", TwoRankOf(J0);
printf "control G2Invariants = %o\n", G2Invariants(C0);
RS := RichelotIsogenousSurfaces(J0);
printf "control #RichelotIsogenousSurfaces = %o  types=%o\n", #RS, [Type(S) : S in RS];
for S in RS do
    if Type(S) eq JacHyp then
        CS := Curve(S);
        printf "control MAGMA codomain: f = %o\n", HyperellipticPolynomials(CS);
        printf "control MAGMA codomain G2Invariants = %o  2-rank = %o\n", G2Invariants(CS), TwoRankOf(S);
    end if;
end for;
gp0 := -2*x^5 - 24*x^4 - 64*x^3 - 32*x^2 - 4*x;   // as produced by code/claude_ov_b4_richelot.gp
Cg0 := HyperellipticCurve(gp0);
printf "control GP codomain G2Invariants = %o  2-rank = %o\n", G2Invariants(Cg0), TwoRankOf(Jacobian(Cg0));
for S in RS do if Type(S) eq JacHyp then
    printf "control GP == MAGMA codomain (isomorphic over Q)? %o ; twist by -1? %o\n",
        IsIsomorphic(Curve(S), Cg0), IsIsomorphic(Curve(S), HyperellipticCurve(-gp0));
end if; end for;

// ---------------- Flynn (6,2) smallest member ----------------------------
printf "=== FLYNN(6,2) SEED n=-1 (u=0, v=1/4, t=1/4) ===\n";
f1 := x^6 + 2*x^5 + 7/2*x^4 + 2*x^3 + 17/16*x^2 + 3/8*x + 1/16;
Report("seed1", f1);
f1i := IntModel(f1);
printf "seed1 integral model: %o\n", f1i;
C1 := HyperellipticCurve(f1i); J1 := Jacobian(C1);
RS1 := RichelotIsogenousSurfaces(J1);
printf "seed1 #RichelotIsogenousSurfaces = %o types=%o\n", #RS1, [Type(S) : S in RS1];
for S in RS1 do
    if Type(S) eq JacHyp then
        CS := Curve(S);
        printf "seed1 MAGMA codomain G2Invariants = %o  2-rank = %o\n", G2Invariants(CS), TwoRankOf(S);
        TT := TorsionSubgroup(S);
        printf "seed1 MAGMA codomain TORSION = %o\n", Invariants(TT);
    end if;
end for;
gp1 := 4368*x^6 + 5600*x^5 + 2744*x^4 - 112*x^3 - 343*x^2 - 42*x + 91;   // as produced by code/claude_ov_b4_richelot.gp
Cg1 := HyperellipticCurve(gp1); Jg1 := Jacobian(Cg1);
printf "seed1 GP codomain G2Invariants = %o  2-rank = %o\n", G2Invariants(Cg1), TwoRankOf(Jg1);
printf "seed1 GP codomain TORSION = %o\n", Invariants(TorsionSubgroup(Jg1));
for S in RS1 do if Type(S) eq JacHyp then
    printf "seed1 GP == MAGMA codomain (isomorphic over Q)? %o ; twist by -1? %o\n",
        IsIsomorphic(Curve(S), Cg1), IsIsomorphic(Curve(S), HyperellipticCurve(-gp1));
end if; end for;
pr := [p : p in PrimesInInterval(11,60)];
DS := Integers()!Discriminant(C1);
printf "seed1 #J(F_p) seed   = %o\n", [ (DS mod p ne 0) select #BaseChange(J1,GF(p)) else -1 : p in pr];
DG := Integers()!Discriminant(Cg1);
printf "seed1 #J(F_p) GP cod = %o\n", [ (DG mod p ne 0 and DS mod p ne 0) select #BaseChange(Jg1,GF(p)) else -1 : p in pr];

// the WRONG (twisted) codomain, to demonstrate the bug that was fixed
gp1bad := -9/4*gp1;
Cb := HyperellipticCurve(gp1bad); Jb := Jacobian(Cb);
printf "seed1 TWISTED-BY-A codomain TORSION = %o (this is the bug that was fixed)\n", Invariants(TorsionSubgroup(Jb));

printf "VALIDATE_DONE\n";
quit;
