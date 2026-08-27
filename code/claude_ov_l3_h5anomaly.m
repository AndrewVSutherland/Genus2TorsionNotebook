// ===========================================================================
// Lane 3 : resolve the H5-vs-census disagreements on the LMFDB RM census.
//
// The Frobenius real-subfield-disc census tests the RM *field*; the Humbert
// equation H_5 tests the RM *order* (H_5 = disc-5 = the MAXIMAL order
// Z[(1+sqrt5)/2]; RM by the non-maximal order of conductor f lands on
// H_{5f^2}, a different Humbert surface).  So a curve with census {5} but
// H5 <> 0 should have a non-maximal RM order, which we detect from Frobenius:
// for RM by an order of discriminant Delta, D_p = a1^2-4(a2-2p) = Delta*m^2,
// so writing D_p = 5*k_p^2 the conductor f divides gcd_p(k_p).
//
// usage: magma -b LABELS:=... code/claude_ov_l3_h5anomaly.m
// ===========================================================================
SetColumns(0);
SetMemoryLimit(4*10^9);
if not assigned IN then IN := "/tmp/claude-1000/-home-claude-torsion-jac/1c22f408-33bb-486e-870f-c5fca84bb522/scratchpad/rmcurves.txt"; end if;
if not assigned ROWS then ROWS := "/tmp/claude-1000/-home-claude-torsion-jac/1c22f408-33bb-486e-870f-c5fca84bb522/scratchpad/anom.txt"; end if;
if not assigned OUT then OUT := "results/claude_ov_l3_h5anomaly.log"; end if;

Z := Integers();  Q := Rationals();
Px<x> := PolynomialRing(Q);
R<J2,J4,J6,J10> := PolynomialRing(Q, 4);
H5 := R ! eval Read("results/claude_ov_l3_humbert5_eqn.txt");

eqns := AssociativeArray();
for l in Split(Read(IN), "\n") do
    s := Split(l, ";");
    if #s lt 2 then continue; end if;
    eqns[s[1]] := <[StringToInteger(t) : t in Split(s[2], ",")],
                   (#s ge 3) select [StringToInteger(t) : t in Split(s[3], ",")] else []>;
end for;

G := Open(OUT, "w");
for lab in [t : t in Split(Read(ROWS), "\n") | #t gt 2] do
    if not IsDefined(eqns, lab) then continue; end if;
    e := eqns[lab];
    f := Px ! e[1]; hh := Px ! e[2];
    g := 4*f + hh^2;
    C := HyperellipticCurve(g);
    I := IgusaClebschInvariants(C);
    v := Evaluate(H5, [Q!I[1],Q!I[2],Q!I[3],Q!I[4]]);
    dsc := Discriminant(g);
    ks := [];  discs := {};
    for p in PrimesInInterval(3, 300) do
        if (Z!Numerator(dsc)) mod p eq 0 then continue; end if;
        fp := PolynomialRing(GF(p)) ! g;
        if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
        L := LPolynomial(HyperellipticCurve(fp));
        a1 := -Coefficient(L,1); a2 := Coefficient(L,2);
        D := (Z!a1)^2 - 4*((Z!a2) - 2*p);
        if D eq 0 then continue; end if;
        d0 := SquarefreeFactorization(D);
        Include(~discs, d0);
        if d0 eq 5 then Append(~ks, Isqrt(D div 5)); end if;
    end for;
    gk := (#ks gt 0) select GCD(ks) else 0;
    line := Sprintf("ANOM %o H5zero=%o discs=%o gcd_k=%o  => RM order disc divides %o",
        lab, v eq 0, Sort(Setseq(discs)), gk, 5*gk^2);
    printf "%o\n", line;  Puts(G, line);
end for;
Flush(G); delete G;
printf "ANOMALY_DONE\n";
quit;
