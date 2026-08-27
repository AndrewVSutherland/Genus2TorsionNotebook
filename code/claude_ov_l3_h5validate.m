// ===========================================================================
// Lane 3 : independent validation of the derived Humbert-5 equation H5
// against the LMFDB census of genus-2 curves with geometric RM.
//
// For each LMFDB RM curve  y^2 + h(x) y = f(x)  we compute
//   (i)  H5(I2,I4,I6,I10)   [the equation derived in claude_ov_l3_humbert5.m]
//   (ii) the Frobenius real-subfield-disc census over p <= 200
// and cross-tabulate.  The equation is correct iff
//   H5 = 0  <=>  census = {5}.
//
// usage: magma -b IN:=... NCH:=8 NMAX:=800 code/claude_ov_l3_h5validate.m
// ===========================================================================
SetColumns(0);
SetMemoryLimit(4*10^9);
if not assigned IN then IN := "/tmp/claude-1000/-home-claude-torsion-jac/1c22f408-33bb-486e-870f-c5fca84bb522/scratchpad/rmcurves.txt"; end if;
if not assigned NCH then NCH := 8; elif Type(NCH) eq MonStgElt then NCH := StringToInteger(NCH); end if;
if not assigned NMAX then NMAX := 800; elif Type(NMAX) eq MonStgElt then NMAX := StringToInteger(NMAX); end if;
if not assigned OUT then OUT := "results/claude_ov_l3_h5validate.log"; end if;

Z := Integers();  Q := Rationals();
Px<x> := PolynomialRing(Q);
R<J2,J4,J6,J10> := PolynomialRing(Q, 4);

H5 := eval Read("results/claude_ov_l3_humbert5_eqn.txt");
H5 := R ! H5;
printf "H5 loaded: #terms=%o\n", #Terms(H5);

function DiscCensus(fI, B)
    dsc := Discriminant(fI);
    discs := {};  np := 0;
    for p in PrimesInInterval(3, B) do
        if (Z!Numerator(dsc)) mod p eq 0 then continue; end if;
        Fp := GF(p);
        fp := PolynomialRing(Fp) ! fI;
        if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
        L := LPolynomial(HyperellipticCurve(fp));
        a1 := -Coefficient(L,1); a2 := Coefficient(L,2);
        n := (Z!a1)^2 - 4*((Z!a2) - 2*p);
        if n ne 0 then
            d := SquarefreeFactorization(n);
            if d ne 1 then Include(~discs, d); end if;
        end if;
        np +:= 1;
    end for;
    return discs, np;
end function;

lines := [l : l in Split(Read(IN), "\n") | #l gt 3];
if #lines gt NMAX then lines := lines[1..NMAX]; end if;
printf "CURVES %o\n", #lines;

dir := "/tmp/claude_ov_l3_h5val";
System("mkdir -p " cat dir cat " && rm -f " cat dir cat "/part*.txt");

for ch in [0..NCH-1] do
    pid := Fork();
    if pid eq 0 then
        F := Open(dir cat "/part" cat IntegerToString(ch) cat ".txt", "w");
        for idx in [ch+1..#lines by NCH] do
            s := Split(lines[idx], ";");
            lab := s[1];
            fc := [StringToInteger(t) : t in Split(s[2], ",")];
            hc := (#s ge 3) select [StringToInteger(t) : t in Split(s[3], ",")] else [];
            f := Px ! fc;  hh := Px ! hc;
            g := 4*f + hh^2;                  // y^2 = 4f + h^2
            if Degree(g) lt 5 or Discriminant(g) eq 0 then
                Puts(F, "ROW " cat lab cat " BAD"); continue;
            end if;
            C := HyperellipticCurve(g);
            I := IgusaClebschInvariants(C);
            v := Evaluate(H5, [Q!I[1], Q!I[2], Q!I[3], Q!I[4]]);
            dc, np := DiscCensus(g, 200);
            Puts(F, Sprintf("ROW %o H5zero=%o ndisc=%o discs=%o nprimes=%o",
                 lab, v eq 0, #dc, Sort(Setseq(dc)), np));
        end for;
        Flush(F); delete F;
        quit;
    end if;
end for;
WaitForAllChildren();

G := Open(OUT, "w");
n := 0; tab := AssociativeArray();
for ch in [0..NCH-1] do
    fn := dir cat "/part" cat IntegerToString(ch) cat ".txt";
    try
        for l in Split(Read(fn), "\n") do
            if #l gt 0 then Puts(G, l); n +:= 1; end if;
        end for;
    catch e ;
    end try;
end for;
Flush(G); delete G;
printf "VALIDATE_DONE rows=%o\n", n;
quit;
