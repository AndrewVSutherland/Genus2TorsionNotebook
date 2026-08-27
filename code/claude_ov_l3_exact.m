// ===========================================================================
// Lane 3 / route B2' : exact stage for sieve candidates.
//
// Input  : a file of lines "HIT w a b c d"  (or "w a b c d")
// Curve  : y^2 = f(x) = x (x-w) (x^2 + a x + b) (x^2 + c x + d)   (monic sextic)
// Output : one EXACT line per candidate with
//            exact CFOrder(D_inf), factor type, TorsionSubgroup invariants,
//            the real-subfield-disc RM census (constant set = RM, scatter = End=Z),
//            the two-prime strict End=Z certificate, and the absolute Igusa
//            invariants (isomorphism fingerprint, for dedup against the two
//            known RM witnesses).
//
// usage: magma -b CANDFILE:=... OUT:=... NCH:=8 MemGB:=4 code/claude_ov_l3_exact.m
// ===========================================================================
SetColumns(0);
if not assigned CANDFILE then CANDFILE := "results/claude_ov_l3_cands.txt"; end if;
if not assigned OUT then OUT := "results/claude_ov_l3_exact.log"; end if;
if not assigned NCH then NCH := 8; elif Type(NCH) eq MonStgElt then NCH := StringToInteger(NCH); end if;
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Z := Integers();  Q := Rationals();
P<x> := PolynomialRing(Q);
PT<T> := PolynomialRing(Q);

// ---- exact CF order of D_inf (pell-cf-order skill, verbatim guards) -------
function SqrtPolyPart(f)
    s := x^3;
    for k in [1..3] do
        d := f - s^2;
        if Degree(d) le 2 then break; end if;
        s := s + (Coefficient(d, 6-k)/(2*Coefficient(s,3)))*x^(3-k);
    end for;
    return s;
end function;

function CFOrder(f, maxsteps)
    s := SqrtPolyPart(f);
    Pi := P!0; Qi := P!1; total := 0;
    for i in [0..maxsteps] do
        if Qi eq 0 then return 0; end if;
        ai := (Pi + s) div Qi;
        total +:= Degree(ai);
        Pn := ai*Qi - Pi;
        if (f - Pn^2) mod Qi ne 0 then return 0; end if;
        Qn := (f - Pn^2) div Qi;
        Pi := Pn; Qi := Qn;
        if i ge 1 and Degree(Qi) le 0 and Qi ne 0 then return total; end if;
    end for;
    return 0;
end function;

// ---- self-test of the CF tool (skill-mandated) ---------------------------
f14 := (x^2+1)*(x^4+5*x^2+4*x+4);
f18 := (x^2-x+1)*(x^4-x^3+9*x^2+8*x-8);
assert CFOrder(f14, 40) eq 14;
assert CFOrder(f18, 40) eq 18;
printf "CF_SELFTEST_OK f14=14 f18=18\n";

// ---- Frobenius data at a good prime --------------------------------------
// returns ok, a1, a2, chi
function FrobData(f, p)
    if p le 2 then return false, 0, 0, PT!0; end if;
    dsc := Discriminant(f);
    if dsc eq 0 then return false, 0, 0, PT!0; end if;
    if (Z!Numerator(dsc)) mod p eq 0 then return false, 0, 0, PT!0; end if;
    if (Z!Denominator(dsc)) mod p eq 0 then return false, 0, 0, PT!0; end if;
    Fp := GF(p);
    fp := PolynomialRing(Fp) ! f;
    if Degree(fp) lt 5 or not IsSquarefree(fp) then return false, 0, 0, PT!0; end if;
    C := HyperellipticCurve(fp);
    L := LPolynomial(C);                       // 1 - a1 T + a2 T^2 - a1 p T^3 + p^2 T^4
    a1 := -Coefficient(L, 1);
    a2 :=  Coefficient(L, 2);
    chi := T^4 - a1*T^3 + a2*T^2 - a1*p*T + p^2;
    return true, Z!a1, Z!a2, chi;
end function;

// squarefree core of the real-quadratic-subfield discriminant
function RealSubDisc(a1, a2, p)
    n := a1^2 - 4*(a2 - 2*p);
    if n eq 0 then return 0; end if;
    d := SquarefreeFactorization(n);
    return d;
end function;

// root-power strictness at p: chi irreducible AND no degree drop for n in [2..12]
function PowerStrict(chi)
    if not IsIrreducible(chi) then return false; end if;
    K := NumberField(chi);
    pi := K.1;
    for nn in [2..12] do
        if Degree(MinimalPolynomial(pi^nn)) lt 4 then return false; end if;
    end for;
    return true;
end function;

function FactorTypeStr(f)
    fac := Factorization(f);                    // degree-only use: SAFE (see conventions)
    ds := Sort([Degree(t[1]) : t in fac | true]);
    dl := [];
    for t in fac do for j in [1..t[2]] do Append(~dl, Degree(t[1])); end for; end for;
    return Sprint(Sort(dl));
end function;

// ---- read candidates ------------------------------------------------------
lines := Split(Read(CANDFILE), "\n");
cands := [];
for l in lines do
    s := [t : t in Split(l, " ") | #t gt 0];
    if #s eq 6 and s[1] eq "HIT" then s := s[2..6]; end if;
    if #s ne 5 then continue; end if;
    ok := true;
    v := [];
    // accept both integers and rationals "n/d" (the B2' sweep emits rationals)
    for t in s do
        try
            u := Split(t, "/");
            if #u eq 1 then
                Append(~v, Q ! StringToInteger(u[1]));
            elif #u eq 2 then
                den := StringToInteger(u[2]);
                if den eq 0 then ok := false; else Append(~v, (Q ! StringToInteger(u[1]))/den); end if;
            else
                ok := false;
            end if;
        catch e ok := false; end try;
    end for;
    if ok and #v eq 5 then Append(~cands, v); end if;
end for;
printf "CANDIDATES %o\n", #cands;

dir := "/tmp/claude_ov_l3_fork";
System("mkdir -p " cat dir cat " && rm -f " cat dir cat "/part*.txt");

for ch in [0..NCH-1] do
    pid := Fork();
    if pid eq 0 then                                          // ---- child ----
        F := Open(dir cat "/part" cat IntegerToString(ch) cat ".txt", "w");
        for idx in [ch+1..#cands by NCH] do
            v := cands[idx];
            w := v[1]; aa := v[2]; bb := v[3]; cc := v[4]; dd := v[5];
            f := x*(x-w)*(x^2+aa*x+bb)*(x^2+cc*x+dd);
            tag := Sprintf("(%o,%o,%o,%o,%o)", w, aa, bb, cc, dd);
            if Discriminant(f) eq 0 then
                Puts(F, "EXACT " cat tag cat " SINGULAR"); continue;
            end if;
            ord := CFOrder(f, 40);
            ftype := &cat[Sprintf("%o.", t) : t in Split(FactorTypeStr(f), " []," )];
            // TorsionSubgroup demands an INTEGRAL model; the B2' sweep emits
            // rational (p1..p4), so scale x -> X/m, y -> Y/m^3 (deg f = 6).
            m := LCM([Denominator(c) : c in Coefficients(f)] cat [Integers()|1]);
            fI := P ! [ Integers() ! (Coefficient(f,i)*m^(6-i)) : i in [0..6] ];
            okc := true;
            try
                C := HyperellipticCurve(fI);
                J := Jacobian(C);
                Tg := TorsionSubgroup(J);
                inv := Invariants(Tg);
                ai := AbsoluteInvariants(C);
            catch e
                okc := false;
            end try;
            if not okc then
                Puts(F, "EXACT " cat tag cat " CURVE_ERROR"); continue;
            end if;
            // RM prescreen + certificate scan in one prime loop
            discs := {};  certs := [];
            for p in PrimesInInterval(3, 200) do
                okp, a1, a2, chi := FrobData(fI, p);
                if not okp then continue; end if;
                D := RealSubDisc(a1, a2, p);
                if D ne 1 and D ne 0 then Include(~discs, D); end if;
                if #certs lt 6 and a1 ne 0 and PowerStrict(chi) then
                    Append(~certs, <p, D, chi>);
                end if;
            end for;
            // two-prime strict certificate: two power-strict primes with a1<>0
            // whose real quadratic subfields DIFFER (=> no RM, splitting fields
            // not both inside an RM tower)
            certok := false; cp := 0; cq := 0;
            for i in [1..#certs] do
                for j in [i+1..#certs] do
                    if certs[i][2] ne certs[j][2] then
                        certok := true; cp := certs[i][1]; cq := certs[j][1]; break i;
                    end if;
                end for;
            end for;
            rm := (#discs le 1);
            fp := &cat[Sprintf("%o;", c) : c in ai];
            cs := &cat[Sprintf("%o@%o|", c[1], c[2]) : c in certs];
            Puts(F, Sprintf("EXACT %o CFord=%o type=%o TORSION=%o ndisc=%o discs=%o RMflag=%o cert=%o certprimes=(%o,%o) certlist=%o absinv=%o",
                 tag, ord, ftype, inv, #discs, Sort(Setseq(discs)), rm, certok, cp, cq, cs, fp));
        end for;
        Flush(F); delete F;
        quit;
    end if;
end for;
WaitForAllChildren();

G := Open(OUT, "w");
n := 0;
for ch in [0..NCH-1] do
    fn := dir cat "/part" cat IntegerToString(ch) cat ".txt";
    try
        for l in Split(Read(fn), "\n") do
            if #l gt 0 then printf "%o\n", l; Puts(G, l); n +:= 1; end if;
        end for;
    catch e ;
    end try;
end for;
Puts(G, Sprintf("EXACT_DONE processed=%o of %o", n, #cands));
Flush(G); delete G;
printf "EXACT_DONE processed=%o of %o\n", n, #cands;
quit;
