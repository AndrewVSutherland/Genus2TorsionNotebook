// ===========================================================================
// Lane 3 : the REPAIRED RM screen, and the definitive validation of the
// derived Humbert-5 equation against LMFDB ground truth.
//
// WHAT WENT WRONG (found by code/claude_ov_l3_h5diag.m):
//   The lab's RM screen -- "the squarefree core of  n_p = a1^2 - 4(a2 - 2p)
//   is constant over good primes  <=>  RM; scattering <=> End = Z" -- is
//   WRONG for every geometrically-RM curve whose RM is NOT defined over Q
//   (LMFDB is_gl2_type = false).  At the primes p that are inert in the
//   quadratic field of definition of the RM, the extra twist forces
//   a1 = Tr(a_p) = 0 exactly; then n_p = 8p - 4a2 and its squarefree core is
//   junk that varies with p.  Those primes have density 1/2, so such a curve
//   LOOKS like End = Z under the naive screen.
//   In the 1000-row run of claude_ov_l3_h5validate.m, 613 of 1000 LMFDB
//   geom_end_alg='RM' curves scattered -- and 3116/4990 = 62.4% of all LMFDB
//   RM curves have is_gl2_type = false.
//
// THE FIX: only primes with a1 <> 0 (and n_p not a perfect square, and L_p
//   irreducible over Q) are informative.  Restricted to those, the core is
//   claimed to be constant = the RM discriminant core.
//
// This script tests the fix and the H5 equation simultaneously, on
//   4990 LMFDB geom_end_alg = 'RM' curves   (positive class)
//   2882 LMFDB geom_end_alg = 'Q'  curves   (End = Z control)
// and cross-tabulates
//   (naive census verdict, repaired census verdict, H5 = 0?, LMFDB truth).
//
// usage: code/claude_magma_slot.sh -b IN:=... NCH:=10 code/claude_ov_l3_h5validate2.m
// ===========================================================================
SetColumns(0);
SetMemoryLimit(4*10^9);
if not assigned IN then IN := "/tmp/claude-1000/-home-claude-torsion-jac/1c22f408-33bb-486e-870f-c5fca84bb522/scratchpad/l3_rmq.txt"; end if;
if not assigned NCH then NCH := 10; elif Type(NCH) eq MonStgElt then NCH := StringToInteger(NCH); end if;
if not assigned NMAX then NMAX := 100000; elif Type(NMAX) eq MonStgElt then NMAX := StringToInteger(NMAX); end if;
if not assigned PB then PB := 300; elif Type(PB) eq MonStgElt then PB := StringToInteger(PB); end if;
if not assigned OUT then OUT := "results/claude_ov_l3_h5validate2.log"; end if;

Z := Integers();  Q := Rationals();
Px<x> := PolynomialRing(Q);
R<J2,J4,J6,J10> := PolynomialRing(Q, 4);
H5 := R ! eval Read("results/claude_ov_l3_humbert5_eqn.txt");
printf "H5 loaded: #terms=%o\n", #Terms(H5);

// naive census (the lab screen) and repaired census, in one pass
function Census(g, B)
    dsc := Z ! Discriminant(g);
    dn := {};  dr := {};  np := 0; na1z := 0; nsq := 0;
    for p in PrimesInInterval(3, B) do
        if dsc mod p eq 0 then continue; end if;
        fp := PolynomialRing(GF(p)) ! g;
        if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
        L := LPolynomial(HyperellipticCurve(fp));
        a1 := Z!(-Coefficient(L,1));  a2 := Z!Coefficient(L,2);
        n := a1^2 - 4*(a2 - 2*p);
        np +:= 1;
        if n eq 0 then continue; end if;
        d := SquarefreeFactorization(n);
        if d eq 1 then nsq +:= 1; continue; end if;       // split over F_p
        Include(~dn, d);                                   // NAIVE screen
        if a1 eq 0 then na1z +:= 1; continue; end if;      // extra-twist prime
        Include(~dr, d);                                   // REPAIRED screen
    end for;
    return dn, dr, np, na1z, nsq;
end function;

lines := [l : l in Split(Read(IN), "\n") | #l gt 5];
if #lines gt NMAX then lines := lines[1..NMAX]; end if;
printf "CURVES %o  PB=%o\n", #lines, PB;

// PER-RUN scratch directory.  A fixed path here cost us a 52-minute run on
// 2026-07-25: a second invocation of this script rm -f'd the part files out
// from under the first one (recovered from /proc/PID/fd, see the header of
// results/claude_ov_l3_h5validate2.log).  Never share a scratch dir.
if not assigned DIR then DIR := Sprintf("/tmp/claude_ov_l3_h5val2_%o", Getpid()); end if;
dir := DIR;
System("mkdir -p " cat dir cat " && rm -f " cat dir cat "/part*.txt");

for ch in [0..NCH-1] do
    pid := Fork();
    if pid eq 0 then
        F := Open(dir cat "/part" cat IntegerToString(ch) cat ".txt", "w");
        for idx in [ch+1..#lines by NCH] do
            s := Split(lines[idx], ";");
            lab := s[1];
            fc := [StringToInteger(t) : t in Split(s[2], ",")];
            // NB Split drops empty fields, so realign from the END
            alg := s[#s-1];  gl2 := s[#s];
            hc := [Z|];
            if #s ge 5 and #s[3] gt 0 then
                hc := [StringToInteger(t) : t in Split(s[3], ",")];
            end if;
            f := Px ! fc;  hh := Px ! hc;
            g := 4*f + hh^2;
            if Degree(g) lt 5 or Discriminant(g) eq 0 then
                Puts(F, "ROW " cat lab cat " BAD"); continue;
            end if;
            C := HyperellipticCurve(g);
            I := IgusaClebschInvariants(C);
            v := Evaluate(H5, [Q!I[1], Q!I[2], Q!I[3], Q!I[4]]);
            dn, dr, np, na1z, nsq := Census(g, PB);
            Puts(F, Sprintf("ROW %o alg=%o gl2=%o H5zero=%o naive=%o rep=%o np=%o na1z=%o nsq=%o",
                 lab, alg, gl2, v eq 0, Sort(Setseq(dn)), Sort(Setseq(dr)), np, na1z, nsq));
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
            if #l gt 0 then Puts(G, l); n +:= 1; end if;
        end for;
    catch e ;
    end try;
end for;
Flush(G); delete G;
printf "VALIDATE2_DONE rows=%o\n", n;
quit;
