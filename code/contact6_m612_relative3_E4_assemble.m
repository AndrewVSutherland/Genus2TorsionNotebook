//////////////////////////////////////////////////////////////////////
// Assemble the exact S3-quotient quartic mp4Q (and C3-octic mp8Q) over
// Q(e) from per-prime residue dumps data/contact6_m612_E4_modp_<p>.txt
// written by code/contact6_m612_relative3_E4_exact.m (OnlyP mode).
//
// CRT + rational reconstruction per coefficient; then:
//   - resolvent consistency  Res_w(mp8, Y - w^2) = +- mp4^2;
//   - exact certification over Q(e): both irreducible,
//     Genus(E4) = 2, Genus(E8) = 4  (pilot values, exact here).
// Dumps mp4Q/mp8Q for stage 2 (rank bounds).
//
// Usage: magma -b Plist:="7,11,13,17,19" code/contact6_m612_relative3_E4_assemble.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals();
K<e> := FunctionField(Q); Kz<zz> := PolynomialRing(K);

if not assigned Plist then Plist := "7,11,13,17"; end if;
plist := [StringToInteger(s) : s in Split(Plist, ",")];

dumps := [* *]; goodp := [];
for p in plist do
    fn := Sprintf("data/contact6_m612_E4_modp_%o.txt", p);
    ok := true;
    try
        rec := eval Read(fn);
        Append(~dumps, rec); Append(~goodp, p);
    catch ee
        printf "missing/failed dump for p=%o; skipping\n", p;
    end try;
end for;
printf "assembling from primes %o\n", goodp;
N := &*goodp; ZN := Integers(N);
printf "CRT modulus ~ 10^%o ; max reconstructible height ~ %o\n",
    Ilog(10, N), Isqrt(N div 2);

// rec[2] = mp4 table (4 coefficients), rec[3] = mp8 table (8 coefficients)
// table[i] = <numlifts, denlifts> with denominator monic normalization
function Assemble(slot, deg)
    co := [K | ];
    for i in [1..deg] do
        // consistent shapes
        nls := {#dumps[j][slot][i][1] : j in [1..#dumps]};
        dls := {#dumps[j][slot][i][2] : j in [1..#dumps]};
        if #nls ne 1 or #dls ne 1 then return false, Kz!0, i; end if;
        nl := Rep(nls); dl := Rep(dls);
        numco := []; denco := [];
        for h in [1..nl] do
            residues := [dumps[j][slot][i][1][h] : j in [1..#dumps]];
            if &and[r eq 0 : r in residues] then Append(~numco, Q!0); continue; end if;
            x := CRT(residues, goodp);
            okr, q := RationalReconstruction(ZN!x);
            if not okr then return false, Kz!0, i; end if;
            Append(~numco, q);
        end for;
        for h in [1..dl] do
            residues := [dumps[j][slot][i][2][h] : j in [1..#dumps]];
            x := CRT(residues, goodp);
            okr, q := RationalReconstruction(ZN!x);
            if not okr then return false, Kz!0, i; end if;
            Append(~denco, q);
        end for;
        numQ := &+[K | numco[h]*e^(h-1) : h in [1..#numco]];
        denQ := &+[K | denco[h]*e^(h-1) : h in [1..#denco]];
        if denQ eq 0 then return false, Kz!0, i; end if;
        Append(~co, numQ/denQ);
    end for;
    return true, Kz!(co cat [K!1]), 0;
end function;

ok4, mp4Q, bad4 := Assemble(2, 4);
if not ok4 then printf "mp4 reconstruction FAILED at coefficient %o (need more primes)\n", bad4-1; quit; end if;
print "mp4Q ="; print mp4Q;
ok8, mp8Q, bad8 := Assemble(3, 8);
if not ok8 then printf "mp8 reconstruction FAILED at coefficient %o (need more primes)\n", bad8-1; quit; end if;
print "mp8Q ="; print mp8Q;

// resolvent consistency
KY<YY> := PolynomialRing(K); KYw<ww> := PolynomialRing(KY);
mp8w := &+[(KY!Coefficient(mp8Q,i))*ww^i : i in [0..8]];
resid := Resultant(mp8w, YY - ww^2);
mp4Y := &+[(KY!Coefficient(mp4Q,i))*YY^i : i in [0..4]];
consistent := resid eq mp4Y^2 or resid eq -mp4Y^2;
printf "RESOLVENT CONSISTENCY: %o\n", consistent;
if not consistent then print "INCONSISTENT — need more primes"; quit; end if;

// exact certification
assert Degree(mp4Q) eq 4 and IsMonic(mp4Q) and IsIrreducible(mp4Q);
E4<a4> := ext<K|mp4Q>;
g4 := Genus(E4);
printf "EXACT E4 GENUS = %o\n", g4; assert g4 eq 2;
assert IsIrreducible(mp8Q);
E8<a8> := ext<K|mp8Q>;
g8 := Genus(E8);
printf "EXACT E8 GENUS = %o\n", g8; assert g8 eq 4;
PrintFile("data/contact6_m612_E4_mp4Q.txt", Sprintf("%m", mp4Q) : Overwrite:=true);
PrintFile("data/contact6_m612_E8_mp8Q.txt", Sprintf("%m", mp8Q) : Overwrite:=true);
print "E4_EXACT_CERTIFIED";
quit;
