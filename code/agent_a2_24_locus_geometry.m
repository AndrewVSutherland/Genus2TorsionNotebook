
//////////////////////////////////////////////////////////////////////
//  [2,24] locus geometry probe (go/no-go for the contact-3 derivation).
//
//  The [2,24] locus on the A(8) chart (r,p,t) is:
//     order-8  (FREE: the chart gives a rational order-8 D8 generically)
//   & 2-rank 2 (W-split, codim 1)
//   & rational 3-torsion (codim 1).
//  => a 1-dimensional variety in the 3-dim (r,p,t) space.
//
//  Count its F_ell-points N(ell) for several primes.  A 1-dim variety
//  has N(ell) ~ c*ell.  The genus shows in the fluctuation:
//     genus 0 (rational)  => N(ell) = ell+1 - (chart multiplicity terms),
//                            clean linear growth, ratio N/ell -> const;
//     higher genus        => O(sqrt(ell)) swings in N(ell) - c*ell.
//  If genus 0 => a rational parametrization exists => Newton-lift will
//  reconstruct low-height rational [2,24] curves.  If high genus =>
//  [2,24] is genuinely hard on this chart.
//
//  Over F_ell, "rational 3-torsion" <=> 3 | #J(F_ell) (finite field).
//  We report, per ell:
//     N_o8      : order-8 present (should be ~ell^3: order-8 is free)
//     N_o8_r2   : + 2-rank 2   (should be ~ell^2: codim 1)
//     N_full    : + 3 | #J     (should be ~ell  : codim 1) == [2,24] locus
//  and N_full/ell (-> constant if 1-dim) plus the raw N_full sequence.
//
//  Usage: magma -b Primes:="7,11,13,17,19,23,29,31" agent_a2_24_locus_geometry.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Z := Integers();

if not assigned Primes then Primes := "7,11,13,17,19,23,29,31";
elif Type(Primes) ne MonStgElt then Primes := Sprint(Primes); end if;
plist := [StringToInteger(s) : s in Split(Primes, ",")];

// A(8) chart over a field k (m=1 gauge), returns f, g8, ellBase over k
function A8f(x, rv, pv, tv)
    e := tv^2 - 2*pv*tv/rv; d := e + 2*pv - rv^2; lambda := rv/tv;
    u := pv + rv*tv - 2*rv;
    v := e + rv^2 - rv*pv - rv^2*tv + 3*pv*tv - rv*tv^2;
    a := rv^2 - lambda;
    b := 2*rv*pv - 2*lambda*(pv + rv*tv) + 2*rv*lambda;
    c := pv^2 + 2*pv*rv^2 - rv^4 - rv^3*tv - rv*pv^2/tv
         - lambda*(rv^2 + e) + 2*lambda*(rv*pv + rv^2*tv - 3*pv*tv + rv*tv^2);
    Qpoly := x^2 + d; q := a*x^2 + b*x + c; g8 := x^2 + u*x + v;
    f := q*(Qpoly^2 + q);
    L := rv*x + (pv - rv^2);
    ellBase := -(q + Qpoly*L);
    return f, g8, ellBase, L;
end function;

function TwoRankFF(fp)
    fac := Factorization(fp);
    degs := [Degree(g[1]) : g in fac];
    // repeated factors (non-squarefree) handled by caller; assume squarefree
    k := #degs; even := 0;
    for mask in [0..2^k-1] do
        ss := 0; for i in [1..k] do if (mask div 2^(i-1)) mod 2 eq 1 then ss +:= degs[i]; end if; end for;
        if ss mod 2 eq 0 then even +:= 1; end if;
    end for;
    return Ilog2(even) - 1;
end function;

printf "LOCUS GEOMETRY PROBE  primes=%o\n", plist;
printf "%-5o %-10o %-10o %-10o %-12o %-12o\n",
    "ell", "N_o8", "N_o8_r2", "N_full", "Nfull/ell", "Nfull-ell";
results := [];
for ell in plist do
    Fl := GF(ell);
    Pl<x> := PolynomialRing(Fl);
    N_o8 := 0; N_o8_r2 := 0; N_full := 0;
    for rr in Fl do
        if rr eq 0 then continue; end if;
        for tt in Fl do
            if tt eq 0 or rr*tt eq 1 then continue; end if;
            for pp in Fl do
                f, g8, ellBase, L := A8f(x, rr, pp, tt);
                if Degree(f) lt 5 then continue; end if;
                if not IsSquarefree(f) then continue; end if;
                // skip even sextic (bielliptic/split): x^1,x^3,x^5 all zero
                if Coefficient(f,1) eq 0 and Coefficient(f,3) eq 0 and Coefficient(f,5) eq 0 then continue; end if;
                // order-8 check: D8 = [g8, v8], want order exactly 8
                J := Jacobian(HyperellipticCurve(f));
                O := J!0;
                // This finite-field probe uses the unscaled A(8) model.
                // The integer-model denominator used in rational searches
                // should not be multiplied into the finite-field Mumford
                // representative.
                v8 := (-ellBase) mod g8;
                ok8 := true;
                try
                    D8 := elt<J | g8, v8>;
                    if 8*D8 ne O or 4*D8 eq O then ok8 := false; end if;
                catch ee ok8 := false; end try;
                if not ok8 then continue; end if;
                N_o8 +:= 1;
                if TwoRankFF(f) lt 2 then continue; end if;
                N_o8_r2 +:= 1;
                if #J mod 3 ne 0 then continue; end if;
                N_full +:= 1;
            end for;
        end for;
    end for;
    ratio := ell gt 0 select RealField(6)!N_full/ell else 0;
    printf "%-5o %-10o %-10o %-10o %-12o %-12o\n",
        ell, N_o8, N_o8_r2, N_full, ratio, N_full - ell;
    Append(~results, <ell, N_o8, N_o8_r2, N_full>);
end for;
print "RESULTS", results;
print "DONE";
quit;
