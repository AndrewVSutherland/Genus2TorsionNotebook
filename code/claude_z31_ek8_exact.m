// claude_z31_ek8_exact.m  (Task B4, exact stage)
// Input: a CAND (r,s) from the C sieve (claude_z31_ek8_sieve.c sweep mode).
// Rebuilds the genus-2 curve /Q with the EK disc-8 Igusa-Clebsch invariants,
// recovers the quadratic twist d whose Jacobian can carry rational 31-torsion
// (via chi_p(+-1) mod 31 sign constraints at many good primes), and runs
// TorsionSubgroup on the surviving twists.
//
// Usage:
//   magma -b Rn:=6000 Rd:=96721 Sn:=3557520 Sd:=96721 claude_z31_ek8_exact.m
// Markers: EXACT_START / NO_CURVE_OVER_Q / DEAD_AT_p / TWIST_CAND d /
//          TORSION d [..] / HIT31 d / EXACT_DONE
SetColumns(0);
SetMemoryLimit(8*10^9);
Q := Rationals();

if not assigned Rn then Rn := 6000; elif Type(Rn) eq MonStgElt then Rn := StringToInteger(Rn); end if;
if not assigned Rd then Rd := 96721; elif Type(Rd) eq MonStgElt then Rd := StringToInteger(Rd); end if;
if not assigned Sn then Sn := 3557520; elif Type(Sn) eq MonStgElt then Sn := StringToInteger(Sn); end if;
if not assigned Sd then Sd := 96721; elif Type(Sd) eq MonStgElt then Sd := StringToInteger(Sd); end if;

r0 := Q!Rn/Q!Rd;  s0 := Q!Sn/Q!Sd;
printf "EXACT_START r=%o s=%o\n", r0, s0;

function EKIG(r, s)
    A1 := 2*r*s^2;
    A  := -(9*r*s + 4*r^2 + 4*r + 1)/3;
    B1 := r*s^2*(3*s + 8*r - 2)/3;
    B  := -(54*r^2*s + 81*r*s - 16*r^3 - 24*r^2 - 12*r - 2)/27;
    B2 := r^2;
    return [-24*B1/A1, -12*A, 96*(A/A1)*B1 - 36*B, -4*A1*B2];
end function;

ig := EKIG(r0, s0);
if ig[4] eq 0 then printf "EXACT_DONE degenerate I10=0\n"; quit; end if;

C0 := HyperellipticCurveFromIgusaClebsch(ig);
if BaseRing(C0) ne Q then printf "NO_CURVE_OVER_Q\nEXACT_DONE\n"; quit; end if;

// integral, reduced model
C0 := ReducedMinimalWeierstrassModel(C0);
f0, h0 := HyperellipticPolynomials(C0);
F0 := 4*f0 + h0^2;                       // y^2 = F0 model (up to square 4)
den := LCM([Denominator(c) : c in Coefficients(F0)]);
PZ<x> := PolynomialRing(Q);
F0 := PZ!(F0*den^2);
printf "MODEL y^2 = %o\n", F0;
disc := Integers()!Discriminant(F0);
bad := [ t[1] : t in Factorization(AbsoluteValue(disc)) ];
printf "BAD_PRIMES %o\n", bad;

// ---- twist-sign constraints from chi_p(+-1) mod 31 at good primes
testps := [ p : p in PrimesInInterval(7, 400) |
            p ne 31 and not p in bad and Denominator(r0) mod p ne 0
            and Denominator(s0) mod p ne 0 ];
constraints := [];   // <p, eps> with eps = +1 (chi(1) side) or -1
dead := false;
for p in testps do
    Cp := ChangeRing(HyperellipticCurve(PolynomialRing(GF(p))!F0), GF(p));
    L := LPolynomial(Cp);
    c1 := Evaluate(L, 1) mod 31; cm1 := Evaluate(L, -1) mod 31;
    if c1 eq 0 and cm1 eq 0 then continue; end if;
    if c1 eq 0 then Append(~constraints, <p, 1>);
    elif cm1 eq 0 then Append(~constraints, <p, -1>);
    else
        printf "DEAD_AT_p %o (neither twist has 31 | chi)\n", p;
        dead := true;
        break;
    end if;
    if #constraints ge 30 then break; end if;
end for;
if dead then printf "EXACT_DONE dead\n"; quit; end if;
printf "N_CONSTRAINTS %o\n", #constraints;

// ---- enumerate candidate d: squarefree, support in bad primes + small primes
supp := Sort(Setseq(Seqset(bad cat [ p : p in PrimesInInterval(2, 20) ])));
nsupp := #supp;
cands := [];
for msk in [0 .. 2^nsupp - 1] do
    d := 1;
    for i in [1 .. nsupp] do
        if BitwiseAnd(msk, 2^(i-1)) ne 0 then d := d * supp[i]; end if;
    end for;
    for dd in [d, -d] do
        ok := true;
        for c in constraints do
            if KroneckerSymbol(dd, c[1]) ne c[2] then ok := false; break; end if;
        end for;
        if ok then Append(~cands, dd); end if;
    end for;
end for;
Sort(~cands, func< a, b | AbsoluteValue(a) - AbsoluteValue(b) >);
printf "TWIST_CANDS %o\n", cands;

// ---- exact torsion on surviving twists (usually 0 or 1 candidate)
for d in cands do
    Cd := HyperellipticCurve(d*F0);
    Cd := ReducedMinimalWeierstrassModel(Cd);
    fd, hd := HyperellipticPolynomials(Cd);
    Fd := 4*fd + hd^2;                       // integral y^2 = Fd model
    // strip the largest square integer content (isomorphic curve)
    ct := GCD([Integers()!c : c in Coefficients(Fd)]);
    sq := 1;
    for t in Factorization(ct) do sq *:= t[1]^(2*(t[2] div 2)); end for;
    Fd := PZ![ (Integers()!c) div sq : c in Coefficients(Fd) ];
    printf "TWIST_MODEL d=%o y^2 = %o\n", d, Fd;
    J := Jacobian(HyperellipticCurve(Fd));
    T := TorsionSubgroup(J);
    inv := Invariants(T);
    printf "TORSION d=%o %o\n", d, inv;
    if #inv ge 1 and inv[#inv] mod 31 eq 0 then
        printf "HIT31 d=%o model y^2 = %o\n", d, Fd;
    end if;
end for;
printf "EXACT_DONE\n";
quit;
