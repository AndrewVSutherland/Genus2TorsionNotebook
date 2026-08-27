// neighborhood.m — Lane D task 2: the (2,2,2,12) hit curve's isogeny
// neighborhood.  RichelotIsogenousSurfaces + TwoPowerIsogenies codomain
// torsion, bad primes of the minimal model vs the forced-bad-reduction
// theorem (must be bad at 3,5,7,11,13), G2-dedupe against the hit.
P<x> := PolynomialRing(Rationals());
// The hit: quintic model from (s,m,n)=(336396,-689185,-166464) on M(2,2,2,6)
A := [1,1,1,2,2];
B := [282322361376, -8243383980, -64241207724, -114724491840, 561915878400];
f0 := &*[A[i] + B[i]*x : i in [1..5]];
C0 := HyperellipticCurve(f0);
C := ReducedMinimalWeierstrassModel(C0);
fC, hC := HyperellipticPolynomials(C);
printf "minimal model: y^2 + (%o)y = %o\n", hC, fC;
Cs := IntegralModel(SimplifiedModel(C));
J := Jacobian(Cs);
inv0 := Invariants(TorsionSubgroup(J));
printf "HIT TORSION (sanity) = %o\n", inv0;
g2hit := G2Invariants(C);

// ---- bad primes ----
// candidate support: primes of the quintic-model disc = primes of the B_i and
// of the pairwise resultants A_i B_j - A_j B_i (roots -A_i/B_i differences)
cand := {};
for i in [1..5] do
    for q in PrimeDivisors(B[i]) do Include(~cand, q); end for;
    for j in [i+1..5] do
        r := A[i]*B[j] - A[j]*B[i];
        if r ne 0 then for q in PrimeDivisors(AbsoluteValue(r)) do Include(~cand, q); end for; end if;
    end for;
end for;
dmin := Integers()!Discriminant(C);
bad := [];
rem := AbsoluteValue(dmin);
for q in Sort(SetToSequence(cand)) do
    v := Valuation(rem, q);
    if v gt 0 then Append(~bad, <q, v>); rem := rem div q^v; end if;
end for;
printf "minimal disc = %o\n", dmin;
printf "bad primes (p, v_p(disc_min)) = %o\n", bad;
printf "unaccounted cofactor of |disc_min| = %o (should be 1)\n", rem;
forced := {3,5,7,11,13};
have := {t[1] : t in bad};
printf "forced-bad-reduction check {3,5,7,11,13} subset bad: %o\n", forced subset have;

// ---- torsion helper on a curve/jacobian with reduction first ----
function TorsInv(Cc)
    Cr := Cc;
    try
        Cr := ReducedMinimalWeierstrassModel(Cc);
    catch e
        Cr := Cc;
    end try;
    Ci := IntegralModel(SimplifiedModel(Cr));
    return Invariants(TorsionSubgroup(Jacobian(Ci))), Cr;
end function;

// ---- Richelot neighbors ----
print "=== RichelotIsogenousSurfaces ===";
rs := RichelotIsogenousSurfaces(J);
printf "count = %o\n", #rs;
for i in [1..#rs] do
    t := Type(rs[i]);
    if t eq JacHyp then
        D := Curve(rs[i]);
        invR, Dr := TorsInv(D);
        printf "richelot_%o: JacHyp torsion %o order %o G2same %o\n",
            i, invR, &*invR, G2Invariants(D) eq g2hit;
        fR, hR := HyperellipticPolynomials(Dr);
        printf "   model: y^2+(%o)y = %o\n", hR, fR;
    else
        printf "richelot_%o: type %o (non-Jacobian codomain!) %o\n", i, t, rs[i];
    end if;
end for;

// ---- TwoPower neighbors ----
print "=== TwoPowerIsogenies ===";
Js, prods, weils := TwoPowerIsogenies(J);
printf "jacobians %o products %o weil_restrictions %o\n", #Js, #prods, #weils;
for i in [1..#Js] do
    D := Curve(Js[i]);
    invT, Dr := TorsInv(D);
    printf "twopower_%o: torsion %o order %o G2same %o\n",
        i, invT, &*invT, G2Invariants(D) eq g2hit;
    fT, hT := HyperellipticPolynomials(Dr);
    printf "   model: y^2+(%o)y = %o\n", hT, fT;
end for;
for i in [1..#prods] do printf "product_%o: %o\n", i, prods[i]; end for;
for i in [1..#weils] do printf "weil_%o: %o\n", i, weils[i]; end for;
print "NEIGHBORHOOD DONE";
quit;
