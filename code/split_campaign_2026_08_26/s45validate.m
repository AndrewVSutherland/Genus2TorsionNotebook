// INDEPENDENT validation of the [11]-split hit, from integer coefficients
// alone (fresh session; no reuse of the construction objects).
SetColumns(0);
SetMemoryLimit(12*10^9);
Q := Rationals();
Px<x> := PolynomialRing(Q);

// the descended integral model, twist d = 2:
f0 := 1274729070449*x^6 + 13837871743558*x^5 + 28290510535140*x^4 - 199834485454360*x^3 - 889545771580140*x^2 - 379037853127192*x + 1991610457904176;
fC := 2*f0;
C := HyperellipticCurve(fC);
printf "minimizing...\n";
Cm := ReducedMinimalWeierstrassModel(C);
fm, hm := HyperellipticPolynomials(Cm);
printf "MINIMAL MODEL: y^2 + (%o)y = %o\n", hm, fm;
J := Jacobian(Cm);
T := TorsionSubgroup(J);
printf "EXACT TORSION: %o (order %o)\n", Invariants(T), #T;

// conductor-ish: bad primes of the minimal model
dm := Integers()!Discriminant(Cm);
printf "minimal disc: %o\n", Factorization(AbsoluteValue(dm));

// geometric splitness evidence: L-polys at good p should factor over Q(sqrt 11)
// pattern: inert p: P_p(T) = Q1(T^..) hmm; direct: at split p the quartic
// L-poly factors into two quadratics over Q iff a_p pattern; print red count
P<TT> := PolynomialRing(Integers());
nred := 0; ntest := 0;
K<w> := QuadraticField(11);
OK := Integers(K);
for p in PrimesInInterval(13, 150) do
    if dm mod p eq 0 then continue; end if;
    fp := PolynomialRing(GF(p))!fC;
    if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
    Jp := Jacobian(HyperellipticCurve(fp));
    Lp := P!Reverse(Coefficients(EulerFactor(Jp)));
    ntest +:= 1;
    fac := Factorization(Lp);
    if #fac gt 1 or fac[1][2] gt 1 then nred +:= 1; end if;
end for;
printf "L-poly reducible over Q at %o/%o good primes (split-p reducibility expected ~half)\n", nred, ntest;

// the K-factors: rebuild E from the X1(11) point and check isogeny facts
RK<xk> := PolynomialRing(K);
s0 := K!(4/5);
r0 := Roots(xk^2 - (s0^3-3*s0^2+4*s0)*xk + s0)[1][1];
c0 := s0*(r0-1); b0 := r0*c0;
E := EllipticCurve([1-c0, -b0, -b0, 0, 0]);
sig := hom< K -> K | -w >;
Es := EllipticCurve([sig(a) : a in aInvariants(E)]);
printf "E(K) torsion: %o;  j(E) in Q: %o\n",
    Invariants(TorsionSubgroup(E)), jInvariant(E) in Q;
// isogeny test via traces at split primes (isogenous => equal traces at all
// primes); one mismatch is a rigorous negative
mismatch := 0; ntr := 0;
for p in PrimesInInterval(5, 400) do
    dec := Decomposition(OK, p);
    if #dec ne 2 then continue; end if;
    pr := dec[1][1];
    okr := true; t1 := 0; t2 := 0;
    try
        t1 := TraceOfFrobenius(Reduction(E, pr));
        t2 := TraceOfFrobenius(Reduction(Es, pr));
    catch e; okr := false; end try;
    if not okr then continue; end if;
    ntr +:= 1;
    if t1 ne t2 then mismatch +:= 1; end if;
end for;
printf "trace mismatches E vs E^sigma at same split prime: %o/%o (>0 => NOT isogenous => J Q-simple)\n", mismatch, ntr;
printf "E has CM: %o\n", HasComplexMultiplication(E);

// direct count identities at fresh primes (not used in the construction)
for p in [101, 103, 107, 109] do
    if dm mod p eq 0 then continue; end if;
    fp := PolynomialRing(GF(p))!fC;
    if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
    nj := #Jacobian(HyperellipticCurve(fp));
    dec := Decomposition(OK, p);
    if #dec eq 2 then
        pr := #Reduction(E, dec[1][1]) * #Reduction(Es, dec[2][1]);
        printf "p=%o split: #J=%o  #E*#Es=%o  match=%o\n", p, nj, pr, nj eq pr;
    else
        pr := #Reduction(E, dec[1][1]);
        printf "p=%o inert: #J=%o  #E(Fp^2)=%o  match=%o\n", p, nj, pr, nj eq pr;
    end if;
end for;
printf "S45VALIDATE_DONE\n";
quit;
