// diagnostic: does SOME quadratic twist of the Mestre model have 11-torsion?
// For each good p compute the L-poly P_p(T) of Jac(C0) mod p and test
// P_p(1), P_p(-1) mod 11.  Every p must admit a valid sign; the sign pattern
// chi_d(p) then determines the twist d.  Also factor the model disc for the
// twist support.  Also independent check: compare #J(F_p) (right sign)
// against the Weil-restriction prediction from E/K.
SetColumns(0);
SetMemoryLimit(12*10^9);
Q := Rationals();
K<w> := QuadraticField(11);
RK<xk> := PolynomialRing(K);
Px<x> := PolynomialRing(Q);
s0 := K!(4/5);
r0 := Roots(xk^2 - (s0^3-3*s0^2+4*s0)*xk + s0)[1][1];
c0 := s0*(r0-1); b0 := r0*c0;
E := EllipticCurve([1-c0, -b0, -b0, 0, 0]);

// the Mestre Q-model from s45twist.log (verbatim coefficients)
f0 := 293905425800403553300959700311514131054287952646953835737528803650981*x^6 + 2877353601236711333671754819708022309227834213337716399097864756854719*x^5 + 4017856528250105355617761697782578648580470297035966457218930861007895*x^4 - 23017157850552281779223650581854187101214773529693805862364667346693940*x^3 + 1554341726238583122216766912112592053687155075346947132712543311925420*x^2 + 35937488045837908039588077521888587921820365339980599384951481645695044*x + 23041739280560127344150227902781441051504782314704526736635612740967756;
dsc := Integers()!Discriminant(f0);
td := TrialDivision(AbsoluteValue(dsc), 10^6);
printf "disc trial division: %o\n", td;

goodp := [ p : p in PrimesInInterval(13, 120) | p ne 11 and dsc mod p ne 0 ];
printf "testing %o good primes\n", #goodp;
OK := Integers(K);
for p in goodp[1..15] do
    Fp := GF(p);
    fp := PolynomialRing(Fp)!f0;
    if Degree(fp) lt 5 or not IsSquarefree(fp) then continue; end if;
    Jp := Jacobian(HyperellipticCurve(fp));
    np1 := #Jp;                       // = P_p(1)
    // P_p(-1) = #J of the nontrivial quadratic twist
    P<T> := PolynomialRing(Integers());
    Lp := P!Reverse(Coefficients(EulerFactor(Jp)));
    npm := Integers()!Evaluate(Lp, -1);
    ok1 := np1 mod 11 eq 0; okm := npm mod 11 eq 0;
    // Weil restriction prediction
    pred := "?";
    dec := Decomposition(OK, p);
    if #dec eq 2 then
        n1 := 0; n2 := 0; okr := true;
        try
            E1p := Reduction(E, dec[1][1]);  n1 := #E1p;
            E2p := Reduction(E, dec[2][1]);  n2 := #E2p;
        catch e; okr := false; end try;
        if okr then pred := Sprintf("split: n1*n2=%o (mod11=%o)", n1*n2, (n1*n2) mod 11); end if;
    else
        // inert: #A(F_p) = #E(F_{p^2}) for the restriction
        okr := true; nE := 0;
        try nE := #Reduction(E, dec[1][1]); catch e; okr := false; end try;
        if okr then pred := Sprintf("inert: #E(F_p2)=%o (mod11=%o)", nE, nE mod 11); end if;
    end if;
    printf "p=%o: P(1)=%o (11|%o)  P(-1)=%o (11|%o)  | RES pred %o\n",
        p, np1, ok1, npm, okm, pred;
end for;
printf "S45DIAG_DONE\n";
quit;
