//////////////////////////////////////////////////////////////////////
//  !!! UNRELIABLE / RECORDED DEAD-END (2026-07-09) !!!
//  This cheap point-count probe OVER-COUNTS: the residual-GCD common-M test
//  admits spurious common roots that are NOT genuine order-3 classes even
//  after filtering M!=0 and kappa=a-M^2!=0.  Symptom: erratic counts across
//  consecutive primes (Delta=2 gave 64,28,28,436 over 7,11,13,17 -- a genuine
//  curve count must vary smoothly ~ell).  Counting the disc=Delta fiber
//  CORRECTLY requires building the Jacobian point [q3, H mod q3] over F_ell
//  and verifying 3D=0, D!=0 (the ContactClassesFF machinery in
//  agent_a2_24_wsplit_contact_cover.m, Mode:=ff_disc).  Kept as a documented
//  negative so the next agent does not repeat the GCD shortcut.
//
//  [2,24] W-split contact cover: dimension + genus PROBE of disc(q3)=Delta
//  fibers, via finite-field point counts (no exact elimination).
//
//  Gate for the "descent vs accept-milestone" decision
//  (notes/agent_a2_24_wsplit_contact_cover.md).  For fixed Delta, the fiber
//  is  { W(r,p,t,beta)=0 and residuals R2=R1=R0=0 have a common M }  in the
//  free coordinates (r,t,beta,U) [p from W, M from residuals].  Count
//  #fiber(F_ell) across primes:
//    - scaling ~ell  => fiber is a CURVE (1-dim);  ~ell^2 => surface.
//    - if a curve, |#C - (ell+1)|/(2*sqrt(ell)) is a lower bound on genus g.
//  Small g on some Delta => descent/Chabauty feasible (Option 2).
//  All populated Delta high-g => accept the simple-Z/24 milestone (Option 3).
//
//  Reuses TriangularSystem VERBATIM from agent_a2_24_wsplit_contact_cover.m.
//  Usage: magma -b Primes:="7,11,13,17,19,23" Deltas:="2,3,5,6,7" \
//         code/agent_a2_24_fiber_genus_probe.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB := 2; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
if not assigned Primes then Primes := "7,11,13,17,19,23"; end if;
if not assigned Deltas then Deltas := "2,3,5,6,7"; end if;
plist := [StringToInteger(s) : s in Split(Primes, ",")];
dlist := [StringToInteger(s) : s in Split(Deltas, ",")];

function TriangularSystem(R, DeltaVal)
    r := R.1; p := R.2; t := R.3; beta := R.4; U := R.5; M := R.6;
    F := FieldOfFractions(R);
    P<X> := PolynomialRing(F);
    rF := F!r; pF := F!p; tF := F!t; bF := F!beta;
    UF := F!U; MF := F!M; DeltaF := F!DeltaVal;
    VF := (UF^2 - DeltaF)/4;
    e := tF^2 - 2*pF*tF/rF;
    d := e + 2*pF - rF^2;
    lambda := rF/tF;
    a := rF^2 - lambda;
    b := 2*rF*pF - 2*lambda*(pF + rF*tF) + 2*rF*lambda;
    c := pF^2 + 2*pF*rF^2 - rF^4 - rF^3*tF - rF*pF^2/tF
         - lambda*(rF^2 + e)
         + 2*lambda*(rF*pF + rF^2*tF - 3*pF*tF + rF*tF^2);
    Qpoly := X^2 + d;
    q := a*X^2 + b*X + c;
    f := q*(Qpoly^2 + q);
    fc := [Coefficient(f, i) : i in [0..6]];
    K := fc[7] - MF^2;
    H2 := (fc[6] - 3*K*UF)/(2*MF);
    H1 := (fc[5] - H2^2 - 3*K*(UF^2 + VF))/(2*MF);
    H0 := (fc[4] - 2*H2*H1 - K*(UF^3 + 6*UF*VF))/(2*MF);
    R2 := fc[3] - (H1^2 + 2*H2*H0) - K*(3*UF^2*VF + 3*VF^2);
    R1 := fc[2] - 2*H1*H0 - K*(3*UF*VF^2);
    R0 := fc[1] - H0^2 - K*VF^3;
    Wexpr := bF^2*(2*d + a + bF^2)^2 - 4*(d^2 + c)*bF^2 - b^2;
    return [R!Numerator(R2), R!Numerator(R1), R!Numerator(R0)], R!Numerator(Wexpr);
end function;

print "FIBER GENUS PROBE  primes:", plist, " deltas:", dlist;
for Delta in dlist do
    printf "\n=== disc(q3) = Delta = %o ===\n", Delta;
    printf "%-5o %-8o %-10o %-10o\n", "ell", "#fiber", "#/ell", "g_lb";
    for ell in plist do
        Fl := GF(ell);
        Rl := PolynomialRing(Fl, 6);
        r := Rl.1; p := Rl.2; t := Rl.3; beta := Rl.4; U := Rl.5; M := Rl.6;
        resids, W := TriangularSystem(Rl, Delta);
        Pm := PolynomialRing(Fl); mm := Pm.1;
        // W depends only on (r,p,t,beta); solve for p per (r,t,beta)
        cnt := 0;
        for rr in Fl do
            if rr eq 0 then continue; end if;
            for tt in Fl do
                if tt eq 0 or rr*tt eq 1 then continue; end if;
                for bb in Fl do
                    // W(rr,p,tt,bb) as a poly in p
                    Wp := Evaluate(W, [Rl!rr, p, Rl!tt, Rl!bb, Rl!0, Rl!0]);
                    // coerce to univariate in p (slot 2)
                    Wpu := Pm ! [MonomialCoefficient(Wp, p^i) : i in [0..Degree(Wp,p)]];
                    if Wpu eq 0 then continue; end if;
                    proots := [rt[1] : rt in Roots(Wpu)];
                    for pp in proots do
                        // a = chart leading coeff = fc[7]; kappa = a - M^2
                        aa := rr^2 - rr/tt;
                        for uu in Fl do
                            gs := [];
                            ok := true;
                            for R in resids do
                                Re := Evaluate(R, [Rl!rr, Rl!pp, Rl!tt, Rl!bb, Rl!uu, M]);
                                Reu := Pm ! [MonomialCoefficient(Re, M^i) : i in [0..Degree(Re,M)]];
                                if Reu eq 0 then ok := false; break; end if;   // identically 0: skip (degenerate)
                                Append(~gs, Reu);
                            end for;
                            if not ok then continue; end if;
                            g := gs[1];
                            for i in [2..#gs] do g := GCD(g, gs[i]); end for;
                            if Degree(g) ge 1 then
                                for rt in Roots(g) do
                                    m0 := rt[1];
                                    if m0 eq 0 then continue; end if;        // H3=0 degenerate
                                    if aa - m0^2 eq 0 then continue; end if; // kappa=0 spurious (perfect square)
                                    cnt +:= 1;                               // genuine contact point
                                end for;
                            end if;
                        end for;
                    end for;
                end for;
            end for;
        end for;
        glb := Abs(cnt - (ell+1)) / (2*Sqrt(RealField(6)!ell));
        printf "%-5o %-8o %-10.3o %-10.2o\n", ell, cnt, RealField(6)!cnt/ell, glb;
    end for;
end for;
print "DONE";
quit;
