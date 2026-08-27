//////////////////////////////////////////////////////////////////////
//  [2,24] W-split contact cover: genus of the U=0, beta-FREE curve slice.
//
//  notes/agent_a2_24_wsplit_contact_cover.md "Next exact move" #2:
//  "stop freezing both U and beta, and instead analyze a one-parameter
//   family such as Delta=3, U=0 with beta free".
//
//  Curve C = { R2=R1=R0=0, Wsplit=0 } in (r,p,t,beta,M), with U=0 and
//  V=(0-Delta)/4 fixed.  4 equations, 5 unknowns => 1-dimensional.
//  Eliminate p (via Wsplit), then M, then t -> plane model in (r,beta);
//  strip boundary/spurious factors (GCD) at each stage; report genus of
//  every substantial factor.  (No known rational point exists on C, so we
//  cannot select "the component through a base point"; we report all.)
//
//  Usage: magma -b Delta:=3 code/agent_a2_24_wsplit_betafree_genus.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB := 4; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);
Q := Rationals();
if not assigned Delta then Delta := 3; elif Type(Delta) eq MonStgElt then Delta := StringToInteger(Delta); end if;

// VERBATIM TriangularSystem from agent_a2_24_wsplit_contact_cover.m
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

R<r,p,t,beta,U,M> := PolynomialRing(Q, 6);
resids, Wsplit := TriangularSystem(R, Delta);
// impose U = 0
sub0 := [r,p,t,beta,R!0,M];
c2 := Evaluate(resids[1], sub0); c1 := Evaluate(resids[2], sub0);
c0 := Evaluate(resids[3], sub0); W := Evaluate(Wsplit, sub0);
printf "Delta=%o U=0 : deg c2=%o c1=%o c0=%o W=%o (Wdeg_p=%o)\n",
    Delta, TotalDegree(c2), TotalDegree(c1), TotalDegree(c0),
    TotalDegree(W), Degree(W, p);

// helper: strip GCD of two polys then Resultant in variable v
function ResStrip(A, B, v)
    g := GCD(A, B);
    if TotalDegree(g) ge 1 then A := A div g; B := B div g; end if;
    return Resultant(A, B, v);
end function;

// stage 1: eliminate p using W (quadratic in p) against each residual
t0 := Cputime();
G2 := ResStrip(W, c2, p); G1 := ResStrip(W, c1, p); G0 := ResStrip(W, c0, p);
printf "elim p [%.1os]: deg G2=%o G1=%o G0=%o\n", Cputime(t0),
    TotalDegree(G2), TotalDegree(G1), TotalDegree(G0);

// stage 2: eliminate M from the (r,t,beta,M) polys
t0 := Cputime();
Ha := ResStrip(G2, G0, M); Hb := ResStrip(G1, G0, M);
printf "elim M [%.1os]: deg Ha=%o Hb=%o\n", Cputime(t0),
    TotalDegree(Ha), TotalDegree(Hb);

// stage 3: eliminate t -> plane curve in (r,beta)
t0 := Cputime();
Plane := ResStrip(Ha, Hb, t);
printf "elim t [%.1os]: plane deg=%o\n", Cputime(t0), TotalDegree(Plane);

P2<a2,b2> := PolynomialRing(Q, 2);
h := hom< R -> P2 | [a2, 0, 0, b2, 0, 0] >;   // r->a2, beta->b2, others 0
plane2 := h(Plane);
printf "plane in (r,beta): deg=%o\n", TotalDegree(plane2);
if plane2 eq 0 then print "PLANE VANISHED"; print "DONE"; quit; end if;
fac := Factorization(plane2);
printf "factors (deg,mult): %o\n", [<TotalDegree(g[1]), g[2]> : g in fac];
for g in fac do
    G := g[1];
    if TotalDegree(G) lt 3 then continue; end if;
    printf "--- factor deg %o mult %o : genus ---\n", TotalDegree(G), g[2];
    try
        C := Curve(AffineSpace(P2), G);
        printf ">>> GENUS = %o\n", Genus(C);
    catch ee printf "genus error: %o\n", ee`Object; end try;
end for;
print "DONE";
quit;
