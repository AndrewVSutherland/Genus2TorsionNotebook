// claude_ov_erratum_control.m -- genus2red ERRATUM propagation (2026-07-25).
//
// CONTROL for the conductor recomputation.  Lane 7 discovered that PARI's
// genus2red omits the 2-part of the conductor (it is documented for p > 2 and
// returns the sentinel exponent -1 at 2).  Before rewriting any recorded
// conductor we must show that the REPLACEMENT method -- Magma's Conductor on
// the reduced minimal Weierstrass model -- reproduces independently published
// conductors exactly, INCLUDING the 2-part.
//
// The control set is 14 LMFDB genus-2 curves taken verbatim from the rows of
// paper/torsion_realizations.tex; the published conductor is the first
// component of the LMFDB label, so it is an independent datum (Booker-
// Sijsling-Sutherland-Voight / Booker-Sutherland), not something this project
// computed.  v_2 of the published conductor ranges over 0,1,2,3,4.
//
// Usage:   magma -b code/claude_ov_erratum_control.m
// Markers: CTRL / CTRLSUMMARY / ERRATUM_CONTROL_DONE
SetColumns(0);
SetMemoryLimit(4*10^9);
Q := Rationals(); P<x> := PolynomialRing(Q); Z := Integers();

// <label, f, h, published conductor>
ctrl := [*
  <"295.a.295.2",     x^5-40*x^3+22*x^2+389*x-608,          x^2+x+1,  295>,
  <"464.a.464.1",     -x^6-2*x^5-2*x^4-x^3,                 x+1,      464>,
  <"464.a.29696.2",   4*x^5+33*x^4+72*x^3+16*x^2+x,         x,        464>,
  <"1416.b.135936.1", -2*x^4-x^3+x+1,                       x^3+x,    1416>,
  <"2600.a.338000.1", 10*x^5+8*x^4-5*x^3-3*x^2+x,           x,        2600>,
  <"1575.c.1",        -27*x^5+23*x^4+8*x^3-4*x^2-x,         x^2+x,    1575>,
  <"3942.b.3",        9*x^5+20*x^4-3*x^3-6*x^2-x,           x^2+x,    3942>,
  <"9450.b.2",        -225*x^5+146*x^4-16*x^3-6*x^2+x,      x^2+x,    9450>,
  <"10512.n.1",       -12*x^5+6*x^4+12*x^3-7*x^2,           x^2+1,    10512>,
  <"12300.e.2",       x^6-7*x^5+8*x^4+16*x^3-5*x^2-5*x,     x^2+x,    12300>,
  <"19044.h.2",       x^6-3*x^5+9*x^4-5*x^3+12*x^2-6*x,     x^2+x,    19044>,
  <"39600.dq.1",      -6*x^5-14*x^4+18*x^3+11*x^2-12*x+2,   x^2+1,    39600>,
  <"8136.c.1",        9*x^6-3*x^5-12*x^4+3*x^3+2*x^2,       x^2+1,    8136>,
  <"28200.e.1",       9*x^6-3*x^5-42*x^4+34*x^3+4*x^2-3*x,  x^2+1,    28200>
*];

nok := 0; nbad := 0; noddok := 0;
for r in ctrl do
    lab := r[1]; f := r[2]; h := r[3]; Npub := r[4];
    C  := HyperellipticCurve(f, h);
    Cm := ReducedMinimalWeierstrassModel(C);
    Dm := Z!Discriminant(Cm);
    v2D := Valuation(AbsoluteValue(Dm), 2);
    N := Conductor(Cm);
    agree := (N eq Npub);
    oddagree := ((N div 2^Valuation(N,2)) eq (Npub div 2^Valuation(Npub,2)));
    if agree then nok +:= 1; else nbad +:= 1; end if;
    if oddagree then noddok +:= 1; end if;
    printf "CTRL %o  Npub=%o=%o  Nmagma=%o=%o  v2(Npub)=%o  v2(disc_min)=%o  ogg_guaranteed=%o  MATCH=%o  oddmatch=%o\n",
        lab, Npub, Factorization(Npub), N, Factorization(N),
        Valuation(Npub,2), v2D, v2D lt 12, agree, oddagree;
end for;
printf "CTRLSUMMARY  matched=%o  mismatched=%o  oddpart_matched=%o  of %o\n", nok, nbad, noddok, #ctrl;
print "ERRATUM_CONTROL_DONE";
quit;
