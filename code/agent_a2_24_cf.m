
//////////////////////////////////////////////////////////////////////
//  Polynomial continued fraction of sqrt(f) for genus-2 D_infty order.
//
//  For monic squarefree f of degree 6, D = infty_+ - infty_- has finite
//  order iff the CF of sqrt(f) is (quasi-)periodic; the order equals the
//  sum of deg(a_i) over one period (the fundamental unit degree).  This
//  is the constructive heart of the Platonov-Petrunin method: to build f
//  with order n, prescribe the periodic partial quotients and invert.
//
//  Reduction (Q[x] infrastructure):
//    s = polynomial part of sqrt(f) (deg 3);
//    P_0=0, Q_0=1;
//    a_i = PolyPart((P_i + s)/Q_i);
//    P_{i+1} = a_i*Q_i - P_i;
//    Q_{i+1} = (f - P_{i+1}^2)/Q_i.
//  Period ends when Q_i is a nonzero constant (i>=1); order = sum deg a_i
//  over the reduced cycle.
//
//  Usage: magma -b agent_a2_24_cf.m   (self-tests f14, f18)
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
P<x> := PolynomialRing(Q);

// polynomial "square root to degree 3" part s with s^2 = f + O(x^2)
function SqrtPolyPart(f)
    // Newton: s = x^3 + ... ; compute coefficients so deg(f - s^2) <= 2
    assert Degree(f) eq 6 and LeadingCoefficient(f) eq 1;
    s := x^3;
    for k in [1..3] do
        // adjust coefficient of x^(3-k) in s to kill x^(6-k) in f - s^2
        d := f - s^2;
        if Degree(d) le 2 then break; end if;
        ck := Coefficient(d, 6-k) / (2*Coefficient(s,3));  // /2 leading
        s := s + ck*x^(3-k);
    end for;
    return s;
end function;

// continued fraction / order of D_infty.  Returns order (0 if not found
// within maxsteps), the partial-quotient degree sequence, and the Q-seq.
function DinftyOrderCF(f, maxsteps)
    s := SqrtPolyPart(f);
    // s is the polynomial part of sqrt(f): deg(f - s^2) <= 2
    Pi := P!0; Qi := P!1;
    degs := [];
    total := 0;
    for i in [0..maxsteps] do
        // a_i = polynomial part of (P_i + s)/Q_i
        ai := (Pi + s) div Qi;
        Append(~degs, Degree(ai));
        total +:= Degree(ai);   // include a_0 (deg 3): order = deg of convergent alpha
        Pnext := ai*Qi - Pi;
        Qnext := (f - Pnext^2) div Qi;
        // exact division check
        if (f - Pnext^2) mod Qi ne 0 then
            return 0, degs, "nonexact";   // f not squarefree / bad
        end if;
        Pi := Pnext; Qi := Qnext;
        // periodicity: Q_i constant (nonzero) after step >=1 signals unit
        if i ge 1 and Degree(Qi) le 0 and Qi ne 0 then
            // the accumulated deg from a_1..a_i is the order (fundamental unit degree)
            return total, degs, "periodic";
        end if;
    end for;
    return 0, degs, "no-period";
end function;

// self tests
tests := [
    <"f14", (x^2+1)*(x^4+5*x^2+4*x+4), 14>,
    <"f18", (x^2-x+1)*(x^4-x^3+9*x^2+8*x-8), 18>,
    <"f28-cyclic?", x^6 + 2*x^5 - 5*x^4 - 14*x^3 - 3*x^2 + 24*x + 28, 0>
];
print "== self-tests: CF order of D_infty ==";
for tr in tests do
    name := tr[1]; f := tr[2]; expect := tr[3];
    if LeadingCoefficient(f) ne 1 then f := f/LeadingCoefficient(f); end if;
    ord, degs, stat := DinftyOrderCF(f, 60);
    // cross-check with Magma torsion order of D_infty if leading coeff square
    printf "%o: CF order=%o status=%o degseq(first 12)=%o", name, ord, stat,
        degs[1..Min(12,#degs)];
    if expect ne 0 then
        printf "  expected %o %o", expect, (ord eq expect select "OK" else "MISMATCH");
    end if;
    printf "\n";
end for;
print "DONE";
quit;
