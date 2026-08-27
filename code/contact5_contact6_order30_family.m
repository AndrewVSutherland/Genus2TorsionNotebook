//////////////////////////////////////////////////////////////////////
//  Explicit one-parameter simultaneous contact-5/contact-6 family.
//
//  Parameter R gives a point (t,Y) on
//
//      Y^2 = 5*t^2 - 6*t + 5
//
//  by the line through (5,10):
//
//      t = (5R^2 - 20R + 19)/(R^2 - 5),
//      Y = -2*(5R^2 - 22R + 25)/(R^2 - 5).
//
//  Put u=t^3 and
//
//      s = t^5 + t^4 + (5/2)t^3 + (1/2)t
//          +/- t*(t-1/2)*(t+1)*Y.
//
//  This lies on the genus-zero core of the simultaneous contact equations.
//  The remaining coefficient q=A-e is recovered from the linear equation
//  F2 - s*F3 = 0.  The resulting curve
//
//      C_R: y^2 = h6^2 - (x-1)^6 = h5^2 - K*x^5
//
//  has rational 6- and 5-torsion contact divisors.  Smooth nonboundary
//  specializations generically have a point of order 30.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned output_file then
    output_file := "data/contact5_contact6_order30_family_samples.txt";
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

function IntegralModel(f)
    L := Z!1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function TorsionOrder(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return &*invs;
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199] do
        try
            fp := ChangeRing(f, GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
                continue;
            end if;
            Lp := LPolynomial(ChangeRing(C, GF(p)));
            fac := Factorization(Lp);
            if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
                return true, p, Lp;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false, 0, P!0;
end function;

function Family(R, branch)
    denR := R^2 - 5;
    if denR eq 0 then
        return false, _, _, _, _, _, _, _, _, _;
    end if;
    t := (5*R^2 - 20*R + 19)/denR;
    Y := -2*(5*R^2 - 22*R + 25)/denR;
    u := t^3;
    s := t^5 + t^4 + (Q!5/2)*t^3 + (Q!1/2)*t
        + branch*t*(t - Q!1/2)*(t + 1)*Y;

    Cc := (u^2 + 1)/(2*u);
    c := (u^2 - 1)/(2*u);
    denq := u^6 + 6*u^4*s - 2*u^4 + 15*u^3*s - u*s^3 + u^2;
    numq := 15*u^5 + 90*u^4 + 20*u^3*s - 6*u^2*s^2 + 231*u^3
        + 2*u^2*s - 15*u*s^2 + 90*u^2 - 20*u*s + 15*u - 2*s;
    if u eq 0 or c eq 0 or denq eq 0 then
        return false, _, _, _, _, _, _, _, _, _;
    end if;
    q := numq/denq;
    A := (s + q)/2;
    e := (s - q)/2;
    B := (15 - s*q)/2;
    d := (B*Cc + 3)/c;
    h6 := x^3 + A*x^2 + B*x + Cc;
    h5 := e*x^2 + d*x + c;
    f := h6^2 - (x-1)^6;
    return true, t, Y, u, s, q, h6, h5, f, ExactQuotient(h5^2 - f, x^5);
end function;

samples := [Q!0, Q!1, Q!2, Q!3, Q!4, Q!6, Q!-1, Q!1/2, Q!5/2];
out := Open(output_file, "w");
fprintf out, "# contact5/contact6 order30 family samples\n";
fprintf out, "# t=(5R^2-20R+19)/(R^2-5), Y=-2*(5R^2-22R+25)/(R^2-5)\n";

for Rval in samples do
    for branch in [-1, 1] do
        ok, t, Y, u, s, q, h6, h5, f, K := Family(Rval, branch);
        fprintf out, "============================================================\n";
        fprintf out, "R=%o branch=%o ok=%o\n", Rval, branch, ok;
        if not ok then
            continue;
        end if;
        fprintf out, "t=%o Y=%o u=%o s=%o q=%o K=%o\n", t, Y, u, s, q, K;
        fprintf out, "h6=%o\n", h6;
        fprintf out, "h5=%o\n", h5;
        fprintf out, "f=%o\n", f;
        fprintf out, "conic_check=%o contact5_check=%o\n",
                Y^2 - (5*t^2 - 6*t + 5), h5^2 - f - K*x^5;
        if Degree(f) ne 5 or Discriminant(f) eq 0 then
            fprintf out, "BAD degree=%o disc=%o\n", Degree(f), Discriminant(f);
            continue;
        end if;
        Ccurve := HyperellipticCurve(f);
        J := Jacobian(Ccurve);
        D5 := J![x, Evaluate(h5,0)];
        D6 := J![x-1, Evaluate(h6,1)];
        fprintf out, "orders D5=%o D6=%o D5+D6=%o\n", Order(D5), Order(D6), Order(D5+D6);
        fI, scale := IntegralModel(f);
        G, phi := TorsionSubgroup(Jacobian(HyperellipticCurve(fI)));
        invs := Invariants(G);
        simple, pcert, Lp := SimpleCertificate(fI);
        fprintf out, "integral_scale=%o\n", scale;
        fprintf out, "integral_f=%o\n", fI;
        fprintf out, "torsion_invariants=%o torsion_order=%o\n", invs, TorsionOrder(invs);
        fprintf out, "simple_certificate=%o p=%o Lp=%o\n", simple, pcert, Lp;
    end for;
end for;

delete out;
print "Wrote", output_file;
quit;
