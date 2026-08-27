//////////////////////////////////////////////////////////////////////
//  Finite-prime diagnostic for extra independent 3-torsion on the
//  contact5/contact6 order-30 family.
//
//  The family already has one rational 3-torsion class, coming from
//  the order-6 contact divisor.  A second independent rational 3-class
//  forces 9 | #J(F_p) at every good p != 3.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned output_file then
    output_file := "data/contact5_contact6_order30_extra3_residue_scan.txt";
end if;

Z := Integers();

function FamilyFinite(R, eps, F)
    PF<x> := PolynomialRing(F);
    denR := R^2 - F!5;
    if denR eq 0 then
        return false, PF!0;
    end if;
    t := (5*R^2 - 20*R + 19)/denR;
    Y := -2*(5*R^2 - 22*R + 25)/denR;
    u := t^3;
    if u eq 0 then
        return false, PF!0;
    end if;
    s := t^5 + t^4 + (F!5/F!2)*t^3 + (F!1/F!2)*t
        + (F!eps)*t*(t - F!1/F!2)*(t + 1)*Y;

    Cc := (u^2 + 1)/(2*u);
    c := (u^2 - 1)/(2*u);
    denq := u^6 + 6*u^4*s - 2*u^4 + 15*u^3*s - u*s^3 + u^2;
    if c eq 0 or denq eq 0 then
        return false, PF!0;
    end if;
    numq := 15*u^5 + 90*u^4 + 20*u^3*s - 6*u^2*s^2 + 231*u^3
        + 2*u^2*s - 15*u*s^2 + 90*u^2 - 20*u*s + 15*u - 2*s;
    q := numq/denq;
    A := (s + q)/2;
    B := (15 - s*q)/2;
    h6 := x^3 + A*x^2 + B*x + Cc;
    f := h6^2 - (x-1)^6;
    return true, f;
end function;

out := Open(output_file, "w");
fprintf out, "# extra independent 3-torsion residue diagnostic\n";
fprintf out, "# necessary condition at good p != 3: 9 divides #J(F_p)\n\n";

for p in [7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73] do
    F := GF(p);
    fprintf out, "p=%o\n", p;
    for eps in [-1, 1] do
        good := 0;
        bad := 0;
        allowed9 := 0;
        allowed27 := 0;
        examples := [];
        for a in [0..p-1] do
            ok, f := FamilyFinite(F!a, eps, F);
            if not ok or Degree(f) ne 5 or Discriminant(f) eq 0 then
                bad +:= 1;
                continue;
            end if;
            good +:= 1;
            C := HyperellipticCurve(f);
            N := Z!Evaluate(LPolynomial(C), 1);
            if N mod 9 eq 0 then
                allowed9 +:= 1;
                if #examples lt 8 then
                    Append(~examples, <a,N>);
                end if;
            end if;
            if N mod 27 eq 0 then
                allowed27 +:= 1;
            end if;
        end for;
        fprintf out, "  eps=%o good=%o bad=%o allowed9=%o allowed27=%o examples=%o\n",
                eps, good, bad, allowed9, allowed27, examples;
    end for;
end for;

delete out;
print "Wrote", output_file;
quit;
