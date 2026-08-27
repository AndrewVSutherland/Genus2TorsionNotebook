//////////////////////////////////////////////////////////////////////
//  Rational-point diagnostics for the two genus-3 extra-2 square
//  conditions in the contact5/contact6 order-30 family.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned output_file then
    output_file := "data/contact5_contact6_order30_extra2_rational_points.txt";
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

curves := [
    <"minus", -(3*x^2 - 12*x + 13)*
        (365*x^6 - 4044*x^5 + 18249*x^4 - 42664*x^3
         + 54039*x^2 - 34764*x + 8755)>,
    <"plus", -(7*x^2 - 32*x + 37)*
        (305*x^6 - 3504*x^5 + 16761*x^4 - 42784*x^3
         + 61611*x^2 - 47664*x + 15595)>
];

function IsRationalSquare(q)
    if q lt 0 then
        return false, Q!0;
    end if;
    n := Numerator(q);
    d := Denominator(q);
    sn := IsSquare(n);
    sd := IsSquare(d);
    if sn and sd then
        return true, Sqrt(Q!q);
    end if;
    return false, Q!0;
end function;

function RationalSearch(f, H)
    pts := [];
    seen := {};
    for den in [1..H] do
        for num in [-H..H] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            r := Q!num/den;
            key := Sprint(r);
            if key in seen then
                continue;
            end if;
            Include(~seen, key);
            ok, y := IsRationalSquare(Evaluate(f, r));
            if ok then
                Append(~pts, <r,y>);
                if y ne 0 then
                    Append(~pts, <r,-y>);
                end if;
            end if;
        end for;
    end for;
    return pts;
end function;

function FpPointCount(f, p)
    F := GF(p);
    PF<xx> := PolynomialRing(F);
    fp := PF!f;
    count := 0;
    for a in F do
        v := Evaluate(fp, a);
        if v eq 0 then
            count +:= 1;
        elif IsSquare(v) then
            count +:= 2;
        end if;
    end for;
    lc := F!LeadingCoefficient(f);
    if IsSquare(lc) then
        count +:= 2;
    end if;
    return count;
end function;

out := Open(output_file, "w");
fprintf out, "# Rational-point diagnostics for extra-2 genus-3 curves\n\n";

for item in curves do
    label := item[1];
    f := item[2];
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    fprintf out, "============================================================\n";
    fprintf out, "%o\n", label;
    fprintf out, "f = %o\n", f;
    fprintf out, "factorization = %o\n", Factorization(f);
    disc := Discriminant(f);
    disc_int := Z!Abs(Numerator(disc));
    fprintf out, "degree = %o genus = %o discriminant = %o\n",
            Degree(f), Genus(C), disc;
    fprintf out, "bad_prime_factors = %o\n", PrimeDivisors(disc_int);
    fprintf out, "leading_coeff = %o square_over_Q = %o\n",
            LeadingCoefficient(f), IsSquare(LeadingCoefficient(f));

    fprintf out, "small prime point counts:\n";
    for p in [7,11,13,17,19,23,29,31,37,41,43,47] do
        if disc_int mod p eq 0 then
            fprintf out, "  p=%o bad\n", p;
        else
            fprintf out, "  p=%o count=%o\n", p, FpPointCount(f, p);
        end if;
    end for;

    for H in [20,100,500,2000] do
        pts := RationalSearch(f, H);
        fprintf out, "search_height %o affine_points %o\n", H, pts;
    end for;

    try
        rb := RankBounds(J);
        fprintf out, "RankBounds(J) = %o\n", rb;
    catch e
        fprintf out, "RankBounds failed: %o\n", e`Object;
    end try;

    try
        T, phi := TorsionSubgroup(J);
        fprintf out, "J torsion invariants = %o\n", Invariants(T);
    catch e
        fprintf out, "TorsionSubgroup failed: %o\n", e`Object;
    end try;
end for;

delete out;
print "Wrote", output_file;
quit;
