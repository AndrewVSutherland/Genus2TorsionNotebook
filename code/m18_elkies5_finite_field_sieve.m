//////////////////////////////////////////////////////////////////////
//  Finite-field sieve for combining Elkies 5-torsion with the
//  M_1(8) family
//
//      C: y^2 = q(x) * (x^4 + q(x)),
//      q = a*x^2 + b*x + c.
//
//  A rational specialization with full rational 2-torsion and a rational
//  5-torsion point has good reductions for which the sextic is split and
//  5 divides #J(F_p).  If no good split residue mod p has 5 | #J(F_p),
//  then every rational example must reduce to the boundary at p.
//
//  Typical run from torsion_jac:
//      magma code/m18_elkies5_finite_field_sieve.m
//////////////////////////////////////////////////////////////////////

prime_list := [3,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61];

function ParameterTriples(F)
    pts := [];
    for a in F do
        for b in F do
            for c in F do
                if a eq 0 and b eq 0 and c eq 0 then
                    continue;
                end if;
                Append(~pts, [a,b,c]);
            end for;
        end for;
    end for;
    return pts;
end function;

function CurvePolynomial(F, vals)
    P<x> := PolynomialRing(F);
    a := vals[1];
    b := vals[2];
    c := vals[3];
    q := a*x^2 + b*x + c;
    f := q * (x^4 + q);
    return f;
end function;

function IsGoodFullSplit(f)
    if Degree(f) ne 6 or Discriminant(f) eq 0 then
        return false;
    end if;
    roots := Roots(f);
    return #roots eq 6 and &and [ rt[2] eq 1 : rt in roots ];
end function;

for p in prime_list do
    F := GF(p);
    pts := ParameterTriples(F);
    good := 0;
    full_split := 0;
    five_possible := [];
    examples_split := [];

    for vals in pts do
        f := CurvePolynomial(F, vals);
        if Degree(f) eq 6 and Discriminant(f) ne 0 then
            good +:= 1;
        end if;
        if not IsGoodFullSplit(f) then
            continue;
        end if;

        full_split +:= 1;
        if #examples_split lt 5 then
            Append(~examples_split, vals);
        end if;

        C := HyperellipticCurve(f);
        J := Jacobian(C);
        if (#J mod 5) eq 0 then
            Append(~five_possible, <vals, #J>);
            print "  candidate", vals, "#J", #J;
        end if;
    end for;

    print "p", p, "triples", #pts, "good", good, "full_split", full_split,
          "five_possible", #five_possible;
    if #examples_split gt 0 then
        print "  split examples", examples_split;
    end if;
    if full_split gt 0 and #five_possible eq 0 then
        print "  NO GOOD FULL-SPLIT M1(8)+5 RESIDUES MOD", p;
        break;
    end if;
end for;

quit;
