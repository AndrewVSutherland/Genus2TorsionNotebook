//////////////////////////////////////////////////////////////////////
// Exact class-by-class test on the distinguished contact-6 Richelot dual.
//
// For a contact-6 source
//
//   y^2 = x * B * C,
//
// the distinguished dual is (up to the standard Richelot twist sign)
//
//   y^2 = Delta * R1 * R2 * R3.
//
// This file tests divisibility by 2 of the three rational dual-kernel
// classes [Ri,0] directly.  The resultant-square tests
//
//   R1 => DB*DC square, R2 => DC square, R3 => DB square
//
// are necessary filters only; the calls below are exact.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
P<x> := PolynomialRing(Q);

if not assigned run_known_seed then run_known_seed := true;
elif Type(run_known_seed) eq MonStgElt then
    run_known_seed := run_known_seed in {"true","True","1","yes"};
end if;

function IntegralSquareScale(f)
    d := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
    return P!(d^2*f), d;
end function;

function HasTwoIndependent3Directions(inv)
    return #[n : n in inv | (Integers()!n) mod 3 eq 0] ge 2;
end function;

function DualData(a,b)
    Delta := a*b+3*a+3*b+5;
    R1 := (b^2-2*a-3)*x^2
          +(2*a*b+6*a+6*b+10)*x +(a^2-2*b-3);
    R2 := 2*x^2-(a+3);
    R3 := 2-(b+3)*x^2;
    DB := (a-3)^2-8*(b+3);
    DC := (b-3)^2-8*(a+3);
    return Delta,R1,R2,R3,DB,DC;
end function;

function ExactDualClassRecord(a,b,twist_sign)
    Delta,R1,R2,R3,DB,DC := DualData(a,b);
    if Delta eq 0 then
        return false,<twist_sign,[],[],false,false,false,P!0>;
    end if;
    g := twist_sign*Delta*R1*R2*R3;
    if Degree(g) notin {5,6} or Discriminant(g) eq 0 then
        return false,<twist_sign,[],[],false,false,false,P!0>;
    end if;
    gI,d := IntegralSquareScale(g);
    J := Jacobian(HyperellipticCurve(gI));
    try
        G,mp := TorsionSubgroup(J);
        inv := Invariants(G);
        classes := [J![r/LeadingCoefficient(r),0] : r in [R1,R2,R3]];
        halves := [IsDivisibleBy(T,2) : T in classes];
        exact612 := inv eq [6,12];
        return true,<twist_sign,inv,halves,
                     HasTwoIndependent3Directions(inv),exact612,
                     &or halves,gI>;
    catch e
        return false,<twist_sign,[],[],false,false,false,gI>;
    end try;
end function;

procedure PrintDualClassAudit(label,a,b)
    Delta,R1,R2,R3,DB,DC := DualData(a,b);
    sqB := DB ne 0 and IsSquare(DB);
    sqC := DC ne 0 and IsSquare(DC);
    sqMix := DB*DC ne 0 and IsSquare(DB*DC);
    print "DUAL_CLASS_SOURCE",label,"a",a,"b",b,"Delta",Delta;
    print "NECESSARY_SQUARES","DB",DB,sqB,"DC",DC,sqC,
          "DBDC",DB*DC,sqMix;
    for eps in [1,-1] do
        ok,rec := ExactDualClassRecord(a,b,eps);
        print "DUAL_TWIST","sign",eps,"ok",ok,
              "torsion",rec[2],"halves_R1_R2_R3",rec[3],
              "has_33",rec[4],"exact_612",rec[5],
              "any_half",rec[6];
    end for;
end procedure;

if run_known_seed then
    print "CONTACT6_M612_DUAL_CLASS_EXACT";
    PrintDualClassAudit("known_simple_66",Q!(133/39),Q!(-7/13));
end if;

if not assigned NoMain then quit; end if;
