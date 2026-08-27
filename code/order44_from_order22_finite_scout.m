//////////////////////////////////////////////////////////////////////
// Bounded, class-specific finite-field scout for cyclic order 44.
//
// The order-22 families have a marked class
//
//     E = W - infinity_+,       2E = -(infinity_+-infinity_-),
//
// where W=(r,0).  Move W to infinity by
//
//     X=1/(x-r),  Y=y*X^3,  g(X)=X^6*f(r+1/X).
//
// Then infinity_+ maps to (0,1), and E maps to the negative of
// D=(0,1)-infinity.  Thus E is divisible by 2 iff D is.  On the odd
// quintic model Magma represents D as J![X,1].
//
// This script enumerates only the honest affine open s != 0,+/-1.
// The four base points s=0,+/-1,infinity, as well as singular open
// reductions, are explicitly KEPT as conservative boundary/bad-
// reduction cases; they are never used to claim a local obstruction.
//
// Typical bounded run:
//   magma -b primes:="3,5,7,13,17,19,23" \
//       code/order44_from_order22_finite_scout.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Z := Integers();

if assigned primes then
    if Type(primes) eq MonStgElt then
        prime_list := [StringToInteger(a) : a in Split(primes, ",") | #a gt 0];
    else
        prime_list := primes;
    end if;
else
    // p=11 is omitted because the marked order is divisible by 11.
    prime_list := [3,5,7,13,17,19,23,29,31];
end if;

if assigned max_print then
    if Type(max_print) eq MonStgElt then
        max_print := StringToInteger(max_print);
    end if;
else
    max_print := 8;
end if;

function IsDivisibleBy2Finite(D,G,phi)
    a := D @@ phi;
    coords := Eltseq(a);
    invs := Invariants(G);
    for i in [1..#coords] do
        if GCD(2,invs[i]) ne 1 and coords[i] mod 2 ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function FlynnPolynomial(F,P,t)
    x := P.1;
    return x^6 + 2*x^5 + (2*t+3)*x^4 + 2*x^3
           + (t^2+1)*x^2 + 2*t*(1-t)*x + t^2;
end function;

function DSPolynomial(F,P,u)
    x := P.1;
    return x^6 - 4*x^5 + 8*(1+u)*x^4 - (10+32*u)*x^3
           + 8*(1+6*u+2*u^2)*x^2
           - 4*(1+6*u+16*u^2)*x + 64*u^2+1;
end function;

function ParameterNumerator(s,eps)
    return -s^2*(s^2+1)*(s^4-s^2+1) + 2*eps*s^5;
end function;

function MoveMarkedRootToInfinity(f,r)
    P := Parent(f);
    X := P.1;
    // X^6*f(r+1/X), written without Laurent polynomials.
    g := &+[Coefficient(f,i)*(r*X+1)^i*X^(6-i) : i in [0..6]];
    return P!g;
end function;

function TestFamilyPrime(label,p)
    F := GF(p);
    P<X> := PolynomialRing(F);

    open_pairs := 0;
    good_models := 0;
    singular_open := 0;
    marked22 := 0;
    collapsed_order := 0;
    divisible := 0;
    examples := [];

    // P^1(F_p) has p+1 points.  The parameterization's displayed open
    // removes 0,+/-1 and infinity.  We keep all four base boundary
    // points conservatively, without pretending that the eps sheets
    // have a regular model there.
    boundary_base_points_kept := 4;

    for s in F do
        if s eq 0 or s^2 eq 1 then
            continue;
        end if;
        den := (s^2-1)^2;
        for epsZ in [-1,1] do
            eps := F!epsZ;
            open_pairs +:= 1;
            num := ParameterNumerator(s,eps);
            if label eq "Flynn" then
                par := num/den;
                r := s^2;
                f := FlynnPolynomial(F,P,par);
            else
                par := num/(4*den);
                r := 1+s^2;
                f := DSPolynomial(F,P,par);
            end if;

            if Degree(f) ne 6 or Evaluate(f,r) ne 0 or Discriminant(f) eq 0 then
                singular_open +:= 1;
                continue;
            end if;
            g := MoveMarkedRootToInfinity(f,r);
            if Degree(g) ne 5 or Coefficient(g,0) ne 1 or Discriminant(g) eq 0 then
                singular_open +:= 1;
                continue;
            end if;

            good_models +:= 1;
            C := HyperellipticCurve(g);
            J := Jacobian(C);
            D := J![X,F!1];
            ord := Order(D);
            if ord ne 22 then
                collapsed_order +:= 1;
                continue;
            end if;
            marked22 +:= 1;
            G,phi := AbelianGroup(J);
            if IsDivisibleBy2Finite(D,G,phi) then
                divisible +:= 1;
                if #examples lt max_print then
                    Append(~examples,
                        <Z!s,epsZ,Z!par,Invariants(G),Eltseq(D@@phi)>);
                end if;
            end if;
        end for;
    end for;

    print "family",label,"p",p,
          "P1_base_points",p+1,
          "boundary_base_points_kept",boundary_base_points_kept,
          "open_sheet_pairs",open_pairs,
          "good_models",good_models,
          "singular_open_kept",singular_open,
          "marked_order22",marked22,
          "collapsed_order_kept",collapsed_order,
          "marked_divisible_by_2",divisible;
    if #examples gt 0 then
        print "  examples <s,eps,source_parameter,J_invariants,D_coordinates>",examples;
    elif marked22 gt 0 then
        print "  NO OPEN MARKED HALVES; boundary and bad reductions remain KEPT";
    else
        print "  NO ORDER-22 OPEN TESTS; this prime gives no scout information";
    end if;
    return <open_pairs,good_models,singular_open,marked22,collapsed_order,divisible>;
end function;

print "ORDER44_FROM_ORDER22_FINITE_SCOUT";
print "criterion: exact marked D=(0,1)-infinity is in 2*J(F_p)";
print "projective policy: s=0,+/-1,infinity and singular reductions are KEPT";

for p in prime_list do
    require IsPrime(p) and p ne 2 and p ne 11:
        "primes must be odd and different from 11";
    _ := TestFamilyPrime("Flynn",p);
    _ := TestFamilyPrime("DaowsudSchmidt",p);
end for;

quit;
