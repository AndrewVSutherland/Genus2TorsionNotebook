//////////////////////////////////////////////////////////////////////
// Projective parameter masks for the three explicit rational curves on
// the full A(2,2,2,8) cover used in target_22224_a2228_curves_plus3.m.
//
// A rational 3-torsion class injects into J(F_p) at every prime p != 3
// of good reduction.  Thus a projective parameter [n:d] is retained if
// the specialized branch tuple is singular modulo p, or if
//             3 | #J_[n:d](F_p).
//
// The homogeneous branch tuples below remove all denominator ambiguity,
// including [1:0].  Output files are suitable for the CRT height sieve.
//
// Run from torsion_jac, for example:
//
//   magma -b PrimeList:=5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61 \
//     code/target_22224_a2228_curves_plus3_masks.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned PrimeList then
    PrimeList := [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61];
elif Type(PrimeList) eq MonStgElt then
    PrimeList := [StringToInteger(s) : s in Split(PrimeList,",")];
end if;

if not assigned output_prefix then
    output_prefix := "results/target_22224_a2228_curves_plus3_masks";
end if;
if not assigned FamilyList then
    FamilyList := [1,2,3,4,5,6];
elif Type(FamilyList) eq MonStgElt then
    FamilyList := [StringToInteger(s) : s in Split(FamilyList,",")];
end if;

Z := Integers();

function HomogeneousTuple(family,n,d)
    if family eq 1 then
        // Common degree four.  On d != 0 this is d^4 times Filip1(n/d).
        return [
            -(n^2+n*d+d^2)^2,
            -4*n*(n+d)^2*d,
            4*n*(n+d)*d^2,
            4*n^2*(n+d)*d
        ];
    elif family eq 2 then
        // Multiply Filip2(n/d) by n*d*D, where
        // D=n^4+2n^3d-n^2d^2-2nd^3+d^4.
        N := n^4-2*n^3*d-n^2*d^2+2*n*d^3+d^4;
        D := n^4+2*n^3*d-n^2*d^2-2*n*d^3+d^4;
        return [-n*d*N,-d^2*D,n*d*D,n^2*D];
    elif family eq 3 then
        // Common degree four for Adam's recovered full-cover curve.
        return [
            -n*(n+d)^2*d,
            (n^2-d^2)*d^2,
            -n^2*(n^2-d^2),
            n*(n-d)^2*d
        ];
    elif family eq 4 then
        // FilipTupleT3.  Multiply the affine tuple by d^2*D, where
        // D=(n^2-2nd-d^2)(n^2+d^2).
        N := n*(n+d)^2*(n-d);
        D := (n^2-2*n*d-d^2)*(n^2+d^2);
        return [-N*d^2,-n^2*D,d^2*D,n*d*D];
    elif family in {5,6} then
        // Product-surface fibers (1,r,s,rs), with
        // s=(B-A*r)/(A-B*r).  Families 5 and 6 correspond to the
        // auxiliary parameters 1/2 and 2/3, hence (A,B)=(16,9),(144,25).
        A := family eq 5 select 16 else 144;
        B := family eq 5 select 9 else 25;
        DD := A*d-B*n;
        NN := B*d-A*n;
        return [d*DD,n*DD,d*NN,n*NN];
    end if;
    error "family must lie in [1..6]";
end function;

function CurvePolynomial(vals)
    F := Universe(vals);
    P<X> := PolynomialRing(F);
    return X*&*[X+z^2:z in vals];
end function;

function SmoothTuple(vals)
    sq := [z^2:z in vals];
    return &and[z ne 0:z in sq] and #Set(sq) eq 4;
end function;

function ThreeDividesJacobian(vals)
    f := CurvePolynomial(vals);
    C := HyperellipticCurve(f);
    nJ := Z!Evaluate(LPolynomial(C),1);
    return nJ mod 3 eq 0,nJ;
end function;

for p in PrimeList do
    assert IsPrime(p) and p notin {2,3};
    F := GF(p);
    filename := Sprintf("%o_p%o.tsv",output_prefix,p);
    out := Open(filename,"w");
    fprintf out,"family\tn\td\tparameter\tstatus\tjacobian_order\n";
    for family in FamilyList do
        boundary := []; good_allowed := []; good_killed := [];
        pts := [<r,F!1>:r in F] cat [<F!1,F!0>];
        for pt in pts do
            n,d := Explode(pt);
            vals := HomogeneousTuple(family,n,d);
            label := d eq 0 select "infinity" else Sprint(Z!n);
            if not SmoothTuple(vals) then
                Append(~boundary,<Z!n,Z!d>);
                fprintf out,"%o\t%o\t%o\t%o\tboundary\t0\n",
                    family,Z!n,Z!d,label;
                continue;
            end if;
            has3,nJ := ThreeDividesJacobian(vals);
            if has3 then
                Append(~good_allowed,<Z!n,Z!d>);
                fprintf out,"%o\t%o\t%o\t%o\tgood_allowed\t%o\n",
                    family,Z!n,Z!d,label,nJ;
            else
                Append(~good_killed,<Z!n,Z!d>);
            end if;
        end for;
        print "P1_MASK","p",p,"family",family,
              "boundary",boundary,"good_allowed",good_allowed,
              "good_killed_count",#good_killed,
              "total_allowed",#boundary+#good_allowed;
    end for;
    delete out;
    print "P1_MASK_FILE",filename;
end for;

quit;
