//////////////////////////////////////////////////////////////////////
//  Richelot-neighbor test for the target [6,12].
//
//  Each seed below has exact rational torsion [6,6] and factor type
//  [1,2,2].  Its rational factorization supplies a pointwise-rational
//  maximal isotropic subgroup of J[2].  A degree-4 Richelot isogeny
//  preserves the full rational 3-primary subgroup (Z/3)^2.  Thus a
//  Richelot codomain whose 2-primary torsion is [2,4] has torsion [6,12].
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

function IntegralSquareScale(f)
    d := LCM([Denominator(Coefficient(f, i)) : i in [0..Degree(f)]]);
    return P!(d^2*f), d;
end function;

function TorsionData(C)
    f, h := HyperellipticPolynomials(C);
    if Degree(h) ge 0 then
        error "expected y^2=f model";
    end if;
    fI, d := IntegralSquareScale(P!f);
    CI := HyperellipticCurve(fI);
    JI := Jacobian(CI);
    A, mp := TorsionSubgroup(JI);
    return Invariants(A), fI;
end function;

function IsTarget612(inv)
    return inv eq [6,12];
end function;

function InvOrder(inv)
    return #inv eq 0 select 1 else &*inv;
end function;

function InvExponent(inv)
    return #inv eq 0 select 1 else inv[#inv];
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
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

seeds := [
    <"simple_h14_a133_39_bm7_13",
     11389248*x^5 - 18252000*x^4 + 42399396*x^3
       - 10288044*x^2 + 29659500*x>,
    <"core_a_m19_9_b3_2",
     944784*x^5 - 1781676*x^4 + 1644624*x^3
       - 791856*x^2 + 186624*x>,
    <"core_a_m43_25_b1_8",
     10000000000*x^5 - 29479000000*x^4 + 34512000000*x^3
       - 18866560000*x^2 + 4096000000*x>,
    <"core_a_m15_8_b5_9",
     191102976*x^5 - 495590400*x^4 + 535237632*x^3
       - 278769600*x^2 + 60466176*x>
];

print "CONTACT6_M612_RICHELOT_SEED_TEST";
print "seed_count", #seeds;

hits := 0;
for rec in seeds do
    label := rec[1];
    f := rec[2];
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    invS, fS := TorsionData(C);
    simple, pcert, Lp := SimpleCertificate(f);
    print "SOURCE", label, "torsion", invS,
          "factor_degrees", Sort([Degree(ff[1]) : ff in Factorization(f)]),
          "simple_cert", simple, "pcert", pcert;
    if simple then
        print "SOURCE_LPOLY", Lp;
    end if;

    try
        Rs := RichelotIsogenousSurfaces(J);
        print "RICHELOT_COUNT", label, #Rs;
        for i in [1..#Rs] do
            obj := Rs[i];
            print "RICHELOT_TYPE", label, i, Type(obj);
            if Type(obj) ne JacHyp then
                print "RICHELOT_NONJAC", label, i, obj;
                continue;
            end if;
            invR, fR := TorsionData(Curve(obj));
            print "RICHELOT", label, i, "torsion", invR,
                  "order", InvOrder(invR), "exponent", InvExponent(invR),
                  "factor_degrees", Sort([Degree(ff[1]) : ff in Factorization(fR)]);
            print "RICHELOT_CURVE", label, i, fR;
            if IsTarget612(invR) then
                hits +:= 1;
                print "HIT_6_12", label, i, fR;
            end if;
        end for;
    catch e
        print "RICHELOT_ERROR", label, e`Object;
    end try;

    if simple then
        try
            Js, products, weil_restrictions := TwoPowerIsogenies(J);
            print "TWOPOWER_COUNT", label, #Js, #products, #weil_restrictions;
            for i in [1..#Js] do
                invT, fT := TorsionData(Curve(Js[i]));
                print "TWOPOWER", label, i, "torsion", invT,
                      "order", InvOrder(invT), "exponent", InvExponent(invT);
                print "TWOPOWER_CURVE", label, i, fT;
                if IsTarget612(invT) then
                    hits +:= 1;
                    print "HIT_6_12_TWOPOWER", label, i, fT;
                end if;
            end for;
        catch e
            print "TWOPOWER_ERROR", label, e`Object;
        end try;
    end if;
end for;

print "DONE hits", hits;
quit;
