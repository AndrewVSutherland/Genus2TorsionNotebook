//////////////////////////////////////////////////////////////////////
// Direct rational-halving test on every recorded exact [6,6] contact-6
// curve.  We test all three nonzero rational 2-classes on the factor-type
// [1,2,2] model
//
//   f = x * ((b+3)x^2+(a-3)x+2)
//         * (2x^2+(b-3)x+(a+3)).
//
// A rational half of any one of these classes enlarges the 2-primary
// subgroup from [2,2] to [2,4].  Since the source has full rational
// 3-primary [3,3], this would force torsion [6,12].
//
// Run from torsion_jac:
//   magma code/contact6_m612_direct_halves_known_hits.m \
//     > data/contact6_m612_direct_halves_known_hits.txt
//////////////////////////////////////////////////////////////////////

SetColumns(0);
load "../../halving_mumford_tools.m";

Q := Rationals();
P<x> := PolynomialRing(Q);

function IntegralSquareScale(f)
    d := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
    return P!(d^2*f), d;
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67] do
        fp := ChangeRing(f,GF(p));
        if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
            continue;
        end if;
        Lp := LPolynomial(ChangeRing(C,GF(p)));
        fac := Factorization(Lp);
        if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
            return true,p,Lp;
        end if;
    end for;
    return false,0,P!0;
end function;

seeds := [
    <"simple_h14", Q!133/39, Q!-7/13>,
    <"core_h6", Q!-19/9, Q!3/2>,
    <"core_h10_2", Q!-43/25, Q!1/8>,
    <"core_h10_3", Q!-15/8, Q!5/9>
];

print "CONTACT6_M612_DIRECT_HALVES_KNOWN_HITS";
print "seed_count",#seeds;

totalHalves := 0;
verifiedHalves := 0;
targetHits := 0;

for seed in seeds do
    label,a,b := Explode(seed);
    qA := x;
    qB := x^2+((a-3)/(b+3))*x+2/(b+3);
    qC := x^2+((b-3)/2)*x+(a+3)/2;
    f := qA*((b+3)*x^2+(a-3)*x+2)
          *(2*x^2+(b-3)*x+(a+3));

    assert Degree(f) eq 5 and Discriminant(f) ne 0;
    CI := HyperellipticCurve(f);
    J := Jacobian(CI);
    h6 := 1+a*x+b*x^2+x^3;
    D6 := J![x-1,Evaluate(h6,Q!1)];
    assert Order(D6) eq 6;
    TG,mp := TorsionSubgroup(J);
    inv := Invariants(TG);
    assert inv eq [6,6];

    fi,scale := IntegralSquareScale(f);
    simple,pcert,Lp := SimpleCertificate(f);
    print "SOURCE",label,"a",a,"b",b,"torsion",inv,
          "simple_cert",simple,"pcert",pcert;
    print "SOURCE_FACTORS",Factorization(f);
    print "SOURCE_INTEGRAL",fi,"scale",scale;
    if simple then print "SOURCE_LPOLY",Lp; end if;

    twoPolys := [<"x_vs_infinity",qA>,
                 <"quadratic_B",qB>,
                 <"quadratic_C",qC>];
    twoClasses := [];
    for rec in twoPolys do
        name := rec[1];
        u := rec[2];
        T := J![u,P!0];
        assert T ne J!0 and 2*T eq J!0;
        Append(~twoClasses,T);

        halves := HalvesFromMumford(u,P!0,f);
        print "TWO_CLASS",label,name,"u",u,"halves",#halves;
        totalHalves +:= #halves;
        for i in [1..#halves] do
            h := halves[i];
            H := DivisorFromHalfTuple(J,h);
            ok := 2*H eq T;
            ord := PointOrderBounded(H,24);
            print " HALF",i,"M",h[1],"N",h[2],"G",h[3],
                  "verified",ok,"order",ord;
            if ok then
                verifiedHalves +:= 1;
                assert ord eq 4;
                TG2,mp2 := TorsionSubgroup(J);
                inv2 := Invariants(TG2);
                print " VERIFIED_TARGET_TORSION",inv2;
                if inv2 eq [6,12] then
                    targetHits +:= 1;
                    print "HIT_6_12",label,name,fi;
                end if;
            end if;
        end for;
    end for;
    assert twoClasses[1] ne twoClasses[2]
           and twoClasses[1] ne twoClasses[3]
           and twoClasses[2] ne twoClasses[3];
    assert twoClasses[1]+twoClasses[2]+twoClasses[3] eq J!0;
    assert 3*D6 eq twoClasses[2];
    print "MARKED_RELATION",label,"3D_equals_quadratic_B",true;
    print "END_SOURCE",label;
end for;

print "DONE","total_halves",totalHalves,
      "verified_halves",verifiedHalves,"target_hits",targetHits;
quit;
