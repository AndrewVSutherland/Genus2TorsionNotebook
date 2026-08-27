//////////////////////////////////////////////////////////////////////
// Fast exact certificate for geometric simplicity and the absence of
// extra geometric endomorphisms of the genus-6 trigonal quotient.
//
// The two Frobenius polynomials are regenerated from the curve by
// contact30_c3root_genus6_frobenius.m at Prime=13 and Prime=17.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
P<x> := PolynomialRing(Q);

P13 := x^12+2*x^11+2*x^10+15*x^9+302*x^8+341*x^7+437*x^6
     +4433*x^5+51038*x^4+32955*x^3+57122*x^2+742586*x+4826809;
P17 := x^12+2*x^11-14*x^10-18*x^9+200*x^8+31*x^7-5037*x^6
     +527*x^5+57800*x^4-88434*x^3-1169294*x^2+2839714*x+24137569;

// The branch divisor of the trigonal map has no collision at 13 or 17.
ZI<z> := PolynomialRing(Integers());
P10 := 675105*z^10 - 24788818*z^9 + 404115156*z^8
     -3847306032*z^7 + 23658817344*z^6 - 98072909568*z^5
     +277203253248*z^4 - 526965620736*z^3
     +644275372032*z^2 - 457272459264*z + 143107555328;
Bbranch := (z-2)*(z-5)*(3*z-14)*(31*z^2-220*z+352)*P10;
branchdisc := Discriminant(Bbranch);
assert Set([ fe[1] : fe in Factorization(branchdisc) ])
    eq {2,3,5,7,11,269};
assert branchdisc mod 13 ne 0 and branchdisc mod 17 ne 0;
print "branch_discriminant_factorization",Factorization(branchdisc);

function RealWeilPolynomial(PP,q)
    rem := PP;
    cs := [];
    for d in [6..0 by -1] do
        basis := P!(x^6*(x+q/x)^d);
        c := Coefficient(rem,6+d);
        Append(~cs,c);
        rem -:= c*basis;
    end for;
    assert rem eq 0 and cs[1] eq 1;
    QY<y> := PolynomialRing(Q);
    return &+[ cs[i+1]*y^(6-i) : i in [0..6] ];
end function;

CMdiscs := [];
RealDiscs := [];
for pair in [ <13,P13>, <17,P17> ] do
    p := Integers()!pair[1];
    PP := pair[2];
    assert IsIrreducible(PP);
    assert (Integers()!Coefficient(PP,6)) mod p ne 0; // ordinary

    G,roots := GaloisGroup(PP);
    assert Order(G) eq 2^6*Factorial(6);
    idnum,iddeg := TransitiveGroupIdentification(G);
    assert idnum eq 293 and iddeg eq 12;

    // These finite checks are useful independently of the full-Weyl-group
    // argument, which proves irreducibility after every finite extension.
    assert &and[ IsIrreducible(WeilPolynomialOverFieldExtension(PP,n))
                 : n in [2..12] ];

    h := RealWeilPolynomial(PP,p);
    assert IsIrreducible(h);
    Kr<a> := NumberField(h);
    Kc<b> := NumberField(PP);
    dr := Discriminant(MaximalOrder(Kr));
    dc := Discriminant(MaximalOrder(Kc));
    Append(~RealDiscs,dr);
    Append(~CMdiscs,dc);

    // For the root field of the full hyperoctahedral group, the maximal
    // real field is its sole proper nontrivial subfield.  Verify this by
    // Magma's exact subfield computation as well.
    sf := Subfields(Kc);
    assert #sf eq 2 and Set([ Degree(tup[1]) : tup in sf ]) eq {6,12};

    print "prime",p;
    print "Frobenius",PP;
    print "ordinary",true;
    print "Galois_order",Order(G),"transitive_id",293,12;
    print "real_Weil_polynomial",h;
    print "CM_subfield_degrees",[6,12];
    print "real_field_discriminant",dr,Factorization(dr);
    print "CM_field_discriminant",dc,Factorization(dc);
end for;

assert RealDiscs[1] ne RealDiscs[2];
assert CMdiscs[1] ne CMdiscs[2];
N13 := Integers()!Evaluate(P13,1);
N17 := Integers()!Evaluate(P17,1);
assert N13 eq 5716043 and N17 eq 25773047 and GCD(N13,N17) eq 1;
print "real_fields_nonisomorphic_by_discriminant",true;
print "CM_fields_nonisomorphic_by_discriminant",true;
print "jacobian_orders",N13,N17,"gcd",GCD(N13,N17);
print "J_Q_torsion_trivial",true;
print "certificate_conclusion", "End_barQ(J)=Z and J is absolutely simple";

quit;
