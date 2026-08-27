//////////////////////////////////////////////////////////////////////
//  Special fibers at the finite p=3 boundary residues of the
//  contact-7 family.
//
//  This records the factorization of the mod-3 special fiber and, when
//  the fiber has one double rational root and squarefree residual cubic,
//  the order of the generalized Jacobian of the nodal limit:
//
//      #J_gen(F_3) = #Jac(y^2=residual_cubic)(F_3) * #torus(F_3).
//////////////////////////////////////////////////////////////////////

SetColumns(0);

F := GF(3);
P<x> := PolynomialRing(F);

function Fbar(a,b)
    return x^5 + (b^2 - F!7)*x^4 + (2*a*b + F!21)*x^3
        + (a^2 - F!7*b - F!35)*x^2
        + (-F!7*a + 2*b + F!35)*x
        + (2*a - F!35/F!4);
end function;

classes := [ <0,1>, <1,0>, <1,1>, <2,2> ];

for ab in classes do
    a := F!ab[1];
    b := F!ab[2];
    f := Fbar(a,b);
    print "CLASS", ab, "f", f, "factor", Factorization(f);
    print " squarefree_factorization", SquarefreeFactorization(f);

    fac := Factorization(f);
    for item in fac do
        if Degree(item[1]) eq 1 and item[2] ge 2 then
            root := -Coefficient(item[1],0)/Coefficient(item[1],1);
            g := ExactQuotient(f, item[1]^2);
            print " double_root", root, "residual", g,
                  "residual_factor", Factorization(g), "g(root)", Evaluate(g,root);

            if Degree(g) in {3,4} and Discriminant(g) ne 0 then
                C := HyperellipticCurve(g);
                Ecount := #Jacobian(C);
                split := IsSquare(Evaluate(g,root));
                torus_order := split select (#F - 1) else (#F + 1);
                print " normalization_J_count", Ecount,
                      "node_split", split,
                      "torus_order", torus_order,
                      "generalized_order", Ecount*torus_order,
                      "mod5", (Ecount*torus_order) mod 5;
            end if;
        end if;
    end for;
end for;

quit;
