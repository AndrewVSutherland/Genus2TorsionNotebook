//////////////////////////////////////////////////////////////////////
// Generic-fiber audit for the projections of the exact r=4 components
// to their (U,v) plane images.  This deliberately works over
// Q(v)[U]/(plane_factor), so the resulting zero-dimensional degree is the
// function-field degree of the projection after structural saturation.
//
// Run under an external timeout:
//   magma -b code/audit_contact6_m36_266_generic_fiber.m
//////////////////////////////////////////////////////////////////////

load "code/audit_contact6_m36_266_plane_extract.m";

Fv<t> := FunctionField(Q);
Pu<z> := PolynomialRing(Fv);

print "GENERIC_FIBER_AUDIT_START";
for item in fac do
    plane_factor := item[1];
    d := TotalDegree(plane_factor);
    gu := Pu!Evaluate(plane_factor, <Q!0,Q!0,z,t>);
    print "GENERIC_FACTOR_BEGIN", d;
    print "factor_degree_over_Qv", Degree(gu),
          "irreducible_over_Qv", IsIrreducible(gu);
    K<alpha> := ext<Fv | gu>;
    S<beta,mu> := PolynomialRing(K,2,"lex");
    rho := hom<R -> S | beta,mu,S!(K!alpha),S!(K!t)>;
    Jraw := ideal<S | rho(F1),rho(F2),rho(F3)>;
    J := Saturation(Jraw, ideal<S | rho(structural)>);
    GJ := GroebnerBasis(J);
    dim, component_degrees := Dimension(J);
    print "generic_fiber_dimension", dim,
          "component_degrees", component_degrees;
    print "generic_fiber_basis_length", #GJ;
    print "generic_fiber_basis_summary",
          [<Degree(q,1),Degree(q,2),TotalDegree(q),#Terms(q)> : q in GJ];
    for q in GJ do
        if TotalDegree(q) le 1 then
            print "linear_recovery", q;
        end if;
    end for;
    print "GENERIC_FACTOR_END", d;
end for;
print "GENERIC_FIBER_AUDIT_DONE";
