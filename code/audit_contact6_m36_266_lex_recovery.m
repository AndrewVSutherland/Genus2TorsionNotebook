//////////////////////////////////////////////////////////////////////
// Exact recovery audit using the rational lex basis modulo each plane
// factor.  Unlike a fresh function-field Groebner basis, this only reduces
// the already-computed lex relations in Q(v)[U]/(factor).
//
// Run under an external timeout:
//   magma -b code/audit_contact6_m36_266_lex_recovery.m
//////////////////////////////////////////////////////////////////////

load "code/audit_contact6_m36_266_plane_extract.m";

Fv<t> := FunctionField(Q);
Pu<z> := PolynomialRing(Fv);
print "LEX_RECOVERY_AUDIT_START";

for item in fac do
    plane_factor := item[1];
    d := TotalDegree(plane_factor);
    gu := Pu!Evaluate(plane_factor, <Q!0,Q!0,z,t>);
    K<alpha> := ext<Fv | gu>;
    S<beta,mu> := PolynomialRing(K,2,"lex");
    rho := hom<L -> S | beta,mu,S!(K!alpha),S!(K!t)>;
    reduced_relations := [rho(q) : q in GL | rho(q) ne 0];
    print "LEX_FACTOR_BEGIN", d;
    print "nonzero_relation_count", #reduced_relations;
    print "relation_summary",
          [<Degree(q,1),Degree(q,2),TotalDegree(q),#Terms(q)>
           : q in reduced_relations];
    linear := [q : q in reduced_relations | TotalDegree(q) le 1];
    print "linear_relation_count", #linear;
    for q in linear do
        print "linear_relation", q;
    end for;
    print "LEX_FACTOR_END", d;
end for;

print "LEX_RECOVERY_AUDIT_DONE";
