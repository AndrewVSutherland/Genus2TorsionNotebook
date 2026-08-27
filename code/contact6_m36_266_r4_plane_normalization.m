//////////////////////////////////////////////////////////////////////
// Plane normalization diagnostics for the exact [2,6,6] r=4 slice.
//////////////////////////////////////////////////////////////////////
SetColumns(0);
if not assigned analysis_mode then analysis_mode := "all"; end if;
if not assigned factor_index then
    factor_index := 0;
elif Type(factor_index) eq MonStgElt then
    factor_index := StringToInteger(factor_index);
end if;
p := 0;
mode := "lex";
load "code/contact6_m36_266_r4_magma_decompose.m";
plane_candidates := [g : g in GL | Degree(g,1) eq 0 and Degree(g,2) eq 0];
require #plane_candidates eq 1: "expected exactly one elimination polynomial in U,v";
plane_factors := Factorization(plane_candidates[1]);
require #plane_factors eq 2 and &and[t[2] eq 1 : t in plane_factors]: "expected two squarefree irreducible plane factors";
print "PLANE_NORMALIZATION_START";
print "factor_degrees", [TotalDegree(t[1]) : t in plane_factors];
A2<u,w> := AffineSpace(K,2);
indices := factor_index eq 0 select [1..#plane_factors] else [factor_index];
for i in indices do
    plane_factor := plane_factors[i][1];
    print "FACTOR", i, "total_degree", TotalDegree(plane_factor), "terms", #Terms(plane_factor);
    print "FACTOR_EXPLICIT", i, plane_factor;
    cp := Evaluate(plane_factor, <K!0,K!0,u,w>);
    Cplane := Curve(A2,cp);
    Ffrac := FunctionField(Cplane);
    if analysis_mode in {"all","fiber"} then
        S<bbb,MMM> := PolynomialRing(Ffrac,2,"grevlex");
        rho := hom<R -> S | bbb,MMM,Ffrac.1,Ffrac.2>;
        Iraw_fiber := ideal<S | [rho(g) : g in [F1,F2,F3]]>;
        Ifiber := Saturation(
            Iraw_fiber, ideal<S | rho(structural*smooth*coprime)>
        );
        Ifiber := Saturation(
            Ifiber, ideal<S | [rho(g) : g in AutoA]>
        );
        Ifiber := Saturation(
            Ifiber, ideal<S | [rho(g) : g in AutoB]>
        );
        print "GENERIC_FIBER", i, "dimension", Dimension(Ifiber), "quotient_dimension", QuotientDimension(Ifiber), "basis_length", #Basis(Ifiber), "basis_summary", [<Degree(g,1),Degree(g,2),#Terms(g)> : g in Basis(Ifiber)];
        print "GENERIC_FIBER_BASIS", i, Basis(Ifiber);
    end if;
    if analysis_mode in {"all","genus"} then
        AF, to_AF := AlgorithmicFunctionField(Ffrac);
        print "NORMALIZATION_FUNCTION_FIELD", i, AF;
        print "NORMALIZATION_DEFINING_POLYNOMIAL", i, DefiningPolynomial(AF);
        print "NORMALIZATION_GENUS", i, Genus(AF);
        try
            PC := ProjectiveClosure(Cplane);
            print "PLANE_CURVE_GENUS_CROSSCHECK", i, Genus(PC);
        catch e
            print "PLANE_CURVE_GENUS_CROSSCHECK_ERROR", i, e`Object;
        end try;
    end if;
end for;
print "PLANE_NORMALIZATION_DONE";
