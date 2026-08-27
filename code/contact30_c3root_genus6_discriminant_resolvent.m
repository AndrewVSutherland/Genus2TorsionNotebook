//////////////////////////////////////////////////////////////////////
// Hyperelliptic discriminant resolvent and S3-closure geometry for the
// genus-6 trigonal quotient of the contact-30 C3-root cover.
//
//   magma -b code/contact30_c3root_genus6_discriminant_resolvent.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(8*10^9);

if not assigned ComputeRank then
    ComputeRank := false;
elif Type(ComputeRank) eq MonStgElt then
    ComputeRank := StringToLower(ComputeRank) in
        { "true", "t", "yes", "y", "1" };
end if;

Q := Rationals();
Qz<z> := PolynomialRing(Q);

P10 :=
    675105*z^10 - 24788818*z^9 + 404115156*z^8
    - 3847306032*z^7 + 23658817344*z^6
    - 98072909568*z^5 + 277203253248*z^4
    - 526965620736*z^3 + 644275372032*z^2
    - 457272459264*z + 143107555328;
q2 := 31*z^2-220*z+352;
h := (z-2)*(z-5)*q2*P10;
assert IsSquarefree(h) and Degree(h) eq 14;
H := HyperellipticCurve(h);
assert Genus(H) eq 6;

print "CONTACT30 GENUS6 DISCRIMINANT RESOLVENT";
print "model y^2 =",h;
print "genus",Genus(H);
print "branch_factor_degrees",
    [ <Degree(fe[1]),fe[2]> : fe in Factorization(h) ];

// The S3 Galois closure is obtained by adjoining sqrt(discriminant) to
// the non-Galois cubic root curve.  Its ramification over the trigonal
// curve has degree 14, so Riemann--Hurwitz gives genus 18.
gD := 6;
ram_X_over_D := 14;
gX := 1+(2*(2*gD-2)+ram_X_over_D)/2;
assert gX eq 18;
print "S3_galois_closure_genus",gX;
print "jacobian_isogeny_dimensions_JX_JH_JD2",[18,6,6,6];

pts := Points(H : Bound := 10000);
print "point_search_bound",10000,"count",#pts;
print "point_search",pts;

if ComputeRank then
    JH := Jacobian(H);
    try
        time lo,hi := RankBounds(JH);
        print "jacobian_rank_bounds",lo,hi;
    catch e
        print "jacobian_rank_bounds_failed",e`Object;
    end try;
end if;

quit;
