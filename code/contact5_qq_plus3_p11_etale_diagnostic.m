//////////////////////////////////////////////////////////////////////
//  Local geometry of the remaining contact-5 [2,20]+3 branch.
//
//  The quadratic-quadratic extra-2 parameter is
//
//    t(r) = -(r^6-2r^5+2r^4-4r^3+4r^2-8r+8)
//             / ((r^2-2)^2 (r^2-2r+2)).
//
//  At p=11 all good fibers fail the necessary 3-divisibility test.
//  The only live residues are r=3,5,7,8, where t=-3 and the curve is
//  singular.  Earlier brute-force lifting found contact points modulo
//  11^2, but did not determine whether the contact cover is locally
//  open over the r-line.
//
//  This script enumerates the cubic-contact points over F_11 and
//  computes the Jacobian of the three residual contact equations.  In
//  particular, the 3x3 vertical block in (m,U,V) decides whether the
//  implicit-function theorem gives a contact point for every nearby r.
//  Derivatives are taken by substituting a+e over F_11(e).
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
Z := Integers();
PQ<rq> := PolynomialRing(Q);

numq := rq^6 - 2*rq^5 + 2*rq^4 - 4*rq^3 + 4*rq^2 - 8*rq + 8;
denq := (rq^2 - 2)^2*(rq^2 - 2*rq + 2);
boundary_num := 3*denq - numq;

print "CONTACT5_QQ_P11_ETALE_START";
print "boundary_numerator_factorization_Q", Factorization(boundary_num);

function TQQ(r)
    den := (r^2 - 2)^2*(r^2 - 2*r + 2);
    if den eq 0 then
        return false, Parent(r)!0;
    end if;
    num := r^6 - 2*r^5 + 2*r^4 - 4*r^3 + 4*r^2 - 8*r + 8;
    return true, -num/den;
end function;

function FamilyCoefficients(t)
    R := Parent(t);
    P<x> := PolynomialRing(R);
    b := (t^2 - 1)/2;
    h := 1 + t*x + b*x^2;
    f := h^2 - ((t + 1)^4/4)*x^5;
    return f, [ Coefficient(f,i) : i in [0..5] ];
end function;

// Unknown cubic and quadratic are
//   ell = m*x^3 + N*x^2 + Rr*x + S,
//   q   = x^2 + U*x + V.
// The top three coefficient equations determine N,Rr,S; these are the
// three remaining equations for ell^2-f=m^2*q^3.
function ContactEquations(r, m, U, V)
    ok, t := TQQ(r);
    if not ok or m eq 0 then
        error "contact equations evaluated at a pole";
    end if;
    f, coeffs := FamilyCoefficients(t);
    N := (3*m^2*U + coeffs[6])/(2*m);
    Rr := (3*m^2*(U^2 + V) + coeffs[5] - N^2)/(2*m);
    S := (m^2*(U^3 + 6*U*V) + coeffs[4] - 2*N*Rr)/(2*m);
    e2 := Rr^2 + 2*N*S - coeffs[3] - 3*m^2*(U^2*V + V^2);
    e1 := 2*Rr*S - coeffs[2] - 3*m^2*U*V^2;
    e0 := S^2 - coeffs[1] - m^2*V^3;
    return [e2,e1,e0];
end function;

function ContactSolutionsAtR(r, F)
    sols := [];
    for m in F do
        if m eq 0 then
            continue;
        end if;
        for U in F do
            for V in F do
                eqs := ContactEquations(r,m,U,V);
                if &and[e eq 0 : e in eqs] then
                    Append(~sols, <m,U,V>);
                end if;
            end for;
        end for;
    end for;
    return sols;
end function;

// Directional derivative at a finite-field point, computed over F(e).
function DirectionalColumn(base, direction, K, e)
    vals := [ K!base[i] + K!direction[i]*e : i in [1..4] ];
    eqs := ContactEquations(vals[1], vals[2], vals[3], vals[4]);
    return [ Evaluate(Derivative(eeq), K!0) : eeq in eqs ];
end function;

p := 11;
F := GF(p);
PF<rf> := PolynomialRing(F);
boundary_F := PF!boundary_num;
print "boundary_numerator_factorization_F11", Factorization(boundary_F);
print "boundary_gcd_derivative_F11", GCD(boundary_F, Derivative(boundary_F));

KF<e> := RationalFunctionField(F);
live := [3,5,7,8];
total_solutions := 0;
vertical_rank_counts := AssociativeArray();
total_rank_counts := AssociativeArray();

for a in live do
    r := F!a;
    ok, t := TQQ(r);
    assert ok and t eq -F!3;
    sols := ContactSolutionsAtR(r,F);
    total_solutions +:= #sols;
    print "RESIDUE", a, "t", t, "contact_solutions", #sols,
          "boundary_derivative", Evaluate(Derivative(boundary_F),r);

    residue_vranks := [];
    residue_tranks := [];
    for sol in sols do
        base := [r,sol[1],sol[2],sol[3]];
        cols := [];
        for j in [1..4] do
            dir := [ F!0 : i in [1..4] ];
            dir[j] := F!1;
            Append(~cols, DirectionalColumn(base,dir,KF,e));
        end for;
        // The construction above concatenates columns; transpose the
        // 4x3 row-major interpretation to obtain the intended 3x4 matrix.
        Jac := Transpose(Matrix(F,4,3,&cat[ [F!c : c in col] : col in cols ]));
        VJac := Submatrix(Jac,1,2,3,3);
        vr := Rank(VJac);
        tr := Rank(Jac);
        Append(~residue_vranks,vr);
        Append(~residue_tranks,tr);
        if IsDefined(vertical_rank_counts,vr) then
            vertical_rank_counts[vr] +:= 1;
        else
            vertical_rank_counts[vr] := 1;
        end if;
        if IsDefined(total_rank_counts,tr) then
            total_rank_counts[tr] +:= 1;
        else
            total_rank_counts[tr] := 1;
        end if;
        print "  SOL", sol, "vertical_rank", vr, "total_rank", tr,
              "r_column", [Jac[i,1] : i in [1..3]];
    end for;
    print "RESIDUE_RANK_SUMMARY", a,
          "vertical", {* x : x in residue_vranks *},
          "total", {* x : x in residue_tranks *};
end for;

seed_r := -Q!1/3;
ok_seed, seed_t := TQQ(seed_r);
assert ok_seed;
print "KNOWN_2220_SEED", "r", seed_r, "r_mod_11", Z!(F!seed_r),
      "t", seed_t, "expected_t", -Q!8233/7225,
      "matches", seed_t eq -Q!8233/7225;

print "TOTAL_CONTACT_SOLUTIONS", total_solutions;
print "VERTICAL_RANK_COUNTS", vertical_rank_counts;
print "TOTAL_RANK_COUNTS", total_rank_counts;
if IsDefined(vertical_rank_counts,3) and vertical_rank_counts[3] eq total_solutions then
    print "VERDICT p11_cover_etale_and_locally_surjective_on_all_live_contact_points";
else
    print "VERDICT p11_cover_has_ramified_or_singular_live_contact_points";
end if;
print "CONTACT5_QQ_P11_ETALE_DONE";

quit;
