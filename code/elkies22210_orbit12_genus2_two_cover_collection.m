/*
The fake 2-Selmer covering collection for the genus-2 discriminant
quotient at s=59/49.

For speed this script asks Magma to prove the required quintic-field
class-group bounds under GRH.  The five returned classes are nevertheless
checked by the complete local-solubility tests in TwoCoverDescent (no
PrimeBound or PrimeCutoff is supplied).

For a class delta in A=Q[theta]/(F(theta)/14), its genus-5 fake-cover
quotient is the complete intersection of three quadrics obtained from

  delta*(u0+u1 theta+...+u4 theta^4)^2 = c0+c1 theta.

The full two-cover of the genus-2 curve has genus 17; this genus-5 curve
is its quotient by the lifted hyperelliptic involution.

Usage:

  magma -b code/elkies22210_orbit12_genus2_two_cover_collection.m
*/

SetColumns(0);
SetSeed(1);
SetClassGroupBounds("GRH");

Q := Rationals();
PX<X> := PolynomialRing(Q);
F := 14*X^5-775*X^4-4060*X^3+338036*X^2
     +656768*X-44043264;
C := HyperellipticCurve(F);

time fake_set, descent_map := TwoCoverDescent(C);
assert #fake_set eq 5;
A := Domain(descent_map);
theta := A.1;

// The descent-map convention on this odd model sends infinity to 14
// and an affine x-coordinate a to a-theta.
known_elements := [A!14,A!(16-theta),A!(118/7-theta)];
known_classes := [descent_map(a) : a in known_elements];
assert #Setseq(Seqset(known_classes)) eq 3;
assert &and[c in fake_set : c in known_classes];

P4<u0,u1,u2,u3,u4> := ProjectiveSpace(Q,4);
R := CoordinateRing(P4);
us := [R.i : i in [1..5]];

function PrimitiveQuadric(g)
    cs := Coefficients(g);
    den := LCM([Denominator(c) : c in cs]);
    gg := den*g;
    nums := [Abs(Integers()!c) : c in Coefficients(gg) | c ne 0];
    cont := #nums eq 0 select 1 else GCD(nums);
    return (1/cont)*gg;
end function;

function FakeCoverQuotient(delta,A,theta,P4,R,us)
    coeff_forms := [R!0 : k in [0..4]];
    for i,j in [0..4] do
        aij := delta*theta^(i+j);
        cc := Eltseq(aij);
        while #cc lt 5 do Append(~cc,0); end while;
        for k in [0..4] do
            coeff_forms[k+1] +:= (Q!cc[k+1])*us[i+1]*us[j+1];
        end for;
    end for;
    coeff_forms := [PrimitiveQuadric(g) : g in coeff_forms];
    cover := Curve(P4,[coeff_forms[3],coeff_forms[4],coeff_forms[5]]);
    return cover,coeff_forms;
end function;

printf "ELKIES22210_ORBIT12_GENUS2_TWO_COVER_COLLECTION\n";
printf "WARNING_CLASS_GROUP_BOUNDS GRH\n";
printf "fake_two_selmer_set_size %o\n",#fake_set;
printf "known_fake_classes %o\n",known_classes;

unknown_count := 0;
cover_number := 0;
for cl in fake_set do
    cover_number +:= 1;
    known := cl in known_classes;
    delta := known select known_elements[Index(known_classes,cl)]
                   else cl @@ descent_map;
    if not known then unknown_count +:= 1; end if;
    nd := Q!Norm(delta);
    okrho,rho := IsSquare(14*nd);
    assert okrho;

    cover,forms := FakeCoverQuotient(delta,A,theta,P4,R,us);
    assert IsNonsingular(cover);
    assert Genus(cover) eq 5;

    printf "cover %o known %o\n",cover_number,known;
    printf "  selmer_class %o\n",cl;
    printf "  delta %o\n",delta;
    printf "  sqrt_14_norm_delta %o\n",rho;
    printf "  c2 %o\n  c3 %o\n  c4 %o\n",
           forms[3],forms[4],forms[5];
    printf "  genus %o degree %o\n",Genus(cover),Degree(cover);

    // A small search is diagnostic only.  The two unknown covers are
    // everywhere locally soluble by construction, so failure to find a
    // small point is not a proof that they are pointless.
    search_bound := known select 10 else 100;
    time small_points := PointSearch(cover,search_bound : Dimension:=1);
    printf "  point_search_bound %o count %o points %o\n",
           search_bound,#small_points,small_points;
end for;
assert unknown_count eq 2;
printf "known_cover_classes 3\n";
printf "unexplained_everywhere_locally_soluble_cover_classes %o\n",unknown_count;
printf "DONE\n";

quit;
