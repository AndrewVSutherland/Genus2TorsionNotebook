// [2,22] hunt, route 1: scan Flynn's order-11 family
//   F_t(x) = x^6+2x^5+(2t+3)x^4+2x^3+(t^2+1)x^2+2t(1-t)x+t^2
// (D_inf = inf+ - inf- has order 11 for generic t; see
// notes/order22_from_order11.md) for rational t whose sextic factors with
// 2-rank >= 2 WITHOUT a rational branch point -- the [2,2,2] quadratic
// splitting route that notes/order222_from_order11.md never examined
// (its closures assumed a rational branch point).  Since
// [2,22] = (Z/2)^2 x Z/11, a member with 2-rank 2 and the 11-class intact
// IS a [2,22] curve -- no halving needed.
//
// Funnel: factor type prefilter (multiset of factor degrees with >= 2
// even-subset classes) -> exact TorsionSubgroup -> print.  Known-positive
// validation first: t = t_-(2) = -36 must give a linear factor and exact
// torsion [22].
//
// Run: magma -b H:=500 code/flynn11_tworank_scan.m

SetColumns(0);
SetSeed(1);
if not assigned H then H := 500; elif Type(H) eq MonStgElt then H := StringToInteger(H); end if;
if not assigned MemGB then MemGB := 3; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
P<x> := PolynomialRing(Q);

function Flynn(t)
    return x^6 + 2*x^5 + (2*t+3)*x^4 + 2*x^3 + (t^2+1)*x^2 + 2*t*(1-t)*x + t^2;
end function;

// candidate factor-degree multisets with 2-rank >= 2 possible
good_types := { {* 2,2,2 *}, {* 1,1,2,2 *}, {* 1,1,1,1,2 *}, {* 1,1,1,1,1,1 *} };

function ExactTorsionInvs(f)
    C := HyperellipticCurve(f);
    try
        Cm := ReducedMinimalWeierstrassModel(C);
        C := SimplifiedModel(Cm);
    catch e ; end try;
    return Invariants(TorsionSubgroup(Jacobian(C)));
end function;

// VALIDATION: t=-36 (= t_-(2)) has branch point x=4, exact torsion [22]
fval := Flynn(-36);
assert Discriminant(fval) ne 0;
assert Evaluate(fval, 4) eq 0;
invval := ExactTorsionInvs(fval);
printf "VALIDATION t=-36 factor_degrees=%o torsion=%o\n",
       {* Degree(fp[1])^^fp[2] : fp in Factorization(fval) *}, invval;
assert invval eq [22];
printf "VALIDATION_PASS\n";

tested := 0; cands := 0;
for b in [1..H] do
    for a in [-H..H] do
        if GCD(a,b) ne 1 then continue; end if;
        t := Q!a/b;
        if t eq 0 then continue; end if;
        f := Flynn(t);
        if Degree(f) ne 6 or Discriminant(f) eq 0 then continue; end if;
        tested +:= 1;
        if tested mod 100000 eq 0 then
            printf "PROGRESS tested=%o cands=%o\n", tested, cands;
        end if;
        dtype := {* Degree(fp[1])^^fp[2] : fp in Factorization(f) *};
        if dtype notin good_types then continue; end if;
        cands +:= 1;
        invs := ExactTorsionInvs(f);
        printf "CAND t=%o type=%o torsion=%o\n", t, dtype, invs;
        if invs eq [2,22] then
            printf "HIT t=%o TORSION [2,22]  f=%o\n", t, f;
        elif 11 in {i mod 11 eq 0 select 11 else 0 : i in invs} then
            printf "NEAR t=%o torsion=%o\n", t, invs;
        end if;
    end for;
end for;
printf "SEARCH_DONE H=%o tested=%o cands=%o\n", H, tested, cands;
quit;
