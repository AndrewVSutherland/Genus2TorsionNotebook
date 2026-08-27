//////////////////////////////////////////////////////////////////////
//  Elkies [2,2,2,10] source-halving cover over finite fields.
//
//  The Clebsch--Klein source is
//
//      sum r_i = sum r_i^3 = 0,
//      C_r : y^2 = x * product_i (x-r_i^2).
//
//  On the open set r_i != 0 with distinct r_i^2, C_r has full
//  rational 2-torsion and the marked 5-torsion.  The 15 nonzero
//  2-classes split into two S_5-orbits:
//
//      orbit 0i : {0,r_i^2},       size 5;
//      orbit ij : {r_i^2,r_j^2},   size 10.
//
//  Send alpha_j to infinity and alpha_i to zero.  Monicizing the
//  resulting odd quintic and applying the (x-T) descent shows that
//  {alpha_i,alpha_j} is divisible by 2 iff the four quantities
//
//   (alpha_k-alpha_i) * product_{l != i,j,k}(alpha_j-alpha_l)
//
//  are squares.  This script counts the two marked covers and also
//  validates the criterion against AbelianGroup(J) on a few fibers.
//
//  Typical run:
//      magma -b validate_limit:=10 \
//          code/elkies22210_source_halving_finite.m
//////////////////////////////////////////////////////////////////////

// For command-line ranges use, for example,
//   magma -b prime_min:=5 prime_max:=101 ...
// This avoids Magma's conversion of a command-line sequence assignment
// such as prime_list:="[11,13]" into a string.
if assigned prime_min or assigned prime_max then
    if not assigned prime_min then prime_min := 5;
    elif Type(prime_min) eq MonStgElt then
        prime_min := StringToInteger(prime_min);
    end if;
    if not assigned prime_max then prime_max := 47;
    elif Type(prime_max) eq MonStgElt then
        prime_max := StringToInteger(prime_max);
    end if;
    prime_list := [p : p in [Maximum(5,prime_min)..prime_max] | IsPrime(p)];
elif not assigned prime_list then
    prime_list := [5,7,11,13,17,19,23,29,31,37,41,43,47];
end if;

if not assigned validate_limit then
    validate_limit := 8;
elif Type(validate_limit) eq MonStgElt then
    validate_limit := StringToInteger(validate_limit);
end if;

function IsOpenTuple(rs)
    if &or [r eq 0 : r in rs] then
        return false;
    end if;
    return #Setseq(Seqset([r^2 : r in rs])) eq 5;
end function;

function PairHalves(branches, i, j)
    // This is literally -D*(c-B)/(c-A), with
    // A=branches[i], B=branches[j], D=Product_c(A-c), after
    // cancelling the factor A-c.  Keeping the common squareclass is
    // essential; cross-ratio quotients alone give false positives.
    rem := [k : k in [1..6] | k ne i and k ne j];
    for k in rem do
        z := branches[k]-branches[j];
        for ell in rem do
            if ell ne k then
                z *:= branches[i]-branches[ell];
            end if;
        end for;
        if z eq 0 or not IsSquare(z) then
            return false;
        end if;
    end for;
    return true;
end function;

function IsDivisibleBy2Finite(J, D)
    G, phi := AbelianGroup(J);
    a := D @@ phi;
    coords := Eltseq(a);
    invs := Invariants(G);
    for i in [1..#coords] do
        if GCD(2,invs[i]) ne 1 and coords[i] mod 2 ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

procedure ValidatePairCriterion(F, rs, pair)
    P<x> := PolynomialRing(F);
    branches := [F!0] cat [r^2 : r in rs];
    f := &*[x-a : a in branches];
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    i := pair[1];
    j := pair[2];
    T := J![(x-branches[i])*(x-branches[j]),F!0];
    criterion := PairHalves(branches,i,j);
    finite_group := IsDivisibleBy2Finite(J,T);
    assert criterion eq finite_group;
end procedure;

print "ELKIES22210_SOURCE_HALVING_FINITE";
print "prime_list", prime_list, "validate_limit", validate_limit;

for p in prime_list do
    F := GF(p);
    total_ck := 0;
    open_ck := 0;
    marked_0i := 0;
    marked_ij := 0;
    any_0i := 0;
    any_ij := 0;
    total_0i_pairs := 0;
    total_ij_pairs := 0;
    validated := 0;
    validated_positive_0i := 0;
    validated_positive_ij := 0;
    sample_0i := [];
    sample_ij := [];

    // On the open set r1 is nonzero, so projective scaling sets r1=1.
    // Brute force is retained in characteristic 3 only; the requested
    // prime list starts at 5, where the quadratic solve is valid.
    r1 := F!1;
    for r2 in F do
        for r3 in F do
            s := r1+r2+r3;
            A := r1^3+r2^3+r3^3;
            if Characteristic(F) eq 3 then
                r4s := [z : z in F];
            elif s eq 0 then
                r4s := [];
            else
                PR<R4> := PolynomialRing(F);
                r4s := [rt[1] : rt in Roots(3*s*R4^2+3*s^2*R4+(s^3-A))];
            end if;

            for r4 in r4s do
                r5 := -(r1+r2+r3+r4);
                rs := [r1,r2,r3,r4,r5];
                if &+rs ne 0 or &+[r^3 : r in rs] ne 0 then
                    continue;
                end if;
                total_ck +:= 1;
                if not IsOpenTuple(rs) then
                    continue;
                end if;
                open_ck +:= 1;
                branches := [F!0] cat [r^2 : r in rs];

                h0 := PairHalves(branches,1,2);
                h1 := PairHalves(branches,2,3);
                if h0 then
                    marked_0i +:= 1;
                    if #sample_0i lt 3 then Append(~sample_0i,rs); end if;
                end if;
                if h1 then
                    marked_ij +:= 1;
                    if #sample_ij lt 3 then Append(~sample_ij,rs); end if;
                end if;

                count0 := 0;
                for i in [2..6] do
                    if PairHalves(branches,1,i) then count0 +:= 1; end if;
                end for;
                count1 := 0;
                for i in [2..5] do
                    for j in [i+1..6] do
                        if PairHalves(branches,i,j) then count1 +:= 1; end if;
                    end for;
                end for;
                total_0i_pairs +:= count0;
                total_ij_pairs +:= count1;
                if count0 gt 0 then any_0i +:= 1; end if;
                if count1 gt 0 then any_ij +:= 1; end if;

                if validated lt validate_limit then
                    ValidatePairCriterion(F,rs,<1,2>);
                    ValidatePairCriterion(F,rs,<2,3>);
                    validated +:= 1;
                end if;
                if h0 and validated_positive_0i eq 0 then
                    ValidatePairCriterion(F,rs,<1,2>);
                    validated_positive_0i := 1;
                end if;
                if h1 and validated_positive_ij eq 0 then
                    ValidatePairCriterion(F,rs,<2,3>);
                    validated_positive_ij := 1;
                end if;
            end for;
        end for;
    end for;

    print "PRIME",p,
          "ck",total_ck,"open",open_ck,
          "marked_0i",marked_0i,"marked_ij",marked_ij,
          "any_0i",any_0i,"any_ij",any_ij,
          "pair_0i",total_0i_pairs,"pair_ij",total_ij_pairs,
          "validated",validated,
          "validated_positive",<validated_positive_0i,validated_positive_ij>;
    print "SAMPLES_0i",sample_0i;
    print "SAMPLES_ij",sample_ij;
    if open_ck gt 0 and marked_0i eq 0 then
        print "NO_GOOD_OPEN_MARKED_0i_MOD",p;
    end if;
    if open_ck gt 0 and marked_ij eq 0 then
        print "NO_GOOD_OPEN_MARKED_ij_MOD",p;
    end if;
end for;

print "DONE";
quit;
