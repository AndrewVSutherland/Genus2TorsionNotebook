//////////////////////////////////////////////////////////////////////
//  Contact-5/contact-6 order-30 family: halve its rational 2-class.
//
//  The source family is the one in
//      code/contact5_contact6_order30_family.m.
//  For a smooth specialization its odd quintic factors as q2*q3, with
//  q2 = h6-(x-1)^3 quadratic.  The class T2=[q2,0] is the visible
//  rational 2-torsion class.  If T2 is divisible by 2 over Q, its half
//  has order 4 and combines with the visible order-15 class to give a
//  rational point of order 60.
//
//  Modes:
//
//    Finite-field diagnostic only:
//      magma -b mode:="finite" primes:="7,11,13,17,19,23" \
//          code/contact30_to_60_halving_search.m
//
//    Boundary-aware rational search:
//      magma -b mode:="search" height:=200 prime_bound:=97 \
//          max_exact:=500 MemGB:=8 \
//          code/contact30_to_60_halving_search.m
//
//    Exact one-parameter diagnostic (useful for smoke tests):
//      magma -b mode:="search" parameter_num:=5 parameter_den:=2 \
//          primes:="7" max_exact:=2 crosscheck_cover:=true \
//          verify_torsion_on_nonhits:=true \
//          code/contact30_to_60_halving_search.m
//
//  The finite sieve never kills bad/boundary reductions.  At good odd
//  primes away from 3 and 5 it keeps precisely the residues for which
//  T2 is in 2*J(F_p).  The rational stage first uses a compact exact
//  elimination for the halving cover, then verifies every putative lift
//  with IsDivisibleBy and TorsionSubgroup.  Any order-60 hit receives the
//  same D4/root-power geometric-simplicity test as the repository's full
//  verifier, not merely an irreducible-Frobenius test.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then
    mode := "search";
end if;
if mode notin {"finite", "search"} then
    error "mode must be \"finite\" or \"search\"";
end if;

if not assigned height then
    height := 200;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned prime_bound then
    prime_bound := 97;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if not assigned max_exact then
    max_exact := 500;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;
if not assigned max_print then
    max_print := 50;
elif Type(max_print) eq MonStgElt then
    max_print := StringToInteger(max_print);
end if;
if not assigned progress_interval then
    progress_interval := 10000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;
if not assigned MemGB then
    MemGB := 8;
elif Type(MemGB) eq MonStgElt then
    MemGB := StringToInteger(MemGB);
end if;
SetMemoryLimit(MemGB*10^9);

if not assigned crosscheck_cover then
    crosscheck_cover := false;
elif Type(crosscheck_cover) eq MonStgElt then
    crosscheck_cover := crosscheck_cover in {"true", "True", "1", "yes"};
end if;
if not assigned verify_torsion_on_nonhits then
    verify_torsion_on_nonhits := false;
elif Type(verify_torsion_on_nonhits) eq MonStgElt then
    verify_torsion_on_nonhits :=
        verify_torsion_on_nonhits in {"true", "True", "1", "yes"};
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);
PT<T> := PolynomialRing(Q);

if assigned primes then
    if Type(primes) eq MonStgElt then
        prime_list := [StringToInteger(s) : s in Split(primes, ",") | #s gt 0];
    else
        prime_list := primes;
    end if;
else
    prime_list := [p : p in PrimesUpTo(prime_bound) | p notin {2,3,5}];
end if;
prime_list := [Z!p : p in prime_list | p notin {2,3,5}];

simplicity_primes :=
    [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97];

function RationalParametersOfHeight(B)
    vals := [];
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) eq 1 then
                Append(~vals, Q!num/Q!den);
            end if;
        end for;
    end for;
    return vals;
end function;

function IntegralModel(f)
    L := Z!1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function GoodQuintic(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function MakeMonic(u)
    return u/LeadingCoefficient(u);
end function;

function TorsionOrder(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return &*invs;
end function;

function TorsionExponent(invs)
    e := 1;
    for n in invs do
        e := LCM(e, n);
    end for;
    return e;
end function;

function PointOfExactOrder(G, phi, target)
    for g in G do
        og := Order(g);
        if og mod target eq 0 then
            h := (og div target)*g;
            if Order(h) eq target then
                return true, phi(h), h;
            end if;
        end if;
    end for;
    return false, Codomain(phi)!0, G!0;
end function;

//////////////////////////////////////////////////////////////////////
//  The rational order-30 family and its distinguished quadratic.
//////////////////////////////////////////////////////////////////////

function FamilyQ(rpar, branch)
    denR := rpar^2 - 5;
    if denR eq 0 then
        return false, Q!0, Q!0, Q!0, Q!0, Q!0,
               P!0, P!0, P!0, P!0, P!0;
    end if;

    t := (5*rpar^2 - 20*rpar + 19)/denR;
    Y := -2*(5*rpar^2 - 22*rpar + 25)/denR;
    upar := t^3;
    if upar eq 0 then
        return false, Q!0, Q!0, Q!0, Q!0, Q!0,
               P!0, P!0, P!0, P!0, P!0;
    end if;
    spar := t^5 + t^4 + (Q!5/2)*t^3 + (Q!1/2)*t
          + branch*t*(t-Q!1/2)*(t+1)*Y;

    Cc := (upar^2+1)/(2*upar);
    c := (upar^2-1)/(2*upar);
    if c eq 0 then
        return false, Q!0, Q!0, Q!0, Q!0, Q!0,
               P!0, P!0, P!0, P!0, P!0;
    end if;
    denq := upar^6 + 6*upar^4*spar - 2*upar^4 + 15*upar^3*spar
          - upar*spar^3 + upar^2;
    if denq eq 0 then
        return false, Q!0, Q!0, Q!0, Q!0, Q!0,
               P!0, P!0, P!0, P!0, P!0;
    end if;
    numq := 15*upar^5 + 90*upar^4 + 20*upar^3*spar
          - 6*upar^2*spar^2 + 231*upar^3 + 2*upar^2*spar
          - 15*upar*spar^2 + 90*upar^2 - 20*upar*spar
          + 15*upar - 2*spar;
    qpar := numq/denq;
    A := (spar+qpar)/2;
    e := (spar-qpar)/2;
    B := (15-spar*qpar)/2;
    d := (B*Cc+3)/c;
    h6 := x^3 + A*x^2 + B*x + Cc;
    h5 := e*x^2 + d*x + c;
    f := h6^2 - (x-1)^6;
    q2 := h6 - (x-1)^3;
    K := P!0;
    if Degree(f) eq 5 and (h5^2-f) mod x^5 eq P!0 then
        K := ExactQuotient(h5^2-f, x^5);
    end if;
    return true, t, Y, upar, spar, qpar, h6, h5, f, K, q2;
end function;

function FamilyFp(a, branch, F)
    PF<z> := PolynomialRing(F);
    denR := a^2 - F!5;
    if denR eq 0 then
        return false, PF!0, PF!0;
    end if;
    t := (F!5*a^2 - F!20*a + F!19)/denR;
    Y := -F!2*(F!5*a^2 - F!22*a + F!25)/denR;
    upar := t^3;
    if upar eq 0 then
        return false, PF!0, PF!0;
    end if;
    spar := t^5 + t^4 + (F!5/F!2)*t^3 + (F!1/F!2)*t
          + F!branch*t*(t-F!1/F!2)*(t+1)*Y;
    Cc := (upar^2+1)/(F!2*upar);
    c := (upar^2-1)/(F!2*upar);
    if c eq 0 then
        return false, PF!0, PF!0;
    end if;
    denq := upar^6 + F!6*upar^4*spar - F!2*upar^4
          + F!15*upar^3*spar - upar*spar^3 + upar^2;
    if denq eq 0 then
        return false, PF!0, PF!0;
    end if;
    numq := F!15*upar^5 + F!90*upar^4 + F!20*upar^3*spar
          - F!6*upar^2*spar^2 + F!231*upar^3
          + F!2*upar^2*spar - F!15*upar*spar^2
          + F!90*upar^2 - F!20*upar*spar + F!15*upar - F!2*spar;
    qpar := numq/denq;
    A := (spar+qpar)/F!2;
    B := (F!15-spar*qpar)/F!2;
    h6 := z^3 + A*z^2 + B*z + Cc;
    f := h6^2 - (z-1)^6;
    q2 := h6 - (z-1)^3;
    return true, f, q2;
end function;

// The projective parameter value R=infinity has t=5 and Y=-10.
// It must be included in every local mask: a rational parameter whose
// denominator is divisible by p reduces here, rather than becoming an
// automatically bad/unknown residue.
function FamilyFpInfinity(branch, F)
    PF<z> := PolynomialRing(F);
    t := F!5;
    Y := -F!10;
    upar := t^3;
    if upar eq 0 then
        return false, PF!0, PF!0;
    end if;
    spar := t^5 + t^4 + (F!5/F!2)*t^3 + (F!1/F!2)*t
          + F!branch*t*(t-F!1/F!2)*(t+1)*Y;
    Cc := (upar^2+1)/(F!2*upar);
    c := (upar^2-1)/(F!2*upar);
    if c eq 0 then
        return false, PF!0, PF!0;
    end if;
    denq := upar^6 + F!6*upar^4*spar - F!2*upar^4
          + F!15*upar^3*spar - upar*spar^3 + upar^2;
    if denq eq 0 then
        return false, PF!0, PF!0;
    end if;
    numq := F!15*upar^5 + F!90*upar^4 + F!20*upar^3*spar
          - F!6*upar^2*spar^2 + F!231*upar^3
          + F!2*upar^2*spar - F!15*upar*spar^2
          + F!90*upar^2 - F!20*upar*spar + F!15*upar - F!2*spar;
    qpar := numq/denq;
    A := (spar+qpar)/F!2;
    B := (F!15-spar*qpar)/F!2;
    h6 := z^3 + A*z^2 + B*z + Cc;
    f := h6^2 - (z-1)^6;
    q2 := h6 - (z-1)^3;
    return true, f, q2;
end function;

//////////////////////////////////////////////////////////////////////
//  Compact exact halving cover.
//
//  Write u=x^2+u1*x+u0 and g=f/u.  After putting
//      k=(1/M)^2, z=N/M,
//  coefficient comparison for
//      u*(x+z)^2-k*g = (x^2+r*x+s)^2
//  leaves two quadratic equations F1=F0=0 in z.  Their resultant in z
//  is a univariate polynomial in k.  Only nonzero rational-square k
//  roots are retained, and every reconstructed tuple is checked by the
//  original square-quartic identity before reaching the Jacobian.
//////////////////////////////////////////////////////////////////////

function IsSquareQ(a)
    a := Q!a;
    if a eq 0 then
        return true, Q!0;
    end if;
    if a lt 0 then
        return false, Q!0;
    end if;
    ok, root := IsSquare(a);
    if ok then
        return true, Q!root;
    end if;
    return false, Q!0;
end function;

function SpecializeK(Fz, kval, QZ)
    zz := QZ.1;
    ans := QZ!0;
    for i in [0..Degree(Fz)] do
        ans +:= Q!Evaluate(Coefficient(Fz, i), kval)*zz^i;
    end for;
    return ans;
end function;

function SquareQuarticData(S)
    if Degree(S) ne 4 then
        return false, Q!0, P!0;
    end if;
    s4 := Coefficient(S, 4);
    if s4 eq 0 then
        return false, Q!0, P!0;
    end if;
    s3 := Coefficient(S, 3);
    s2 := Coefficient(S, 2);
    G := x^2 + (s3/(2*s4))*x + (4*s4*s2-s3^2)/(8*s4^2);
    return S eq s4*G^2, Q!s4, G;
end function;

function CompactHalvingCandidates(u, f)
    if Degree(u) ne 2 or LeadingCoefficient(u) eq 0 then
        return false, [* *], "bad_u";
    end if;
    u := MakeMonic(P!u);
    if f mod u ne 0 then
        return false, [* *], "u_not_factor";
    end if;
    g := ExactQuotient(f, u);
    if Degree(g) ne 3 then
        return false, [* *], "bad_g";
    end if;

    u1 := Q!Coefficient(u, 1);
    u0 := Q!Coefficient(u, 0);
    g3 := Q!Coefficient(g, 3);
    g2 := Q!Coefficient(g, 2);
    g1 := Q!Coefficient(g, 1);
    g0 := Q!Coefficient(g, 0);

    QK<k> := PolynomialRing(Q);
    KZ<z> := PolynomialRing(QK);
    r := z + (u1-g3*k)/2;
    s := ((u1+g3*k)*z + u0-g2*k-(u1-g3*k)^2/4)/2;
    F1 := 2*r*s - (u1*z^2+2*u0*z-g1*k);
    F0 := s^2 - (u0*z^2-g0*k);

    res := Resultant(F1, F0);
    if res eq 0 then
        return false, [* *], "zero_resultant";
    end if;

    QZ<zz> := PolynomialRing(Q);
    out := [* *];
    for kr in Roots(res) do
        kval := Q!kr[1];
        if kval eq 0 then
            continue;
        end if;
        square, sqrtk := IsSquareQ(kval);
        if not square or sqrtk eq 0 then
            continue;
        end if;
        f1q := SpecializeK(F1, kval, QZ);
        f0q := SpecializeK(F0, kval, QZ);
        common := GCD(f1q, f0q);
        if Degree(common) lt 1 then
            continue;
        end if;
        for zr in Roots(common) do
            zval := Q!zr[1];
            for sign in [-1,1] do
                Mval := Q!sign/sqrtk;
                Nval := zval*Mval;
                ell := u*(Mval*x+Nval);
                S := ExactQuotient(ell^2-f, u);
                okSquare, lambda, G := SquareQuarticData(S);
                if okSquare then
                    Append(~out, <kval,zval,Mval,Nval,G,ell,lambda>);
                end if;
            end for;
        end for;
    end for;
    return true, out, "ok";
end function;

//////////////////////////////////////////////////////////////////////
//  Finite-field gate for divisibility of the distinguished T2.
//////////////////////////////////////////////////////////////////////

function IsDoubleInFiniteJacobian(J, D)
    G, phi := AbelianGroup(J);
    d := D @@ phi;
    H := sub<G | [2*G.i : i in [1..Ngens(G)]]>;
    return d in H, Invariants(G);
end function;

function FiniteResidueStatus(p, branch)
    F := GF(p);
    allowed := {};
    bad := {};
    good := 0;
    halvable := 0;
    exponent60 := 0;
    errors := 0;

    // Encode infinity by the integer p, so the mask covers P^1(F_p).
    for aa in [0..p] do
        if aa eq p then
            ok, f, q2 := FamilyFpInfinity(branch, F);
        else
            ok, f, q2 := FamilyFp(F!aa, branch, F);
        end if;
        if not ok or not GoodQuintic(f) or Degree(q2) ne 2 or f mod q2 ne 0 then
            Include(~bad, aa);
            continue;
        end if;
        good +:= 1;
        try
            C := HyperellipticCurve(f);
            J := Jacobian(C);
            T2 := J![q2/LeadingCoefficient(q2), F!0];
            isdouble, invs := IsDoubleInFiniteJacobian(J, T2);
            if isdouble then
                Include(~allowed, aa);
                halvable +:= 1;
            end if;
            if TorsionExponent(invs) mod 60 eq 0 then
                exponent60 +:= 1;
            end if;
        catch e
            // An implementation failure is unknown, never a local kill.
            Include(~bad, aa);
            errors +:= 1;
        end try;
    end for;
    return allowed, bad, good, halvable, exponent60, errors;
end function;

function ResidueKey(p, branch)
    return IntegerToString(p) cat ":" cat IntegerToString(branch);
end function;

function ResidueOfRational(rpar, p)
    num := Numerator(rpar);
    den := Denominator(rpar);
    if den mod p eq 0 then
        return true, p;
    end if;
    F := GF(p);
    return true, Z!(F!num/F!den);
end function;

function PassesFiniteHalving(rpar, branch, allowed, badres, plist)
    bad_primes := [];
    good_checked := 0;
    for p in plist do
        ok, residue := ResidueOfRational(rpar, p);
        if not ok then
            Append(~bad_primes, p);
            continue;
        end if;
        key := ResidueKey(p, branch);
        if residue in badres[key] then
            Append(~bad_primes, p);
            continue;
        end if;
        good_checked +:= 1;
        if residue notin allowed[key] then
            return false, p, good_checked, bad_primes;
        end if;
    end for;
    return true, 0, good_checked, bad_primes;
end function;

//////////////////////////////////////////////////////////////////////
//  Full geometric-simplicity certificate.
//////////////////////////////////////////////////////////////////////

function FrobeniusPolynomial(C, p)
    ef := EulerFactor(C, p);
    d := Degree(ef);
    return &+[Q!Coefficient(ef, i)*T^(d-i) : i in [0..d]];
end function;

function FullSimplicityCertificate(C, plist)
    // First try the Leprévost D4 criterion.
    for p in plist do
        try
            Phi := FrobeniusPolynomial(C, p);
            fac := Factorization(Phi);
            if #fac ne 1 or Degree(fac[1][1]) ne 4 or fac[1][2] ne 1 then
                continue;
            end if;
            Gal := GaloisGroup(Phi);
            try
                desc := TransitiveGroupDescription(Gal);
            catch e2
                desc := "not transitive";
            end try;
            if Order(Gal) eq 8 and desc eq "D(4)" then
                return true, "D4", p, Phi;
            end if;
        catch e
            continue;
        end try;
    end for;

    // Fallback: an absolutely simple reduction, checked through root
    // powers 2,...,12 as in verify_simple_torsion_candidate.m.
    for p in plist do
        try
            Phi := FrobeniusPolynomial(C, p);
            fac := Factorization(Phi);
            if Degree(Phi) ne 4 or #fac ne 1 or
               Degree(fac[1][1]) ne 4 or fac[1][2] ne 1 then
                continue;
            end if;
            K<pi> := NumberField(Phi);
            power_ok := true;
            for n in [2..12] do
                if Degree(MinimalPolynomial(pi^n)) lt 4 then
                    power_ok := false;
                    break;
                end if;
            end for;
            if power_ok then
                return true, "root_power", p, Phi;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false, "none", 0, PT!0;
end function;

//////////////////////////////////////////////////////////////////////
//  Precompute and report the finite masks.
//////////////////////////////////////////////////////////////////////

allowed := AssociativeArray();
badres := AssociativeArray();

print "contact30 -> 60 distinguished-2 halving search";
print "mode", mode, "height", height, "primes", prime_list,
      "max_exact", max_exact, "MemGB", MemGB,
      "crosscheck_cover", crosscheck_cover;
print "FINITE_DIAGNOSTIC_BEGIN";
for p in prime_list do
    for branch in [-1,1] do
        good_allowed, bad, good, halvable, exp60, errors :=
            FiniteResidueStatus(p, branch);
        key := ResidueKey(p, branch);
        allowed[key] := good_allowed;
        badres[key] := bad;
        print "FINITE", "p", p, "branch", branch,
              "total", p+1, "good", good, "halvable", halvable,
              "allowed", #good_allowed, "bad", #bad,
              "exponent60", exp60, "errors", errors;
    end for;
end for;
print "FINITE_DIAGNOSTIC_END";

if mode eq "finite" then
    quit;
end if;

//////////////////////////////////////////////////////////////////////
//  Boundary-aware rational search and exact verification.
//////////////////////////////////////////////////////////////////////

if assigned parameter_num or assigned parameter_den then
    if not (assigned parameter_num and assigned parameter_den) then
        error "parameter_num and parameter_den must be supplied together";
    end if;
    if Type(parameter_num) eq MonStgElt then
        parameter_num := StringToInteger(parameter_num);
    end if;
    if Type(parameter_den) eq MonStgElt then
        parameter_den := StringToInteger(parameter_den);
    end if;
    if parameter_den eq 0 then
        error "parameter_den must be nonzero";
    end if;
    params := [Q!parameter_num/Q!parameter_den];
else
    params := RationalParametersOfHeight(height);
end if;

parameter_count := #params;
branch_checked := 0;
family_ok := 0;
smooth := 0;
finite_survivors := 0;
cover_tests := 0;
cover_positive := 0;
cover_errors := 0;
exact_tests := 0;
contact30_verified := 0;
halvable := 0;
cover_mismatches := 0;
torsion_tests := 0;
order60_hits := 0;
simple_hits := 0;
printed := 0;
hits := [* *];
kill_counts := AssociativeArray();
bad_signature_counts := AssociativeArray();
stop := false;

for rpar in params do
    if stop then
        break;
    end if;
    for branch in [-1,1] do
        branch_checked +:= 1;
        if progress_interval gt 0 and branch_checked mod progress_interval eq 0 then
            print "PROGRESS", "branches", branch_checked,
                  "smooth", smooth, "finite", finite_survivors,
                  "cover_positive", cover_positive,
                  "exact", exact_tests, "halvable", halvable,
                  "order60", order60_hits, "simple", simple_hits;
        end if;

        pass, killp, good_checked, bad_primes :=
            PassesFiniteHalving(rpar, branch, allowed, badres, prime_list);
        if not pass then
            if IsDefined(kill_counts, killp) then
                kill_counts[killp] +:= 1;
            else
                kill_counts[killp] := 1;
            end if;
            continue;
        end if;

        ok, t, Y, upar, spar, qpar, h6, h5, f, Kcontact, q2 :=
            FamilyQ(rpar, branch);
        if not ok then
            continue;
        end if;
        family_ok +:= 1;
        if not GoodQuintic(f) or Degree(q2) ne 2 or f mod q2 ne 0 then
            continue;
        end if;
        smooth +:= 1;
        finite_survivors +:= 1;

        sig := Join([IntegerToString(p) : p in bad_primes], ",");
        if IsDefined(bad_signature_counts, sig) then
            bad_signature_counts[sig] +:= 1;
        else
            bad_signature_counts[sig] := 1;
        end if;

        q2 := MakeMonic(q2);
        cover_tests +:= 1;
        cover_ok := false;
        cover_data := [* *];
        cover_status := "not_run";
        try
            cover_computed, cover_data, cover_status :=
                CompactHalvingCandidates(q2, f);
            cover_ok := cover_computed and #cover_data gt 0;
        catch e
            cover_errors +:= 1;
            cover_status := "exception";
        end try;
        if cover_ok then
            cover_positive +:= 1;
        end if;

        if printed lt max_print then
            print "SURVIVOR", "R", rpar, "branch", branch,
                  "good_checked", good_checked, "bad_primes", bad_primes,
                  "cover", cover_ok, "cover_count", #cover_data,
                  "cover_status", cover_status;
            printed +:= 1;
        end if;

        // In production, exact Jacobian work is needed only for a positive
        // compact cover.  crosscheck_cover=true deliberately also checks
        // cover-negative cases, to validate the elimination on small runs.
        if not cover_ok and not crosscheck_cover then
            continue;
        end if;
        if exact_tests ge max_exact then
            continue;
        end if;

        exact_tests +:= 1;
        try
            fI, L := IntegralModel(f);
            C := HyperellipticCurve(fI);
            J := Jacobian(C);

            D5 := J![x, L*Evaluate(h5, 0)];
            D6 := J![x-1, L*Evaluate(h6, 1)];
            ord5 := Order(D5);
            ord6 := Order(D6);
            ord30 := Order(D5+D6);
            if ord5 eq 5 and ord6 eq 6 and ord30 eq 30 then
                contact30_verified +:= 1;
            else
                print "CONTACT_SANITY_FAIL", "R", rpar, "branch", branch,
                      "orders", [ord5,ord6,ord30];
                continue;
            end if;

            T2 := J![q2, Q!0];
            if Order(T2) ne 2 then
                print "T2_SANITY_FAIL", "R", rpar, "branch", branch,
                      "order", Order(T2);
                continue;
            end if;
            isdiv, half := IsDivisibleBy(T2, 2);
            if isdiv then
                halvable +:= 1;
            end if;
            if isdiv ne cover_ok then
                cover_mismatches +:= 1;
                print "COVER_MISMATCH", "R", rpar, "branch", branch,
                      "cover", cover_ok, "IsDivisibleBy", isdiv;
            end if;

            // Reconstruct the first compact half on the integral model and
            // check it directly.  Under y_integral=L*y, ell is multiplied by L.
            if cover_ok then
                cd := cover_data[1];
                Ghalf := P!cd[5];
                ellI := P!(L*cd[6]);
                Hcover := J![Ghalf, -ellI mod Ghalf];
                if 2*Hcover ne T2 or Order(Hcover) ne 4 then
                    print "RECONSTRUCTION_FAIL", "R", rpar, "branch", branch;
                    continue;
                end if;
            end if;

            if not isdiv and not verify_torsion_on_nonhits then
                print "EXACT", "R", rpar, "branch", branch,
                      "contact30", true, "T2_divisible", false;
                continue;
            end if;

            Gtors, phi := TorsionSubgroup(J);
            torsion_tests +:= 1;
            invs := Invariants(Gtors);
            ordtors := TorsionOrder(invs);
            exptors := TorsionExponent(invs);
            print "EXACT", "R", rpar, "branch", branch,
                  "contact30", true, "T2_divisible", isdiv,
                  "half_order", isdiv select Order(half) else 0,
                  "torsion", invs, "order", ordtors, "exponent", exptors,
                  "bad_primes", bad_primes;

            if isdiv and exptors mod 60 eq 0 then
                order60_hits +:= 1;
                ok60, P60, g60 := PointOfExactOrder(Gtors, phi, 60);
                simple, cert_kind, pcert, Phi :=
                    FullSimplicityCertificate(C, simplicity_primes);
                if simple then
                    simple_hits +:= 1;
                end if;
                Append(~hits, <rpar,branch,invs,fI,half,ok60,P60,
                                simple,cert_kind,pcert,Phi>);
                print "ORDER60_HIT", "R", rpar, "branch", branch,
                      "torsion", invs, "simple", simple,
                      "certificate", cert_kind, "p", pcert, "Phi", Phi;
                print "  f =", fI;
                print "  half_T2 =", half;
                if ok60 then
                    print "  P60 =", P60;
                    print "  abstract_g60 =", g60;
                end if;
                if #hits ge max_hits then
                    stop := true;
                    break;
                end if;
            end if;
        catch e
            print "EXACT_ERROR", "R", rpar, "branch", branch, e`Object;
        end try;
    end for;
end for;

print "DONE contact30 -> 60";
print "parameters", parameter_count,
      "branch_checked", branch_checked,
      "family_ok", family_ok,
      "smooth", smooth,
      "finite_survivors", finite_survivors,
      "cover_tests", cover_tests,
      "cover_positive", cover_positive,
      "cover_errors", cover_errors,
      "exact_tests", exact_tests,
      "contact30_verified", contact30_verified,
      "halvable", halvable,
      "cover_mismatches", cover_mismatches,
      "torsion_tests", torsion_tests,
      "order60_hits", order60_hits,
      "simple_hits", simple_hits;
print "KILL_COUNTS";
for p in Sort([k : k in Keys(kill_counts)]) do
    print p, kill_counts[p];
end for;
print "BAD_SIGNATURE_COUNTS";
for key in Sort([k : k in Keys(bad_signature_counts)]) do
    print key, bad_signature_counts[key];
end for;
for hit in hits do
    print "HIT_RECORD", hit;
end for;

quit;
