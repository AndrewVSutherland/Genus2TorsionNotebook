//////////////////////////////////////////////////////////////////////
//  Contact-6 scaffold for M(3,6), then target [6,6].
//
//  Work on an odd genus-2 model y^2=f(x).  Put
//
//      h = 1 + a*x + b*x^2 + x^3,
//      f = h^2 - (x-1)^6.
//
//  Then f has degree <= 5 and f(0)=0.  If the curve is smooth and
//  h(1) != 0, the function h-y has divisor 6P-6*infinity for
//  P=(1,h(1)), so D=P-infinity is a marked order-6 class generically.
//
//  This is not the split M(2,2,2,6) family.  It is the direct contact-6
//  M(6) scaffold.  A [3,6] specialization is detected by requiring an
//  independent rational 3-torsion class; the target [6,6] is detected
//  exactly by TorsionSubgroup, equivalently by halving the independent
//  3-torsion class.
//
//  Typical runs:
//      magma -b mode:=finite prime_bound:=43 code/contact6_m36_search.m
//      magma -b mode:=search height:=20 prime_bound:=43 code/contact6_m36_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then
    mode := "search";
end if;
if not assigned height then
    height := 20;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned prime_bound then
    prime_bound := 43;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if not assigned max_exact then
    max_exact := 200;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned max_hits then
    max_hits := 50;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;
if not assigned progress_interval then
    progress_interval := 100000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            q := Q!num/den;
            key := Sprint(q);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, q);
            end if;
        end for;
    end for;
    return vals;
end function;

function Contact6Polynomial(a, b)
    h := 1 + a*x + b*x^2 + x^3;
    f := h^2 - (x-1)^6;
    return f, h;
end function;

function IntegralModel(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function GoodPolynomial(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function ThreeRankFromInvariants(invs)
    return #[n : n in invs | (Z!n) mod 3 eq 0];
end function;

function TwoRankFromInvariants(invs)
    return #[n : n in invs | (Z!n) mod 2 eq 0];
end function;

function AnyDivisible(invs, d)
    for n in invs do
        if (Z!n) mod d eq 0 then
            return true;
        end if;
    end for;
    return false;
end function;

function Has36(invs)
    return ThreeRankFromInvariants(invs) ge 2 and AnyDivisible(invs, 6);
end function;

function Has66(invs)
    return #[n : n in invs | (Z!n) mod 6 eq 0] ge 2;
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
        try
            fp := ChangeRing(f, GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
                continue;
            end if;
            Lp := LPolynomial(ChangeRing(C, GF(p)));
            fac := Factorization(Lp);
            if #fac eq 1 and fac[1][2] eq 1 and Degree(fac[1][1]) eq 4 then
                return true, p, Lp;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false, 0, P!0;
end function;

function Contact6PolynomialFinite(F, a, b)
    PF<X> := PolynomialRing(F);
    h := 1 + a*X + b*X^2 + X^3;
    f := h^2 - (X-1)^6;
    return f, h;
end function;

function FiniteData(F, a, b)
    f, h := Contact6PolynomialFinite(F, a, b);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false, [], 0, 0;
    end if;
    yP := Evaluate(h, F!1);
    if yP eq 0 or Evaluate(f, F!1) ne yP^2 then
        return false, [], 0, 0;
    end if;
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    A, phi := AbelianGroup(J);
    invs := Invariants(A);
    return true, invs, ThreeRankFromInvariants(invs), TwoRankFromInvariants(invs);
end function;

function ResidueKey(a, b, p)
    F := GF(p);
    try
        aa := F!a;
        bb := F!b;
    catch e
        return false, 0;
    end try;
    ok, invs, trank, tworank := FiniteData(F, aa, bb);
    if not ok then
        return false, 0;
    end if;
    return true, Z!aa + p*(Z!bb);
end function;

function AllowedResidues(p)
    F := GF(p);
    allowed36 := { Z | };
    allowed66 := { Z | };
    good := 0;
    rank_counts := AssociativeArray();
    for a in F do
        for b in F do
            ok, invs, trank, tworank := FiniteData(F, a, b);
            if not ok then
                continue;
            end if;
            good +:= 1;
            keyrank := Sprintf("%o,%o", trank, tworank);
            if IsDefined(rank_counts, keyrank) then
                rank_counts[keyrank] +:= 1;
            else
                rank_counts[keyrank] := 1;
            end if;
            key := Z!a + p*(Z!b);
            if Has36(invs) then
                Include(~allowed36, key);
            end if;
            if Has66(invs) then
                Include(~allowed66, key);
            end if;
        end for;
    end for;
    return allowed36, allowed66, good, rank_counts;
end function;

function PassesResidues(a, b, residue_data, target)
    for data in residue_data do
        p := data[1];
        allowed := target eq "66" select data[3] else data[2];
        open, key := ResidueKey(a, b, p);
        if open and key notin allowed then
            return false, p;
        end if;
    end for;
    return true, 0;
end function;

function ExactData(a, b)
    f, h := Contact6Polynomial(a, b);
    if not GoodPolynomial(f) then
        return false, [], 0, P!0, P!0, 0, 0;
    end if;
    yP := Evaluate(h, Q!1);
    if yP eq 0 or Evaluate(f, Q!1) ne yP^2 then
        return false, [], 0, f, h, yP, 0;
    end if;
    fI, L := IntegralModel(f);
    if not GoodPolynomial(fI) then
        return false, [], 0, fI, h, yP, L;
    end if;
    C := HyperellipticCurve(fI);
    J := Jacobian(C);
    D := J![x-1, L*yP];
    ordD := Order(D);
    G, phi := TorsionSubgroup(J);
    return true, Invariants(G), ordD, fI, h, yP, L;
end function;

procedure FiniteMode()
    primes := [p : p in PrimesUpTo(prime_bound) | p notin {2,3}];
    print "Contact-6 M(6) finite residue data for [3,6] and [6,6]";
    print "prime_bound", prime_bound, "primes", primes;
    for p in primes do
        allowed36, allowed66, good, rank_counts := AllowedResidues(p);
        print "p", p, "good", good,
              "allowed36", #allowed36, "allowed66", #allowed66;
        print " rank_counts", rank_counts;
        if #allowed36 le 20 then
            print " allowed36_residues", Sort(Setseq(allowed36));
        end if;
        if #allowed66 le 20 then
            print " allowed66_residues", Sort(Setseq(allowed66));
        end if;
    end for;
end procedure;

procedure SearchMode()
    primes := [p : p in PrimesUpTo(prime_bound) | p notin {2,3}];
    residue_data := [];
    print "Contact-6 M(6) search for [3,6] and [6,6]";
    print "height", height, "prime_bound", prime_bound,
          "max_exact", max_exact, "max_hits", max_hits;
    print "precomputing residue filters";
    for p in primes do
        allowed36, allowed66, good, rank_counts := AllowedResidues(p);
        Append(~residue_data, <p, allowed36, allowed66>);
        print " p", p, "good", good,
              "allowed36", #allowed36, "allowed66", #allowed66;
    end for;

    params := RationalParametersOfHeight(height);
    checked := 0;
    smooth := 0;
    survivors36 := 0;
    survivors66 := 0;
    exact := 0;
    hits36 := [];
    hits66 := [];
    kill36 := AssociativeArray();
    kill66 := AssociativeArray();

    for a in params do
        for b in params do
            checked +:= 1;
            if progress_interval gt 0 and checked mod progress_interval eq 0 then
                print "progress", checked, "smooth", smooth,
                      "survivors36", survivors36, "survivors66", survivors66,
                      "exact", exact, "hits36", #hits36, "hits66", #hits66;
            end if;

            f, h := Contact6Polynomial(a, b);
            if not GoodPolynomial(f) or Evaluate(h, Q!1) eq 0 then
                continue;
            end if;
            smooth +:= 1;

            pass36, pbad36 := PassesResidues(a, b, residue_data, "36");
            if not pass36 then
                if IsDefined(kill36, pbad36) then kill36[pbad36] +:= 1; else kill36[pbad36] := 1; end if;
                continue;
            end if;
            survivors36 +:= 1;

            pass66, pbad66 := PassesResidues(a, b, residue_data, "66");
            if pass66 then
                survivors66 +:= 1;
            else
                if IsDefined(kill66, pbad66) then kill66[pbad66] +:= 1; else kill66[pbad66] := 1; end if;
            end if;

            if exact ge max_exact then
                continue;
            end if;
            ok, invs, ordD, fI, h0, yP, L := ExactData(a, b);
            if not ok then
                continue;
            end if;
            exact +:= 1;
            if ordD ne 6 then
                continue;
            end if;

            if Has36(invs) then
                simple, pcert, Lp := SimpleCertificate(fI);
                Append(~hits36, <a,b,invs,simple,pcert,fI>);
                print "HIT36", "a", a, "b", b, "invs", invs,
                      "simple", simple, "pcert", pcert;
                print "  f", fI;
            end if;
            if Has66(invs) then
                simple, pcert, Lp := SimpleCertificate(fI);
                Append(~hits66, <a,b,invs,simple,pcert,fI>);
                print "HIT66", "a", a, "b", b, "invs", invs,
                      "simple", simple, "pcert", pcert;
                print "  f", fI;
                if #hits66 ge max_hits then
                    break a;
                end if;
            end if;
        end for;
    end for;

    print "Done";
    print "checked", checked;
    print "smooth", smooth;
    print "survivors36", survivors36;
    print "survivors66", survivors66;
    print "exact", exact;
    print "hits36", #hits36;
    print "hits66", #hits66;
    print "kill36", kill36;
    print "kill66", kill66;
    for H in hits36 do
        print "H36", H;
    end for;
    for H in hits66 do
        print "H66", H;
    end for;
end procedure;

if mode eq "finite" then
    FiniteMode();
elif mode eq "search" then
    SearchMode();
else
    error "unknown mode";
end if;

quit;
