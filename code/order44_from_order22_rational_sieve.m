//////////////////////////////////////////////////////////////////////
// Conservative multi-prime rational-parameter sieve for order 44.
//
// For each family/sign/prime, the allowed P^1(F_p) mask contains:
//   * exact open marked-class halves;
//   * singular or collapsed-order open reductions (kept);
//   * s=0,+/-1,infinity (kept).
// Good open points where the marked order-22 class is not in 2J(F_p)
// are the only rejected residues.  Reduced s=a/b passing every mask
// is checked exactly over Q with IsDivisibleBy(D,2).
//
// Typical bounded run:
//   magma -b height:=300 primes:="3,5,7,13,17,19,23,29,31" \
//       code/order44_from_order22_rational_sieve.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

if not assigned height then height := 300; end if;
if Type(height) eq MonStgElt then height := StringToInteger(height); end if;
if assigned primes then
    if Type(primes) eq MonStgElt then
        prime_list := [StringToInteger(a) : a in Split(primes, ",") | #a gt 0];
    else
        prime_list := primes;
    end if;
else
    prime_list := [3,5,7,13,17,19,23,29,31];
end if;

function IsDivisibleBy2Finite(D,G,phi)
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

function FamilyPolynomial(family,P,par)
    x := P.1;
    if family eq "Flynn" then
        return x^6 + 2*x^5 + (2*par+3)*x^4 + 2*x^3
               + (par^2+1)*x^2 + 2*par*(1-par)*x + par^2;
    end if;
    return x^6 - 4*x^5 + 8*(1+par)*x^4 - (10+32*par)*x^3
           + 8*(1+6*par+2*par^2)*x^2
           - 4*(1+6*par+16*par^2)*x + 64*par^2+1;
end function;

function FamilyData(family,P,s,e)
    // The cancelled formula is valid at s=e and exposes the genuine
    // finite member there.  It retains its sole pole at s=-e.
    if s eq -e then return false,_,_,_; end if;
    z := s*(s^2+e*s+1)/(s+e);
    if family eq "Flynn" then
        par := -z^2;
        r := s^2;
    else
        par := -z^2/4;
        r := 1+s^2;
    end if;
    f := FamilyPolynomial(family,P,par);
    return true,f,r,par;
end function;

function MoveMarkedRootToInfinity(f,r)
    P := Parent(f);
    X := P.1;
    return P!(&+[Coefficient(f,i)*(r*X+1)^i*X^(6-i) : i in [0..6]]);
end function;

function IntegralSquareModel(f)
    den := Z!1;
    for c in Coefficients(f) do
        den := LCM(den,Denominator(c));
    end for;
    return Parent(f)!(den^2*f),den;
end function;

function MaskKey(family,e,j)
    return Sprintf("%o_%o_%o",family,e,j);
end function;

function BuildAllowedMask(family,e,p)
    F := GF(p);
    P<X> := PolynomialRing(F);
    allowed := { Z | p,0,1,p-1 }; // p encodes infinity.
    good_rejected := 0;
    exact_open_halves := 0;
    kept_bad := 0;

    for s in F do
        if s eq 0 or s^2 eq 1 then continue; end if;
        ok,f,r,par := FamilyData(family,P,s,F!e);
        assert ok;
        if Degree(f) ne 6 or Evaluate(f,r) ne 0 or Discriminant(f) eq 0 then
            Include(~allowed,Z!s);
            kept_bad +:= 1;
            continue;
        end if;
        g := MoveMarkedRootToInfinity(f,r);
        if Degree(g) ne 5 or Coefficient(g,0) ne 1 or Discriminant(g) eq 0 then
            Include(~allowed,Z!s);
            kept_bad +:= 1;
            continue;
        end if;
        J := Jacobian(HyperellipticCurve(g));
        D := J![X,F!1];
        if Order(D) ne 22 then
            Include(~allowed,Z!s);
            kept_bad +:= 1;
            continue;
        end if;
        G,phi := AbelianGroup(J);
        if IsDivisibleBy2Finite(D,G,phi) then
            Include(~allowed,Z!s);
            exact_open_halves +:= 1;
        else
            good_rejected +:= 1;
        end if;
    end for;
    return allowed,good_rejected,exact_open_halves,kept_bad;
end function;

if assigned only_flynn then
    families := ["Flynn"];
else
    families := ["Flynn","DaowsudSchmidt"];
end if;
signs := [-1,1];
masks := AssociativeArray();

print "ORDER44_FROM_ORDER22_RATIONAL_SIEVE";
print "height",height,"primes",prime_list;
print "MASK_BUILD";
for family in families do
    for e in signs do
        for j in [1..#prime_list] do
            p := prime_list[j];
            require IsPrime(p) and p ne 2 and p ne 11:
                "primes must be odd and different from 11";
            mask,rejected,halves,bad := BuildAllowedMask(family,e,p);
            masks[MaskKey(family,e,j)] := mask;
            print "mask",family,"sign",e,"p",p,
                  "allowed_P1",#mask,"of",p+1,
                  "good_rejected",rejected,
                  "open_half_residues",halves,"bad_open_kept",bad,
                  "allowed_codes",Sort(Setseq(mask));
        end for;
    end for;
end for;

stage_counts := [Z!0 : j in [1..#prime_list]];
reduced_rationals := 0;
family_sign_trials := 0;
sieve_survivors := 0;
exact_smooth := 0;
exact_boundary_skipped := 0;
exact_singular := 0;
exact_order_collapse := 0;
exact_errors := 0;
exact_halves := 0;

PQ<xq> := PolynomialRing(Q);

for b in [1..height] do
    for a in [-height..height] do
        if GCD(Abs(a),b) ne 1 then continue; end if;
        reduced_rationals +:= 1;
        sQ := Q!a/b;
        for family in families do
            for eZ in signs do
                family_sign_trials +:= 1;
                passed := true;
                for j in [1..#prime_list] do
                    p := prime_list[j];
                    if b mod p eq 0 then
                        code := p;
                    else
                        code := Z!((GF(p)!a)/(GF(p)!b));
                    end if;
                    if not code in masks[MaskKey(family,eZ,j)] then
                        passed := false;
                        break;
                    end if;
                    stage_counts[j] +:= 1;
                end for;
                if not passed then continue; end if;
                sieve_survivors +:= 1;

                // s=0 and s=-e are genuine degenerate boundary values.
                // s=e is retained: cancellation gives a smooth member.
                if sQ eq 0 or sQ eq -eZ then
                    exact_boundary_skipped +:= 1;
                    continue;
                end if;
                ok,f,r,par := FamilyData(family,PQ,sQ,Q!eZ);
                if not ok or Degree(f) ne 6 or Evaluate(f,r) ne 0
                   or Discriminant(f) eq 0 then
                    exact_singular +:= 1;
                    continue;
                end if;
                g := MoveMarkedRootToInfinity(f,r);
                if Degree(g) ne 5 or Coefficient(g,0) ne 1 or Discriminant(g) eq 0 then
                    exact_singular +:= 1;
                    continue;
                end if;
                exact_smooth +:= 1;
                try
                    gI,yscale := IntegralSquareModel(g);
                    J := Jacobian(HyperellipticCurve(gI));
                    D := J![xq,Q!yscale];
                    if Order(D) ne 22 then
                        exact_order_collapse +:= 1;
                        continue;
                    end if;
                    div2,H := IsDivisibleBy(D,2);
                    if div2 then
                        exact_halves +:= 1;
                        print "EXACT_ORDER44_HIT",family,"sign",eZ,
                              "s",sQ,"parameter",par,
                              "D_order",Order(D),"half_order",Order(H),
                              "curve",f;
                    end if;
                catch err
                    exact_errors +:= 1;
                    print "EXACT_ERROR",family,"sign",eZ,"s",sQ,err`Object;
                end try;
            end for;
        end for;
    end for;
end for;

print "SIEVE_DONE";
print "reduced_rationals",reduced_rationals,
      "family_sign_trials",family_sign_trials;
print "stage_counts_after_each_prime",[
    <prime_list[j],stage_counts[j]> : j in [1..#prime_list]];
print "sieve_survivors",sieve_survivors,
      "exact_boundary_skipped",exact_boundary_skipped,
      "exact_smooth",exact_smooth,
      "exact_singular",exact_singular,
      "exact_order_collapse",exact_order_collapse,
      "exact_errors",exact_errors,
      "exact_halves",exact_halves;

quit;
