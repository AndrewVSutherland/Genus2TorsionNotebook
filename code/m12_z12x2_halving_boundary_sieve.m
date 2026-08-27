//////////////////////////////////////////////////////////////////////
//  Boundary sieve for halving the independent 2-torsion point in the
//  one-parameter M(12) family a=(1-r)/4.
//
//  The good-reduction sieve modulo 7 leaves only boundary residues.
//  This script records:
//    1. the tangent-equation behavior modulo 7;
//    2. allowed affine residues at several auxiliary primes, where
//       "allowed" means either bad/boundary reduction or good reduction
//       with the independent 2-torsion divisible by 2 over F_p;
//    3. a rational height search using those residue conditions before
//       running exact IsDivisibleBy on the surviving nonsingular curves.
//
//  Typical run from torsion_jac:
//      magma -b height:=300 code/m12_z12x2_halving_boundary_sieve.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 300;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

Q := Rationals();
P<X> := PolynomialRing(Q);

function IsDivisibleBy2Finite(J, D)
    G, phi := AbelianGroup(J);
    a := D @@ phi;
    coords := Eltseq(a);
    invs := Invariants(G);
    for i in [1..#coords] do
        if GCD(2, invs[i]) ne 1 and coords[i] mod 2 ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function OddQuinticForR(r)
    a := (1-r)/4;
    T := a*X^2 - X + r;
    h := (X-r)*(T+1);
    W := h^2 + 4*a*X^2*T*(T+1);

    f5 := P!0;
    for i in [0..Degree(W)] do
        for j in [0..i] do
            f5 +:= Coefficient(W, i)*Binomial(i,j)*2^(i-j)*X^(6-j);
        end for;
    end for;
    return f5;
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function IrreducibleFrobeniusCertificate(f)
    C := HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47] do
        try
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

function IndependentTangentEquationSolutionCountMod7(rv)
    F := GF(7);
    if rv eq F!0 or rv eq F!1 then
        return -1, [];
    end if;

    R<A,B,m,n> := PolynomialRing(F, 4);
    r := F!rv;
    eq0 := (4*r^2-4*r)*B^2 + ((-r/4+1/2)/(r-1))*n^2
           - 1/4*r^3 + 3/4*r^2 - 3/4*r + 1/4;
    eq1 := (8*r^2-8*r)*A*B + ((-1/2*r+1)/(r-1))*m*n
           - n^2 - 5/2*r^3 + 11/2*r^2 - 7/2*r + 1/2;
    eq2 := (4*r^2-4*r)*A^2 + ((-1/4*r+1/2)/(r-1))*m^2
           - 2*m*n + (8*r^2-8*r)*B
           + 1/4*r^4 - 33/4*r^3 + 45/4*r^2 - 11/4*r - 1/2;
    eq3 := -m^2 + (8*r^2-8*r)*A + r^4 - 8*r^3 + 4*r^2 + 4*r - 1;

    sols := [];
    for av in F do
        for bv in F do
            for mv in F do
                for nv in F do
                    if Evaluate(eq0, [av,bv,mv,nv]) eq 0
                       and Evaluate(eq1, [av,bv,mv,nv]) eq 0
                       and Evaluate(eq2, [av,bv,mv,nv]) eq 0
                       and Evaluate(eq3, [av,bv,mv,nv]) eq 0 then
                        Append(~sols, <av,bv,mv,nv>);
                    end if;
                end for;
            end for;
        end for;
    end for;
    return #sols, sols;
end function;

function AllowedResidues(p)
    F := GF(p);
    PF<Y> := PolynomialRing(F);
    allowed := {};

    for rr in F do
        if rr in {F!0, F!1, F!2} then
            Include(~allowed, Integers()!rr);
            continue;
        end if;

        a := (1-rr)/4;
        T := a*Y^2 - Y + rr;
        h := (Y-rr)*(T+1);
        W := h^2 + 4*a*Y^2*T*(T+1);

        f5 := PF!0;
        for i in [0..Degree(W)] do
            for j in [0..i] do
                f5 +:= Coefficient(W, i)*Binomial(i,j)*2^(i-j)*Y^(6-j);
            end for;
        end for;

        if Degree(f5) ne 5 or Discriminant(f5) eq 0 then
            Include(~allowed, Integers()!rr);
            continue;
        end if;

        beta_ind := (2-rr)/(4*(rr-1));
        C := HyperellipticCurve(f5);
        J := Jacobian(C);
        Tind := J![Y-beta_ind, F!0];
        if IsDivisibleBy2Finite(J, Tind) then
            Include(~allowed, Integers()!rr);
        end if;
    end for;
    return allowed;
end function;

print "M(12) Z/12 x Z/4 boundary sieve";
print "height", height;

print "mod 7 tangent-equation counts for independent beta:";
for rv in [0..6] do
    cnt, sols := IndependentTangentEquationSolutionCountMod7(GF(7)!rv);
    if cnt lt 0 then
        print "  r", rv, "coordinate boundary";
    else
        print "  r", rv, "solutions", cnt,
              "sample", sols[1..Minimum(#sols, 4)];
    end if;
end for;

primes := [7,11,19,23,31,43,47,59,67,71];
allowed := AssociativeArray(Integers());
for p in primes do
    allowed[p] := AllowedResidues(p);
    print "p", p, "allowed", #allowed[p], Sort(Setseq(allowed[p]));
end for;

seen := {};
candidates := [];
exact_tests := 0;
hits := [];
total := 0;

for den in [1..height] do
    for num in [-height..height] do
        if GCD(num, den) ne 1 then
            continue;
        end if;
        r := Q!num/den;
        key := Sprint(r);
        if key in seen then
            continue;
        end if;
        Include(~seen, key);
        total +:= 1;

        if r in {Q!0, Q!1, Q!2} then
            continue;
        end if;

        pass := true;
        for p in primes do
            // If p divides the denominator, this is the infinity residue,
            // treated here as boundary/allowed rather than excluded.
            if den mod p eq 0 then
                continue;
            end if;
            res := Integers()!((Integers(p)!num)*(Integers(p)!den)^-1);
            if res notin allowed[p] then
                pass := false;
                break;
            end if;
        end for;
        if not pass then
            continue;
        end if;
        Append(~candidates, r);

        f5 := OddQuinticForR(r);
        if Degree(f5) ne 5 or Discriminant(f5) eq 0 then
            continue;
        end if;

        fI, L := IntegralModelPolynomial(f5);
        if Discriminant(fI) eq 0 then
            continue;
        end if;

        beta_ind := (2-r)/(4*(r-1));
        C := HyperellipticCurve(fI);
        J := Jacobian(C);
        Tind := J![X-beta_ind, Q!0];
        exact_tests +:= 1;
        if IsDivisibleBy(Tind, 2) then
            simple, pcert, Lp := IrreducibleFrobeniusCertificate(f5);
            Append(~hits, <r, f5, simple, pcert, Lp>);
            print "HIT r", r, "simple", simple, "p", pcert, "L", Lp;
        end if;
    end for;

    if den mod 50 eq 0 then
        print "den", den, "total", total, "candidates", #candidates,
              "exact_tests", exact_tests, "hits", #hits;
    end if;
end for;

print "Done";
print "total", total;
print "candidates", #candidates, candidates;
for r in candidates do
    f5 := OddQuinticForR(r);
    print "candidate r", r, "degree", Degree(f5), "discriminant", Discriminant(f5);
end for;
print "exact_tests", exact_tests;
print "hits", #hits;
for H in hits do
    print H;
end for;
quit;
