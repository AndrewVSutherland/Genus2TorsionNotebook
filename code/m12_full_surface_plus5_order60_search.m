//////////////////////////////////////////////////////////////////////
//  Full two-parameter M(12) split-Weierstrass cover plus 5-torsion.
//
//  The M(12) model is
//      y^2 + (x-r)(T+1)y = a*x^2*T*(T+1),
//      T = a*x^2-x+r.
//  To make the distinguished Weierstrass point rational, put
//      a = (1-z^2)/(4*(r+1)).
//  Then w = 2*(r+1)/(1-z) is a root of T+1.  Moving w to infinity
//  gives an odd quintic with a visibly rational point of order 12.
//  An additional rational 5-torsion point therefore gives exponent 60.
//
//  For every good affine residue (r,z) modulo p != 2,5, the script
//  keeps it only if 5 divides #J(F_p).  Boundary residues and parameters
//  whose denominator is divisible by p are passed, so the mask is a
//  rigorous necessary filter rather than a heuristic exclusion.
//
//  Typical runs (from torsion_jac):
//      magma -b mode:="finite" prime_bound:=19 \
//          code/m12_full_surface_plus5_order60_search.m
//      magma -b mode:="search" height:=20 prime_bound:=19 max_exact:=100 \
//          code/m12_full_surface_plus5_order60_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then
    mode := "search";
end if;
if not assigned height then
    height := 15;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned prime_bound then
    prime_bound := 19;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if not assigned max_exact then
    max_exact := 100;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;
if not assigned max_print then
    max_print := 30;
elif Type(max_print) eq MonStgElt then
    max_print := StringToInteger(max_print);
end if;
if not assigned progress_interval then
    progress_interval := 100000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;
if not assigned simplicity_bound then
    simplicity_bound := 97;
elif Type(simplicity_bound) eq MonStgElt then
    simplicity_bound := StringToInteger(simplicity_bound);
end if;
if not assigned certify_all_exact then
    certify_all_exact := false;
elif Type(certify_all_exact) eq MonStgElt then
    certify_all_exact := certify_all_exact eq "true";
end if;
if not assigned MemGB then
    MemGB := 8;
elif Type(MemGB) eq MonStgElt then
    MemGB := StringToInteger(MemGB);
end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);
PX<X> := PolynomialRing(Q);
PT<U> := PolynomialRing(Q);

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num, den) ne 1 then
                continue;
            end if;
            q := Q!num/Q!den;
            key := Sprint(q);
            if key notin seen then
                Include(~seen, key);
                Append(~vals, q);
            end if;
        end for;
    end for;
    return vals;
end function;

function OddQuintic(r, z)
    if r eq -1 or z in {Q!-1, Q!1} then
        return false, PX!0, Q!0, Q!0, Q!0;
    end if;
    a := (1-z^2)/(4*(r+1));
    if a eq 0 then
        return false, PX!0, a, Q!0, Q!0;
    end if;
    T := a*x^2-x+r;
    h := (x-r)*(T+1);
    W := h^2+4*a*x^2*T*(T+1);
    w := 2*(r+1)/(1-z);
    if Evaluate(T+1, w) ne 0 then
        return false, PX!0, a, w, Q!0;
    end if;
    f5 := PX!0;
    for i in [0..Degree(W)] do
        for j in [0..i] do
            f5 +:= Coefficient(W,i)*Binomial(i,j)*w^(i-j)*X^(6-j);
        end for;
    end for;
    Xp := -1/w;
    Yp := Evaluate(h,0)*Xp^3;
    return true, f5, a, Xp, Yp;
end function;

function OddQuinticFinite(rr, zz, F)
    PF<xx> := PolynomialRing(F);
    if rr eq -1 or zz^2 eq 1 then
        return false, PF!0;
    end if;
    a := (1-zz^2)/(F!4*(rr+1));
    if a eq 0 then
        return false, PF!0;
    end if;
    T := a*xx^2-xx+rr;
    h := (xx-rr)*(T+1);
    W := h^2+4*a*xx^2*T*(T+1);
    w := 2*(rr+1)/(1-zz);
    if Evaluate(T+1,w) ne 0 then
        return false, PF!0;
    end if;
    f5 := PF!0;
    for i in [0..Degree(W)] do
        for j in [0..i] do
            f5 +:= Coefficient(W,i)*Binomial(i,j)*w^(i-j)*xx^(6-j);
        end for;
    end for;
    return true, f5;
end function;

function IsGoodQuintic(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function ResidueKey(a, b, p)
    return a*p+b;
end function;

function FiniteMask(p)
    F := GF(p);
    allowed := {};
    bad := {};
    good := 0;
    killed := 0;
    for ri in [0..p-1] do
        for zi in [0..p-1] do
            key := ResidueKey(ri,zi,p);
            ok, f5 := OddQuinticFinite(F!ri,F!zi,F);
            if not ok or not IsGoodQuintic(f5) then
                Include(~bad,key);
                continue;
            end if;
            good +:= 1;
            try
                C := HyperellipticCurve(f5);
                N := Z!Evaluate(LPolynomial(C),1);
                if N mod 5 eq 0 then
                    Include(~allowed,key);
                else
                    killed +:= 1;
                end if;
            catch e
                Include(~bad,key);
                good -:= 1;
            end try;
        end for;
    end for;
    return allowed,bad,good,killed;
end function;

function ResidueOfRational(q,p)
    if Denominator(q) mod p eq 0 then
        return false,0;
    end if;
    F := GF(p);
    return true,Z!(F!Numerator(q)/F!Denominator(q));
end function;

function PassesMasks(r,z,plist,allowed,badres)
    good_checked := 0;
    boundary_primes := [];
    for p in plist do
        okr,rr := ResidueOfRational(r,p);
        okz,zz := ResidueOfRational(z,p);
        if not okr or not okz then
            Append(~boundary_primes,p);
            continue;
        end if;
        key := ResidueKey(rr,zz,p);
        if key in badres[p] then
            Append(~boundary_primes,p);
            continue;
        end if;
        good_checked +:= 1;
        if key notin allowed[p] then
            return false,p,good_checked,boundary_primes;
        end if;
    end for;
    return true,0,good_checked,boundary_primes;
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L,Denominator(Coefficient(f,i)));
    end for;
    return PX!(L^2*f),L;
end function;

function TorsionExponent(invs)
    e := 1;
    for n in invs do
        e := LCM(e,n);
    end for;
    return e;
end function;

function PointOfOrder60(G,phi)
    if #Invariants(G) eq 0 then
        return false,Codomain(phi)!0,G!0;
    end if;
    g := G!0;
    for i in [1..Ngens(G)] do
        g +:= G.i;
    end for;
    og := Order(g);
    if og mod 60 ne 0 then
        return false,Codomain(phi)!0,G!0;
    end if;
    g60 := (og div 60)*g;
    return Order(g60) eq 60,phi(g60),g60;
end function;

function FrobeniusPolynomial(C,p)
    ef := EulerFactor(C,p);
    d := Degree(ef);
    return &+[Q!Coefficient(ef,i)*U^(d-i) : i in [0..d]];
end function;

function FullSimplicityCertificate(C,plist)
    // The D4 test is a direct geometric-simplicity certificate.
    for p in plist do
        try
            Phi := FrobeniusPolynomial(C,p);
            fac := Factorization(Phi);
            if #fac ne 1 or fac[1][2] ne 1 or Degree(fac[1][1]) ne 4 then
                continue;
            end if;
            Gal := GaloisGroup(Phi);
            desc := "";
            try
                desc := TransitiveGroupDescription(Gal);
            catch e2
                desc := "unknown";
            end try;
            if Order(Gal) eq 8 and desc eq "D(4)" then
                return true,"D4",p,Phi;
            end if;
        catch e
            continue;
        end try;
    end for;

    // Fallback: certify an absolutely simple reduction by checking that
    // no Frobenius root power through the surface bound 12 drops degree.
    for p in plist do
        try
            Phi := FrobeniusPolynomial(C,p);
            fac := Factorization(Phi);
            if Degree(Phi) ne 4 or #fac ne 1 or fac[1][2] ne 1 or
               Degree(fac[1][1]) ne 4 then
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
                return true,"root_power",p,Phi;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false,"none",0,PT!0;
end function;

plist := [p : p in PrimesUpTo(prime_bound) | p notin {2,5}];
allowed := AssociativeArray();
badres := AssociativeArray();

print "M(12) full two-parameter split cover plus 5/order-60 search";
print "mode",mode,"height",height,"prime_bound",prime_bound,
      "primes",plist,"max_exact",max_exact;
print "Precomputing rigorous affine residue masks";
for p in plist do
    good_allowed,bad,good,killed := FiniteMask(p);
    allowed[p] := good_allowed;
    badres[p] := bad;
    print "MASK",p,"good",good,"allowed5",#good_allowed,
          "bad_boundary",#bad,"killed",killed;
end for;

if mode eq "finite" then
    print "DONE finite masks";
    quit;
end if;

params := RationalParametersOfHeight(height);
checked := 0;
mask_survivors := 0;
smooth_survivors := 0;
unique_survivors := 0;
exact_tests := 0;
order12_verified := 0;
hits := 0;
simple_hits := 0;
printed := 0;
exact_printed := 0;
kill_counts := AssociativeArray();
seen_curves := {};
candidates := [];

for r in params do
    for z in params do
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "PROGRESS",checked,"mask_survivors",mask_survivors,
                  "smooth",smooth_survivors,"unique",unique_survivors;
        end if;
        if r eq -1 or z in {Q!-1,Q!1} then
            continue;
        end if;
        pass,killp,good_checked,boundary_primes :=
            PassesMasks(r,z,plist,allowed,badres);
        if not pass then
            if IsDefined(kill_counts,killp) then
                kill_counts[killp] +:= 1;
            else
                kill_counts[killp] := 1;
            end if;
            continue;
        end if;
        mask_survivors +:= 1;

        ok,f5,a,Xp,Yp := OddQuintic(r,z);
        if not ok or not IsGoodQuintic(f5) then
            continue;
        end if;
        smooth_survivors +:= 1;
        // The signs z and -z give the same a and hence the same curve;
        // they merely select the other root of T+1 as the base point.
        curve_key := Sprint(<r,a>);
        if curve_key in seen_curves then
            continue;
        end if;
        Include(~seen_curves,curve_key);
        unique_survivors +:= 1;
        Append(~candidates,
               <good_checked,r,z,a,Xp,Yp,f5,boundary_primes>);
        if printed lt max_print then
            print "SURVIVOR","r",r,"z",z,"a",a,
                  "good_checked",good_checked,
                  "boundary_primes",boundary_primes;
            printed +:= 1;
        end if;
    end for;
end for;

print "CANDIDATE_POOL","parameter_survivors",smooth_survivors,
      "unique_curves",unique_survivors;
print "Exact checks are prioritized by the number of genuinely good masks";

// Boundary-rich candidates are far more common.  Test the candidates with
// the most genuine good-reduction constraints first, without discarding any
// candidate from the pool.
for priority in Reverse([0..#plist]) do
    for H in candidates do
        if exact_tests ge max_exact or hits ge max_hits then
            break priority;
        end if;
        good_checked := H[1];
        if good_checked ne priority then
            continue;
        end if;
        r := H[2];
        z := H[3];
        a := H[4];
        Xp := H[5];
        Yp := H[6];
        f5 := H[7];
        boundary_primes := H[8];
        fI,L := IntegralModelPolynomial(f5);
        if not IsGoodQuintic(fI) then
            continue;
        end if;
        exact_tests +:= 1;
        try
            C := HyperellipticCurve(fI);
            J := Jacobian(C);
            D12 := J![X-Xp,L*Yp];
            od := Order(D12);
            if od ne 12 then
                print "WARNING_D12_ORDER","r",r,"z",z,"order",od;
                continue;
            end if;
            order12_verified +:= 1;

            G,phi := TorsionSubgroup(J);
            invs := Invariants(G);
            exp := TorsionExponent(invs);
            if exact_printed lt max_print then
                print "EXACT","r",r,"z",z,"a",a,
                      "torsion",invs,"exponent",exp,
                      "D12_order",od,"good_checked",good_checked,
                      "boundary_primes",boundary_primes;
                exact_printed +:= 1;
            end if;

            if certify_all_exact then
                splist0 := [p : p in PrimesUpTo(simplicity_bound) | p ge 3];
                simple0,stype0,sp0,Phi0 :=
                    FullSimplicityCertificate(C,splist0);
                print "EXACT_SIMPLICITY","r",r,"z",z,
                      "simple",simple0,"type",stype0,"prime",sp0;
                if simple0 then
                    print "  Frobenius =",Phi0;
                end if;
            end if;

            if exp mod 60 eq 0 then
                hits +:= 1;
                ok60,P60,g60 := PointOfOrder60(G,phi);
                splist := [p : p in PrimesUpTo(simplicity_bound) | p ge 3];
                simple,stype,sp,Phi := FullSimplicityCertificate(C,splist);
                if simple then
                    simple_hits +:= 1;
                end if;
                print "ORDER60_HIT","r",r,"z",z,"a",a,
                      "torsion",invs,"exponent",exp,
                      "explicit60",ok60,"simple",simple,
                      "simple_type",stype,"simple_prime",sp;
                print "  f_integral =",fI;
                print "  D12 =",D12;
                if ok60 then
                    print "  P60 =",P60;
                    print "  abstract_g60 =",g60;
                end if;
                if simple then
                    print "  Frobenius =",Phi;
                end if;
            end if;
        catch e
            print "EXACT_ERROR","r",r,"z",z,"message",e`Object;
            continue;
        end try;
    end for;
end for;

print "DONE M12 full-surface plus 5/order 60";
print "checked",checked,"mask_survivors",mask_survivors,
      "smooth_survivors",smooth_survivors,
      "unique_survivors",unique_survivors,"exact_tests",exact_tests,
      "order12_verified",order12_verified,"hits",hits,
      "simple_hits",simple_hits;
print "kill_counts";
for p in Sort([k : k in Keys(kill_counts)]) do
    print p,kill_counts[p];
end for;
quit;
