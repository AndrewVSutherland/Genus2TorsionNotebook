//////////////////////////////////////////////////////////////////////
//  Component-wise p=5 boundary analysis for the halved-contact
//  M(2,12)+extra-3 route toward simple [3,12].
//
//  The affine M(2,12) chart is
//      m = (1-z^2)/(4*(r+1)).
//
//  On good p=5 affine points no finite group is compatible with [3,12].
//  Therefore rational examples must reduce to the p=5 boundary.  This
//  script factors the affine bad divisor, classifies P^1 x P^1 residue
//  classes at p=5, and then runs the same prime-to-5 residue filters
//  component by component before exact torsion checks.
//
//  Typical run:
//      magma -b height:=40 prime_bound:=43 \
//          code/contact6_m36_halveD_p5_boundary_analysis.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned height then
    height := 30;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned prime_bound then
    prime_bound := 43;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if not assigned max_exact then
    max_exact := 300;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned max_print then
    max_print := 80;
elif Type(max_print) eq MonStgElt then
    max_print := StringToInteger(max_print);
end if;
if not assigned progress_interval then
    progress_interval := 500000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals();
Z := Integers();
Qx<x> := PolynomialRing(Q);
F5 := GF(5);
R5<zz,rr> := PolynomialRing(F5, 2);

Hhigh := zz^6*rr^4 + 4*zz^6*rr^3 + 4*zz^6*rr^2
       + 3*zz^4*rr^4 + 3*zz^4*rr^3 + 4*zz^4*rr^2
       + 4*zz^4*rr + zz^2*rr^4 + 3*zz^2*rr^3
       + 2*zz^2*rr^2 + 4*zz^2*rr + zz^2 + 4*rr^2 + 3;

function JoinLabels(labels)
    if #labels eq 0 then
        return "";
    end if;
    ss := Sort(labels);
    out := ss[1];
    for i in [2..#ss] do
        out cat:= "+" cat ss[i];
    end for;
    return out;
end function;

procedure Increment(~A, key, amount)
    if IsDefined(A, key) then
        A[key] +:= amount;
    else
        A[key] := amount;
    end if;
end procedure;

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

function P1Residue(q, p)
    num := Numerator(q);
    den := Denominator(q);
    if (den mod p) eq 0 then
        return true, 0;
    end if;
    K := GF(p);
    return false, Z!(K!num / K!den);
end function;

function P1ValueString(is_inf, val)
    if is_inf then
        return "inf";
    end if;
    return IntegerToString(val);
end function;

function P5Labels(z_inf, zval, r_inf, rval)
    labels := [];
    if z_inf then
        Append(~labels, "Zinf");
    else
        z5 := F5!zval;
        if z5 eq 0 then
            Append(~labels, "Z0");
        end if;
        if z5 eq 1 then
            Append(~labels, "Zplus");
        end if;
        if z5 eq -1 then
            Append(~labels, "Zminus");
        end if;
    end if;

    if r_inf then
        Append(~labels, "Rinf");
    else
        r5 := F5!rval;
        if r5 eq 0 then
            Append(~labels, "R0");
        end if;
        if r5 eq -1 then
            Append(~labels, "Rminus1");
        end if;
    end if;

    if not z_inf and not r_inf then
        hval := Evaluate(Hhigh, <F5!zval, F5!rval>);
        if hval eq 0 then
            Append(~labels, "H");
        end if;
    end if;

    if #labels eq 0 then
        Append(~labels, "open_good");
    end if;
    return labels, JoinLabels(labels);
end function;

function P5Class(z, r)
    zinf, zval := P1Residue(z, 5);
    rinf, rval := P1Residue(r, 5);
    labels, sig := P5Labels(zinf, zval, rinf, rval);
    key := "z=" cat P1ValueString(zinf, zval) cat ",r=" cat
           P1ValueString(rinf, rval);
    return key, labels, sig, zinf, zval, rinf, rval;
end function;

function M212PolynomialFromZR(z, r)
    if r eq -1 or z^2 eq 1 then
        return false, Qx!0, Q!0;
    end if;
    m := (1-z^2)/(4*(r+1));
    if m eq 0 then
        return false, Qx!0, Q!0;
    end if;
    T := m*x^2 - x + r;
    h := (x-r)*(T+1);
    W := h^2 + 4*m*x^2*T*(T+1);
    return true, W, m;
end function;

function M212PolynomialFinite(K, z, r)
    P<X> := PolynomialRing(K);
    if r eq -K!1 or z^2 eq K!1 then
        return false, P!0;
    end if;
    m := (1-z^2)/(4*(r+1));
    if m eq 0 then
        return false, P!0;
    end if;
    T := m*X^2 - X + r;
    h := (X-r)*(T+1);
    W := h^2 + 4*m*X^2*T*(T+1);
    return true, W;
end function;

function GoodHyperellipticPolynomial(f)
    return Degree(f) in {5,6} and Discriminant(f) ne 0;
end function;

function Has312Invariants(invs)
    return #[n : n in invs | (Z!n) mod 3 eq 0] ge 2
       and #[n : n in invs | (Z!n) mod 12 eq 0] ge 1;
end function;

function InvariantsFinitePolynomial(f)
    C := HyperellipticCurve(f);
    A, phi := AbelianGroup(Jacobian(C));
    return Invariants(A);
end function;

function PairKey(a, b)
    return <Z!a, Z!b>;
end function;

function ResidueTables(p)
    K := GF(p);
    allowed := {};
    bad := {};
    for z in K do
        for r in K do
            key := PairKey(Z!z, Z!r);
            ok, W := M212PolynomialFinite(K, z, r);
            if not ok or not GoodHyperellipticPolynomial(W) then
                Include(~bad, key);
                continue;
            end if;
            invs := InvariantsFinitePolynomial(W);
            if Has312Invariants(invs) then
                Include(~allowed, key);
            end if;
        end for;
    end for;
    return allowed, bad;
end function;

function ResidueOfRational(q, p)
    num := Numerator(q);
    den := Denominator(q);
    if (den mod p) eq 0 then
        return false, 0;
    end if;
    K := GF(p);
    return true, Z!(K!num / K!den);
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return Qx!(L^2*f), L;
end function;

function SimpleCertificate(f)
    C := HyperellipticCurve(f);
    for p in [5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
        try
            fp := ChangeRing(f, GF(p));
            if not GoodHyperellipticPolynomial(fp) then
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
    return false, 0, Qx!0;
end function;

function TorsionOrder(invs)
    if #invs eq 0 then
        return 1;
    end if;
    return &*invs;
end function;

// Print the affine bad divisor factorization.
P5x<X> := PolynomialRing(R5);
d := 4*(rr+1);
n := 1-zz^2;
Tnum := n*X^2 - d*X + d*rr;
Tp1num := Tnum + d;
Wclear := d*(X-rr)^2*Tp1num^2 + 4*n*X^2*Tnum*Tp1num;
disc := Discriminant(Wclear);

print "Contact-6/M(2,12)+extra3 p=5 boundary analysis";
print "height", height, "prime_bound", prime_bound, "max_exact", max_exact;
print "Affine cleared sextic discriminant factorization over F_5:";
print Factorization(disc);
print "";

print "P1xP1 residue classes over F_5";
class_counts := AssociativeArray();
component_counts := AssociativeArray();
for zi in [-1..4] do
    zinf := zi eq -1;
    zval := zi;
    for ri in [-1..4] do
        rinf := ri eq -1;
        rval := ri;
        labels, sig := P5Labels(zinf, zval, rinf, rval);
        Increment(~class_counts, sig, 1);
        for lab in labels do
            Increment(~component_counts, lab, 1);
        end for;
        print "  class", "z", P1ValueString(zinf,zval),
              "r", P1ValueString(rinf,rval), "sig", sig;
    end for;
end for;
print "signature class counts";
for k in Sort([k : k in Keys(class_counts)]) do
    print " ", k, class_counts[k];
end for;
print "component class incidences";
for k in Sort([k : k in Keys(component_counts)]) do
    print " ", k, component_counts[k];
end for;
print "";

primes := [p : p in PrimesUpTo(prime_bound) | p notin {2,3,5}];
allowed_by_p := AssociativeArray();
bad_by_p := AssociativeArray();
print "Building prime-to-5 residue tables", primes;
for p in primes do
    allowed, bad := ResidueTables(p);
    allowed_by_p[p] := allowed;
    bad_by_p[p] := bad;
    print " prime", p, "allowed_312", #allowed, "bad", #bad;
end for;

params := RationalParametersOfHeight(height);
checked := 0;
smooth := 0;
unique_ar := {};
open_killed_p5 := 0;
boundary_smooth_by_sig := AssociativeArray();
residue_survivors_by_sig := AssociativeArray();
simple_by_sig := AssociativeArray();
exact_by_sig := AssociativeArray();
torsion_by_sig := AssociativeArray();
class_survivors := AssociativeArray();
hits := [];
residue_survivors := 0;
simple_survivors := 0;
exact_tests := 0;
printed := 0;

for z in params do
    for r in params do
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "progress", checked, "smooth", smooth,
                  "residue_survivors", residue_survivors,
                  "simple", simple_survivors, "exact", exact_tests,
                  "hits", #hits;
        end if;

        ok, W, m := M212PolynomialFromZR(z, r);
        if not ok or not GoodHyperellipticPolynomial(W) then
            continue;
        end if;

        ar_key := Sprint(<m,r>);
        if ar_key in unique_ar then
            continue;
        end if;
        Include(~unique_ar, ar_key);
        smooth +:= 1;

        p5class, labels, sig, zinf, zval, rinf, rval := P5Class(z, r);
        if sig eq "open_good" then
            open_killed_p5 +:= 1;
            continue;
        end if;
        Increment(~boundary_smooth_by_sig, sig, 1);

        killed := false;
        for p in primes do
            okz, rz := ResidueOfRational(z, p);
            okr, rr0 := ResidueOfRational(r, p);
            if not okz or not okr then
                continue;
            end if;
            key := PairKey(rz, rr0);
            if key in bad_by_p[p] then
                continue;
            end if;
            if key notin allowed_by_p[p] then
                killed := true;
                break;
            end if;
        end for;
        if killed then
            continue;
        end if;

        residue_survivors +:= 1;
        Increment(~residue_survivors_by_sig, sig, 1);
        Increment(~class_survivors, p5class cat " | " cat sig, 1);

        WI, scale := IntegralModelPolynomial(W);
        simple, pcert, Lp := SimpleCertificate(WI);
        if not simple then
            if printed lt max_print then
                print "SURVIVOR_NONSIMPLE_UNKNOWN", "p5", p5class, "sig", sig,
                      "z", z, "r", r, "m", m;
                printed +:= 1;
            end if;
            continue;
        end if;

        simple_survivors +:= 1;
        Increment(~simple_by_sig, sig, 1);
        if exact_tests ge max_exact then
            continue;
        end if;

        C := HyperellipticCurve(WI);
        G, phi := TorsionSubgroup(Jacobian(C));
        invs := Invariants(G);
        exact_tests +:= 1;
        Increment(~exact_by_sig, sig, 1);
        tors_key := sig cat " | " cat Sprint(invs);
        Increment(~torsion_by_sig, tors_key, 1);

        if printed lt max_print then
            print "EXACT", "p5", p5class, "sig", sig,
                  "z", z, "r", r, "m", m,
                  "torsion", invs, "pcert", pcert;
            printed +:= 1;
        end if;

        if Has312Invariants(invs) then
            Append(~hits, <z,r,m,sig,p5class,invs,WI,pcert,Lp>);
            print "HIT312", "p5", p5class, "sig", sig,
                  "z", z, "r", r, "m", m, "torsion", invs;
            print " W", WI;
        end if;
    end for;
end for;

print "";
print "DONE p=5 component boundary analysis";
print "checked", checked, "smooth_unique_ar", smooth, "unique_ar", #unique_ar,
      "open_killed_p5", open_killed_p5,
      "residue_survivors", residue_survivors,
      "simple_survivors", simple_survivors,
      "exact_tests", exact_tests, "hits", #hits;
print "boundary_smooth_by_sig";
for k in Sort([k : k in Keys(boundary_smooth_by_sig)]) do
    print " ", k, boundary_smooth_by_sig[k];
end for;
print "residue_survivors_by_sig";
for k in Sort([k : k in Keys(residue_survivors_by_sig)]) do
    print " ", k, residue_survivors_by_sig[k];
end for;
print "simple_by_sig";
for k in Sort([k : k in Keys(simple_by_sig)]) do
    print " ", k, simple_by_sig[k];
end for;
print "exact_by_sig";
for k in Sort([k : k in Keys(exact_by_sig)]) do
    print " ", k, exact_by_sig[k];
end for;
print "torsion_by_sig";
for k in Sort([k : k in Keys(torsion_by_sig)]) do
    print " ", k, torsion_by_sig[k];
end for;
print "class_survivors";
for k in Sort([k : k in Keys(class_survivors)]) do
    print " ", k, class_survivors[k];
end for;
for H in hits do
    print "H", H;
end for;

quit;
