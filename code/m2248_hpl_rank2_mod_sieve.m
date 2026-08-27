//////////////////////////////////////////////////////////////////////
//  Local modular sieve for the fixed-rho HPL rank-2 Legendre lattice.
//
//  We test residue classes
//
//      Q(a,b,T) = a*P0 + b*P1 + T,    T in E0[2],
//
//  on the Legendre model
//
//      E0: V^2 = U(U-1)(U-rho0^2).
//
//  A rational point in the fixed-rho elliptic search must satisfy,
//  modulo every good prime:
//
//      U is a nonzero square,
//      tau^2 = U lifts to the quartic,
//      C = (sigma^2-rho0^2)/(sigma^2-1) is a nonzero square,
//      F1,F2,F3,F4 are nonzero squares.
//
//  This script checks how restrictive those conditions are locally.
//////////////////////////////////////////////////////////////////////

if not assigned prime_bound then
    prime_bound := 101;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;

if not assigned max_report_classes then
    max_report_classes := 40;
elif Type(max_report_classes) eq MonStgElt then
    max_report_classes := StringToInteger(max_report_classes);
end if;

if not assigned max_crt_classes then
    max_crt_classes := 200000;
elif Type(max_crt_classes) eq MonStgElt then
    max_crt_classes := StringToInteger(max_crt_classes);
end if;

Q := Rationals();
Z := Integers();

rho0 := Q!58466134224 / Q!53109477625;
sigma0 := Q!719363573659505664 / Q!749082246897952705;
tau0 := Q!307598400 / Q!352612321;
dB := Q!72946054224 / Q!53109477625;
Y0 := dB*(tau0^2-1);
V0 := tau0*Y0;
qbar0 := tau0;

U0 := tau0^2;
U1 := Q!185261445034489541 / Q!75208331264762500;
V1 := Q!2019033491444667568388406032021 / Q!950637494391639082835996875000;

function CanReduce(q, p)
    return (Z!Denominator(q) mod p) ne 0;
end function;

function CanReduceAll(vals, p)
    for q in vals do
        if not CanReduce(q, p) then
            return false;
        end if;
    end for;
    return true;
end function;

function IsNonzeroSquareK(x)
    if x eq 0 then
        return false, Parent(x)!0;
    end if;
    K := Parent(x);
    for i in [0..#K-1] do
        root := K!i;
        if root^2 eq x then
            return true, root;
        end if;
    end for;
    return false, K!0;
end function;

function CoverValues(rho, sigma, tau)
    F1 := (1 + rho)*(1 + sigma)*(1 + tau);
    F2 := rho*(1 + rho)*(rho + sigma)*(rho + tau);
    F3 := sigma*(1 + sigma)*(rho + sigma)*(sigma + tau);
    F4 := tau*(1 + tau)*(rho + tau)*(sigma + tau);
    return F1, F2, F3, F4;
end function;

function CRT2(r1, m1, r2, m2)
    r1 := Z!r1; m1 := Z!m1; r2 := Z!r2; m2 := Z!m2;
    g := GCD(m1, m2);
    delta := r2 - r1;
    if delta mod g ne 0 then
        return false, 0, 0;
    end if;

    l := LCM(m1, m2);
    m1g := m1 div g;
    m2g := m2 div g;
    if m2g eq 1 then
        k := Z!0;
    else
        inv := InverseMod(m1g mod m2g, m2g);
        k := ((delta div g) * inv) mod m2g;
    end if;
    return true, (r1 + m1*k) mod l, l;
end function;

function Key3(a, b, t)
    return Sprintf("%o,%o,%o", a, b, t);
end function;

function GoodPrimeData(p)
    if p eq 2 then
        return false, _, _, _, _, _, _, _;
    end if;
    vals := [rho0, qbar0, U0, V0, U1, V1];
    if not CanReduceAll(vals, p) then
        return false, _, _, _, _, _, _, _;
    end if;

    K := GF(p);
    rho := K!rho0;
    qbar := K!qbar0;
    if rho eq 0 or rho^2 eq 1 then
        return false, _, _, _, _, _, _, _;
    end if;

    E := EllipticCurve([K!0, -(K!1 + rho^2), K!0, rho^2, K!0]);
    P0 := E![K!U0, K!V0, K!1];
    P1 := E![K!U1, K!V1, K!1];
    Tors := [
        <"O", E!0>,
        <"T0", E![K!0, K!0, K!1]>,
        <"T1", E![K!1, K!0, K!1]>,
        <"Trho", E![rho^2, K!0, K!1]>
    ];
    return true, K, E, rho, qbar, P0, P1, Tors;
end function;

function PrimeSurvivors(p)
    ok, K, E, rho, qbar, P0, P1, Tors := GoodPrimeData(p);
    if not ok then
        return false, 0, 0, 0, [], [];
    end if;

    n0 := Order(P0);
    n1 := Order(P1);

    total_classes := n0*n1*#Tors;
    u_classes := 0; u_branches := 0;
    c_classes := 0; c_branches := 0;
    full_classes := 0; full_branches := 0;
    c_survivors := [];
    full_survivors := [];

    for a in [0..n0-1] do
        for b in [0..n1-1] do
            for ti in [1..#Tors] do
                Qe := a*P0 + b*P1 + Tors[ti][2];
                if Qe eq E!0 then
                    continue;
                end if;

                U := Qe[1]/Qe[3];
                V := Qe[2]/Qe[3];
                okU, tau0p := IsNonzeroSquareK(U);
                if not okU then
                    continue;
                end if;

                class_u := false;
                class_c := false;
                class_full := false;
                for tau in [tau0p, -tau0p] do
                    if tau eq 0 or tau^2 eq 1 then
                        continue;
                    end if;
                    class_u := true;
                    u_branches +:= 1;

                    sigma := qbar^2 * rho / tau;
                    if sigma eq 0 or sigma^2 eq 1 then
                        continue;
                    end if;
                    Cden := sigma^2 - 1;
                    if Cden eq 0 then
                        continue;
                    end if;
                    Cval := (sigma^2 - rho^2)/Cden;
                    okC, croot := IsNonzeroSquareK(Cval);
                    if not okC then
                        continue;
                    end if;
                    class_c := true;
                    c_branches +:= 1;

                    F1, F2, F3, F4 := CoverValues(rho, sigma, tau);
                    ok1, r1 := IsNonzeroSquareK(F1);
                    ok2, r2 := IsNonzeroSquareK(F2);
                    ok3, r3 := IsNonzeroSquareK(F3);
                    ok4, r4 := IsNonzeroSquareK(F4);
                    if ok1 and ok2 and ok3 and ok4 then
                        class_full := true;
                        full_branches +:= 1;
                    end if;
                end for;

                if class_u then
                    u_classes +:= 1;
                end if;
                if class_c then
                    c_classes +:= 1;
                    Append(~c_survivors, <a,b,ti>);
                end if;
                if class_full then
                    full_classes +:= 1;
                    Append(~full_survivors, <a,b,ti>);
                end if;
            end for;
        end for;
    end for;

    print "PRIME", p,
          "E_order", #E,
          "ordP0", n0,
          "ordP1", n1,
          "total_classes", total_classes;
    print "COUNTS", p,
          "u_classes", u_classes,
          "u_branches", u_branches,
          "c_classes", c_classes,
          "c_branches", c_branches,
          "full_classes", full_classes,
          "full_branches", full_branches;

    if #full_survivors le max_report_classes then
        print "FULL_SURVIVORS", p, full_survivors;
    else
        print "FULL_SURVIVORS", p, "count", #full_survivors;
    end if;

    return true, n0, n1, total_classes, c_survivors, full_survivors;
end function;

print "HPL fixed-rho rank-2 modular sieve";
print "prime_bound", prime_bound,
      "max_report_classes", max_report_classes,
      "max_crt_classes", max_crt_classes;
print "rho0", rho0;
print "qbar0", qbar0;
print "P0_UV", U0, V0;
print "P1_UV", U1, V1;

combined := [];
combined_seen := {};
combined_init := false;
mod_a := 1;
mod_b := 1;
combined_truncated := false;
used_primes := [];

for p in PrimesUpTo(prime_bound) do
    ok, n0, n1, total_classes, c_survivors, full_survivors := PrimeSurvivors(p);
    if not ok then
        print "SKIP_PRIME", p;
        continue;
    end if;
    hpl_good := false;
    for s in full_survivors do
        if s[2] eq 0 and s[3] eq 1 and (s[1] eq (1 mod n0) or s[1] eq ((n0-1) mod n0)) then
            hpl_good := true;
            break;
        end if;
    end for;
    if not hpl_good then
        print "SKIP_CRT_HPL_BAD", p, "known_hpl_class_not_full_nonzero";
        continue;
    end if;

    Append(~used_primes, p);

    if not combined_init then
        mod_a := n0;
        mod_b := n1;
        for s in full_survivors do
            key := Key3(s[1], s[2], s[3]);
            if key notin combined_seen then
                Include(~combined_seen, key);
                Append(~combined, s);
            end if;
        end for;
        combined_init := true;
    elif not combined_truncated then
        next_combined := [];
        next_seen := {};
        next_mod_a := LCM(mod_a, n0);
        next_mod_b := LCM(mod_b, n1);

        for old in combined do
            for s in full_survivors do
                if old[3] ne s[3] then
                    continue;
                end if;
                oka, ra, la := CRT2(old[1], mod_a, s[1], n0);
                if not oka then
                    continue;
                end if;
                okb, rb, lb := CRT2(old[2], mod_b, s[2], n1);
                if not okb then
                    continue;
                end if;
                key := Key3(ra, rb, old[3]);
                if key notin next_seen then
                    Include(~next_seen, key);
                    Append(~next_combined, <ra, rb, old[3]>);
                    if #next_combined gt max_crt_classes then
                        combined_truncated := true;
                        break old;
                    end if;
                end if;
            end for;
        end for;

        combined := next_combined;
        combined_seen := next_seen;
        mod_a := next_mod_a;
        mod_b := next_mod_b;
    end if;

    print "CRT_AFTER", p,
          "mod_a", mod_a,
          "mod_b", mod_b,
          "survivors", #combined,
          "truncated", combined_truncated;
    if #combined le max_report_classes and not combined_truncated then
        print "CRT_SURVIVORS", combined;
    end if;
end for;

print "DONE";
print "used_primes", used_primes;
print "final_mod_a", mod_a;
print "final_mod_b", mod_b;
print "final_full_survivors", #combined;
print "final_truncated", combined_truncated;
if #combined le max_report_classes and not combined_truncated then
    print "final_survivors", combined;
end if;
