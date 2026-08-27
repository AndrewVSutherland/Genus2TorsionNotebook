//////////////////////////////////////////////////////////////////////
// Bounded rational-base scan for the degree-2 [5,5] contact family.
//
// Typical run:
//   magma -b integer_height:=30 rational_height:=6 \
//       code/z5x5_degree2_new_base_exact_scan.m
//
// Phase 1 covers all integral triples in the declared box.  Phase 2
// covers all reduced rational triples of coordinate height at most the
// declared bound, omitting triples already covered in phase 1.  Candidates
// must have actual finite 5-rank at least two at two good primes among
// p=3,...,43 and a geometrically-simple Frobenius certificate at a good
// small prime.  Singular or undefined reductions of the displayed equation
// are skipped, never treated as obstructions or intrinsic bad reduction.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned integer_height then
    integer_height := 30;
elif Type(integer_height) eq MonStgElt then
    integer_height := StringToInteger(integer_height);
end if;
if not assigned rational_height then
    rational_height := 6;
elif Type(rational_height) eq MonStgElt then
    rational_height := StringToInteger(rational_height);
end if;

Q := Rationals();
Z := Integers();
Qx<x> := PolynomialRing(Q);
screen_primes := [3, 7, 11, 13, 17];
fallback_screen_primes := [19, 23, 29, 31, 37, 41, 43];
certificate_primes := [3, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43];

function FiveRank(invs)
    return #[n : n in invs | (Z!n) mod 5 eq 0];
end function;

function BasePolynomial(a, b, k)
    return (1+a*x+b*x^2)^2-k*x^5;
end function;

function IntegralSquareModel(f)
    denominator_lcm := 1;
    for i in [0..Degree(f)] do
        denominator_lcm := LCM(denominator_lcm,
                               Denominator(Coefficient(f, i)));
    end for;
    return Qx!(denominator_lcm^2*f), denominator_lcm;
end function;

function ResidueKey(base, p)
    K := GF(p);
    residues := [];
    for q in base do
        if Denominator(q) mod p eq 0 then
            return false, <>;
        end if;
        Append(~residues, Z!(K!Numerator(q)/K!Denominator(q)));
    end for;
    return true, <residues[1], residues[2], residues[3]>;
end function;

function FiniteTable(p)
    K := GF(p);
    P<X> := PolynomialRing(K);
    finite_table := AssociativeArray();
    tested := 0;
    smooth := 0;
    rank_counts := AssociativeArray();
    for a in K do
        for b in K do
            for k in K do
                if k eq 0 then
                    continue;
                end if;
                tested +:= 1;
                f := (1+a*X+b*X^2)^2-k*X^5;
                if Degree(f) ne 5 or Discriminant(f) eq 0 then
                    continue;
                end if;
                smooth +:= 1;
                A, phi := AbelianGroup(Jacobian(HyperellipticCurve(f)));
                invs := Invariants(A);
                rank := FiveRank(invs);
                if IsDefined(rank_counts, rank) then
                    rank_counts[rank] +:= 1;
                else
                    rank_counts[rank] := 1;
                end if;
                finite_table[Sprint(<Z!a, Z!b, Z!k>)] := <rank, invs, #A>;
            end for;
        end for;
    end for;
    printf "SCREEN_SUMMARY p=%o tested=%o smooth=%o allowed=%o rank_counts=%o\n",
           p, tested, smooth,
           #[key : key in Keys(finite_table) | finite_table[key][1] ge 2],
           Sort([<r, rank_counts[r]> : r in Keys(rank_counts)]);
    return finite_table;
end function;

function PowerTwelveIrreducible(L)
    ZW<W> := PolynomialRing(Z);
    ZWT<T> := PolynomialRing(ZW);
    LL := &+[ZW!Coefficient(L, i)*T^i : i in [0..Degree(L)]];
    transform := ZW!Resultant(LL, T^12-W);
    return IsIrreducible(transform), transform;
end function;

function DirectFiniteDatum(f, p)
    try
        fp := ChangeRing(f, GF(p));
        if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
            return false, <>;
        end if;
        A, phi := AbelianGroup(Jacobian(HyperellipticCurve(fp)));
        invs := Invariants(A);
        return true, <FiveRank(invs), invs, #A>;
    catch e
        return false, <>;
    end try;
end function;

function GeometricSimpleCertificate(f)
    for p in certificate_primes do
        try
            fp := ChangeRing(f, GF(p));
            if Degree(fp) ne 5 or Discriminant(fp) eq 0 then
                continue;
            end if;
            L := LPolynomial(HyperellipticCurve(fp));
            simple, transform := PowerTwelveIrreducible(L);
            if simple then
                return true, p, L, transform;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false, 0, Z!0, Z!0;
end function;

function RationalValues(B)
    values := [];
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
                Append(~values, q);
            end if;
        end for;
    end for;
    Sort(~values);
    return values;
end function;

finite_tables := AssociativeArray();
for p in screen_primes do
    finite_tables[p] := FiniteTable(p);
end for;

tested_bases := {};
phase_counts := AssociativeArray();
total_survivors := 0;
simple_survivors := 0;
exact_hits := [];
insufficient_good_primes := 0;

procedure TestBase(phase, base, ~tested_bases, ~phase_counts,
                   ~total_survivors, ~simple_survivors, ~exact_hits,
                   ~insufficient_good_primes)
    key := Sprint(base);
    if key in tested_bases then
        return;
    end if;
    Include(~tested_bases, key);
    counts := phase_counts[phase];
    counts[1] +:= 1;
    phase_counts[phase] := counts;
    if base[3] eq 0 then
        return;
    end if;
    f := BasePolynomial(base[1], base[2], base[3]);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return;
    end if;

    finite_data := [* *];
    good_rank_two_primes := 0;
    for p in screen_primes do
        defined, residue := ResidueKey(base, p);
        if not defined then
            Append(~finite_data, <p, "undefined_base_reduction">);
            continue;
        end if;
        residue_key := Sprint(residue);
        if not IsDefined(finite_tables[p], residue_key) then
            Append(~finite_data, <p, residue, "singular_displayed_reduction">);
            continue;
        end if;
        datum := finite_tables[p][residue_key];
        Append(~finite_data, <p, residue, datum>);
        if datum[1] lt 2 then
            return;
        end if;
        good_rank_two_primes +:= 1;
    end for;
    if good_rank_two_primes lt 2 then
        for p in fallback_screen_primes do
            good, datum := DirectFiniteDatum(f, p);
            if not good then
                Append(~finite_data,
                       <p, "singular_or_undefined_displayed_reduction">);
                continue;
            end if;
            Append(~finite_data, <p, datum>);
            if datum[1] lt 2 then
                return;
            end if;
            good_rank_two_primes +:= 1;
            if good_rank_two_primes ge 2 then
                break;
            end if;
        end for;
    end if;
    if good_rank_two_primes lt 2 then
        insufficient_good_primes +:= 1;
        printf "INSUFFICIENT_GOOD_PRIMES phase=%o base=%o finite=%o\n",
               phase, base, finite_data;
        return;
    end if;
    counts := phase_counts[phase];
    counts[2] +:= 1;
    phase_counts[phase] := counts;
    total_survivors +:= 1;

    simple, pcert, Lcert, transform := GeometricSimpleCertificate(f);
    if not simple then
        printf "NO_SIMPLE_CERT phase=%o base=%o finite=%o\n", phase, base, finite_data;
        return;
    end if;
    counts := phase_counts[phase];
    counts[3] +:= 1;
    phase_counts[phase] := counts;
    simple_survivors +:= 1;

    f_integral, y_scale := IntegralSquareModel(f);
    C := HyperellipticCurve(f_integral);
    G, phi := TorsionSubgroup(Jacobian(C));
    invs := Invariants(G);
    counts := phase_counts[phase];
    counts[4] +:= 1;
    phase_counts[phase] := counts;
    printf "EXACT_TEST phase=%o base=%o finite=%o simple_p=%o L=%o y_scale=%o torsion=%o\n",
           phase, base, finite_data, pcert, Lcert, y_scale, invs;
    if FiveRank(invs) ge 2 then
        Append(~exact_hits, <phase, base, f, invs, pcert, Lcert, transform>);
        printf "EXACT_55_HIT phase=%o base=%o f=%o torsion=%o simple_p=%o L=%o transform=%o\n",
               phase, base, f, invs, pcert, Lcert, transform;
    end if;
end procedure;

phase_counts["integer"] := [0, 0, 0, 0];
integer_values := [Q!n : n in [-integer_height..integer_height]];
for a in integer_values do
    for b in integer_values do
        for k in integer_values do
            TestBase("integer", <a, b, k>, ~tested_bases, ~phase_counts,
                     ~total_survivors, ~simple_survivors, ~exact_hits,
                     ~insufficient_good_primes);
        end for;
    end for;
end for;

phase_counts["rational"] := [0, 0, 0, 0];
rational_values := RationalValues(rational_height);
for a in rational_values do
    for b in rational_values do
        for k in rational_values do
            TestBase("rational", <a, b, k>, ~tested_bases, ~phase_counts,
                     ~total_survivors, ~simple_survivors, ~exact_hits,
                     ~insufficient_good_primes);
        end for;
    end for;
end for;

print "# counts are <enumerated unique triples, two-prime rank survivors, simple survivors, exact tests>";
print "PHASE_SUMMARY integer", phase_counts["integer"];
print "PHASE_SUMMARY rational", phase_counts["rational"];
print "TOTAL_UNIQUE_BASES", #tested_bases;
print "TOTAL_TWO_PRIME_SURVIVORS", total_survivors;
print "TOTAL_SIMPLE_SURVIVORS", simple_survivors;
print "INSUFFICIENT_GOOD_PRIME_COUNT", insufficient_good_primes;
print "EXACT_HIT_COUNT", #exact_hits;
for hit in exact_hits do
    print "HIT", hit;
end for;
