//////////////////////////////////////////////////////////////////////
// Strict certificates for three existing 2-primary torsion banks.
//
// Targets:
//   A(2,2,4,4): data/tor2244_bank.txt
//   M(2,2,2,8): data/tor2228_bank.txt
//   M(2,4,4):   regenerate the small elliptic-fiber-product pilot
//
// For every tested curve this script computes the exact rational torsion
// subgroup first.  Only an exact target match is then sent to the
// geometric-simplicity screen.  A witness is a good prime at which the
// Frobenius quartic is irreducible and
//
//       Degree(MinimalPolynomial(pi^n)) = 4,  n = 2,...,12.
//
// Typical use (from torsion_jac):
//
//   magma -b code/certify_frontier_existing_banks.m
//   magma -b target:=2244 max_bank_rows:=100 \
//       code/certify_frontier_existing_banks.m
//   magma -b target:=2244 start_bank_row:=26629 max_bank_rows:=1 \
//       code/certify_frontier_existing_banks.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned target then
    target := "all";
end if;
if Type(target) ne MonStgElt then
    target := IntegerToString(Integers()!target);
end if;

if not assigned max_bank_rows then
    max_bank_rows := 25;
elif Type(max_bank_rows) eq MonStgElt then
    max_bank_rows := StringToInteger(max_bank_rows);
end if;

if not assigned start_bank_row then
    start_bank_row := 1;
elif Type(start_bank_row) eq MonStgElt then
    start_bank_row := StringToInteger(start_bank_row);
end if;

if not assigned witnesses_per_curve then
    witnesses_per_curve := 2;
elif Type(witnesses_per_curve) eq MonStgElt then
    witnesses_per_curve := StringToInteger(witnesses_per_curve);
end if;

if not assigned m244_height then
    m244_height := 4;
elif Type(m244_height) eq MonStgElt then
    m244_height := StringToInteger(m244_height);
end if;

if not assigned m244_point_bound then
    m244_point_bound := 80;
elif Type(m244_point_bound) eq MonStgElt then
    m244_point_bound := StringToInteger(m244_point_bound);
end if;

if not assigned max_m244_curves then
    max_m244_curves := 25;
elif Type(max_m244_curves) eq MonStgElt then
    max_m244_curves := StringToInteger(max_m244_curves);
end if;

if not assigned diagnose_failures then
    diagnose_failures := false;
elif Type(diagnose_failures) eq MonStgElt then
    diagnose_failures := diagnose_failures in {"true","True","1","yes"};
end if;

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);

PrimeList := [
    3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
    73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,
    157,163,167,173,179,181,191,193,197,199
];

function ReadTupleFile(filename)
    rows := Split(Read(filename), "\n");
    tuples := [];
    for rawrow in rows do
        n := #rawrow;
        if n gt 0 and rawrow[n] eq "\r" then n -:= 1; end if;
        if n ge 2 and rawrow[1] eq "[" and rawrow[n] eq "]" then
            parts := Split(rawrow[2..n-1], ",");
            if #parts eq 4 then
                Append(~tuples, [StringToInteger(s) : s in parts]);
            end if;
        end if;
    end for;
    return tuples;
end function;

function SquareBranchModel(tup)
    a,b,c,d := Explode([Q!z : z in tup]);
    return x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2);
end function;

function CountCurve(fp)
    Fq := BaseRing(Parent(fp));
    count := 0;
    for xx in Fq do
        value := Evaluate(fp,xx);
        if value eq 0 then
            count +:= 1;
        elif IsSquare(value) then
            count +:= 2;
        end if;
    end for;
    if IsOdd(Degree(fp)) then
        count +:= 1;
    elif IsSquare(LeadingCoefficient(fp)) then
        count +:= 2;
    end if;
    return count;
end function;

function RootPowerAtPrime(f, pp)
    if Z!Denominator(LeadingCoefficient(f)) mod pp eq 0 then
        return false, P!0, [];
    end if;
    if Z!Numerator(LeadingCoefficient(f)) mod pp eq 0 then
        return false, P!0, [];
    end if;
    disc := Discriminant(f);
    if Z!Denominator(disc) mod pp eq 0 or Z!Numerator(disc) mod pp eq 0 then
        return false, P!0, [];
    end if;

    Fp := GF(pp); Pp := PolynomialRing(Fp);
    fp := Pp![Fp!c : c in Coefficients(f)];
    if Degree(fp) ne Degree(f) or not IsSquarefree(fp) then
        return false, P!0, [];
    end if;
    Fp2 := GF(pp^2); Pp2 := PolynomialRing(Fp2);
    fp2 := Pp2![Fp2!c : c in Coefficients(f)];

    a1 := pp + 1 - CountCurve(fp);
    a2 := (CountCurve(fp2) - pp^2 - 1 + a1^2) div 2;
    // Independent built-in cross-check of the point-count reconstruction.
    Lp := LPolynomial(HyperellipticCurve(fp));
    assert [Z!Coefficient(Lp,i) : i in [0..4]]
           eq [1,-a1,a2,-a1*pp,pp^2];
    R<T> := PolynomialRing(Q);
    chi := T^4 - a1*T^3 + a2*T^2 - a1*pp*T + pp^2;
    if not IsIrreducible(chi) then
        return false, P!chi, [];
    end if;
    K<pi> := NumberField(chi);
    degrees := [Degree(MinimalPolynomial(pi^n)) : n in [2..12]];
    return &and[d eq 4 : d in degrees], P!chi, degrees;
end function;

function RootPowerCertificates(f, how_many)
    out := [];
    for pp in PrimeList do
        ok, chi, degrees := RootPowerAtPrime(f,pp);
        if ok then
            Append(~out,<pp,chi,degrees>);
            if #out ge how_many then break; end if;
        end if;
    end for;
    return out;
end function;

function ExactTorsion(f)
    C := HyperellipticCurve(f);
    G, mp := TorsionSubgroup(Jacobian(C));
    return Invariants(G), #G;
end function;

function IntegralSquareTwistModel(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L,Denominator(Coefficient(f,i)));
    end for;
    // Multiplication by the rational square L^2 preserves the Q-isomorphism
    // class while giving TorsionSubgroup the integral model it requires.
    return P!(L^2*f);
end function;

procedure PrintSuccess(label, source, f, invs, certs)
    print "SUCCESS_GROUP", label;
    print "SOURCE", source;
    print "CURVE_Y2_EQUALS", f;
    Cs := SimplifiedModel(HyperellipticCurve(f));
    fs,hs := HyperellipticPolynomials(Cs);
    print "SIMPLIFIED_HYPERELLIPTIC_POLYNOMIALS", fs, hs;
    print "EXACT_TORSION", invs;
    print "EXACT_TORSION_ORDER", &*invs;
    for W in certs do
        print "ROOT_POWER_WITNESS", "p", W[1], "chi", W[2],
              "degrees_n2_to_12", W[3];
    end for;
end procedure;

function ScanTupleBank(label, filename, wanted)
    tuples := ReadTupleFile(filename);
    first := Max(1,start_bank_row);
    stop := Min(#tuples, first+max_bank_rows-1);
    print "BANK_START", label, filename, "rows_available", #tuples,
          "row_interval", [first..stop];
    exact_matches := 0;
    for i in [first..stop] do
        tup := tuples[i];
        f := SquareBranchModel(tup);
        if Discriminant(f) eq 0 then continue; end if;

        // Exact rational torsion is deliberately computed before Frobenius.
        invs, ord := ExactTorsion(f);
        print "BANK_ROW", label, i, tup, "exact", invs, "order", ord;
        if invs ne wanted then continue; end if;
        exact_matches +:= 1;

        certs := RootPowerCertificates(f,witnesses_per_curve);
        if #certs ge witnesses_per_curve then
            PrintSuccess(label,Sprintf("%o row %o tuple %o",filename,i,tup),
                         f,invs,certs);
            print "BANK_DONE", label, "tested", i,
                  "exact_matches", exact_matches, "success", true;
            return true;
        end if;
        if diagnose_failures then
            print "NO_ROOT_POWER_CERTIFICATE", label, "row", i,
                  "witnesses_found", #certs;
            try
                print "FAILURE_AUTOMORPHISM_GROUP_ORDER",
                      #AutomorphismGroup(HyperellipticCurve(f));
            catch e
                print "FAILURE_AUTOMORPHISM_GROUP_CHECK_FAILED", e`Object;
            end try;
            try
                print "FAILURE_RATIONAL_DEGREE2_SUBCOVERS",
                      #Degree2Subcovers(HyperellipticCurve(f));
            catch e
                print "FAILURE_DEGREE2_SUBCOVER_CHECK_FAILED", e`Object;
            end try;
        end if;
    end for;
    print "BANK_DONE", label, "tested_interval", [first..stop],
          "exact_matches", exact_matches, "success", false;
    return false;
end function;

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num,den) ne 1 then continue; end if;
            r := Q!num/den;
            key := Sprint(r);
            if key notin seen then
                Include(~seen,key);
                Append(~vals,r);
            end if;
        end for;
    end for;
    return vals;
end function;

function M244EllipticCurves(s,t)
    A := t^2 + 2*s;
    E2 := EllipticCurve([Q!0,Q!(2*A),Q!0,Q!(A^2-4*s^2),Q!0]);
    E3 := EllipticCurve([Q!0,Q!(-4*A),Q!0,Q!(16*s^2),Q!0]);
    return E2,E3;
end function;

function DualPhiE3ToE2(R,s,E2)
    E3 := Curve(R);
    if R eq E3!0 then return E2!0; end if;
    X := R[1]; Y := R[2];
    if X eq 0 then return E2!0; end if;
    return E2![Y^2/(4*X^2),Y*(16*s^2-X^2)/(8*X^2)];
end function;

function M244Polynomial(s,t,P1,R)
    E2 := Curve(P1);
    Q2 := P1 + DualPhiE3ToE2(R,s,E2);
    if P1 eq E2!0 or Q2 eq E2!0 then return false,P!0,_,_; end if;
    if P1[1] eq 0 or Q2[1] eq 0 then return false,P!0,_,_; end if;
    u := P1[2]/(2*P1[1]);
    v := Q2[2]/(2*Q2[1]);
    if u eq 0 or v eq 0 or u^2 eq v^2 then return false,P!0,_,_; end if;
    f := x*(x+u^2)*(x+v^2)*(x^2+(t^2+2*s)*x+s^2);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then return false,P!0,_,_; end if;
    return true,f,u,v;
end function;

function ScanM244Pilot()
    params := RationalParametersOfHeight(m244_height);
    seen := {};
    tested := 0;
    exact_matches := 0;
    print "M244_START", "height", m244_height,
          "point_bound", m244_point_bound,
          "curve_limit", max_m244_curves;
    for s in params do
        for t in params do
            if s eq 0 or t eq 0 or t^2+4*s eq 0 then continue; end if;
            E2,E3 := M244EllipticCurves(s,t);
            if Discriminant(E2) eq 0 or Discriminant(E3) eq 0 then continue; end if;
            pts2 := Points(E2 : Bound := m244_point_bound);
            pts3 := Points(E3 : Bound := m244_point_bound);
            if #pts2 le 4 or #pts3 le 2 then continue; end if;

            for P1 in pts2 do
                for R in pts3 do
                    ok,f,u,v := M244Polynomial(s,t,P1,R);
                    if not ok then continue; end if;
                    key := Sprint(f);
                    if key in seen then continue; end if;
                    Include(~seen,key);
                    tested +:= 1;

                    fI := IntegralSquareTwistModel(f);

                    // Exact rational torsion is deliberately computed first.
                    invs,ord := ExactTorsion(fI);
                    print "M244_CURVE", tested, "s", s, "t", t,
                          "u", u, "v", v, "exact", invs, "order", ord;
                    if invs eq [2,4,4] then
                        exact_matches +:= 1;
                        certs := RootPowerCertificates(fI,witnesses_per_curve);
                        if #certs ge witnesses_per_curve then
                            source := Sprintf(
                                "M(2,4,4) pilot s=%o t=%o u=%o v=%o",s,t,u,v);
                            PrintSuccess("[2,4,4]",source,fI,invs,certs);
                            try
                                print "M244_AUTOMORPHISM_GROUP_ORDER",
                                      #AutomorphismGroup(HyperellipticCurve(fI));
                            catch e
                                print "M244_AUTOMORPHISM_GROUP_CHECK_FAILED", e`Object;
                            end try;
                            try
                                print "M244_RATIONAL_DEGREE2_SUBCOVERS",
                                      #Degree2Subcovers(HyperellipticCurve(fI));
                            catch e
                                print "M244_DEGREE2_SUBCOVER_CHECK_FAILED", e`Object;
                            end try;
                            print "AUTOMORPHISM_NOTE root-power witness rules out geometric elliptic splitting for this curve";
                            print "M244_DONE", "tested", tested,
                                  "exact_matches", exact_matches, "success", true;
                            return true;
                        end if;
                    end if;
                    if tested ge max_m244_curves then
                        print "M244_DONE", "tested", tested,
                              "exact_matches", exact_matches, "success", false;
                        return false;
                    end if;
                end for;
            end for;
        end for;
    end for;
    print "M244_DONE", "tested", tested,
          "exact_matches", exact_matches, "success", false;
    return false;
end function;

success2244 := false;
success2228 := false;
success244 := false;

if target eq "all" or target eq "2244" then
    success2244 := ScanTupleBank(
        "[2,2,4,4]",
        "data/tor2244_bank.txt",
        [2,2,4,4]);
end if;

if target eq "all" or target eq "2228" then
    success2228 := ScanTupleBank(
        "[2,2,2,8]",
        "data/tor2228_bank.txt",
        [2,2,2,8]);
end if;

if target eq "all" or target eq "244" then
    success244 := ScanM244Pilot();
end if;

print "SUMMARY", "success_2244", success2244,
      "success_2228", success2228, "success_244", success244;
quit;
