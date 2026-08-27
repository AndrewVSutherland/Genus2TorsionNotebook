//////////////////////////////////////////////////////////////////////
//  Pull back the M(2,2,2,8) rational 3-torsion equations to section
//  families on the K3 surface.
//
//  The raw cubic-contact condition is
//
//      h(x)^2 - f(x) = m^2*q(x)^3
//
//  with q = x^2 + U*x + V and h = m*x^3 + N*x^2 + R*x + S.
//  Since S^2 = m^2*V^3 and gcd(q,f)=1, write V=v^2 and S=m*v^3.
//  Put L=1/m.  Eliminating N,R,S gives three equations in L,U,v:
//
//      A  = 2L^2*e1 + 6(U^2+v^2) - (L^2+3U)^2,
//      F1 = (L^2+3U)A + 8v^3 - 4e2L^2 - 4U^3 - 24Uv^2,
//      F2 = A^2 + 16(L^2+3U)v^3 - 16e3L^2 - 48(U^2v^2+v^4),
//      F3 = A*v^3 - 2e4L^2 - 6Uv^4,
//
//  where f = x^5 + e1*x^4 + e2*x^3 + e3*x^2 + e4*x.
//
//  A rational point with L*v*(U^2-4v^2) != 0 and gcd(q,f)=1 gives a
//  rational 3-torsion class.
//
//  Typical runs from torsion_jac:
//
//      magma -b mode:="summary" code/m2228_three_torsion_equations.m
//
//      magma -b mode:="local" code/m2228_three_torsion_equations.m
//
//      magma -b mode:="summary" equations_output_file:=data/m2228_3torsion_eqs.txt \
//          code/m2228_three_torsion_equations.m
//////////////////////////////////////////////////////////////////////

if not assigned mode then
    mode := "summary";
end if;

load "code/m2248_section_multiples.m";

Q := Rationals();

local_primes := [13,17,19,23,29,31];

function IntegralPrimitivePolynomial(f)
    coeffs := Coefficients(f);
    if #coeffs eq 0 then
        return f;
    end if;

    LCMden := LCM([ Denominator(c) : c in coeffs ]);
    g := Parent(f)!(LCMden*f);
    nums := [ Integers()!c : c in Coefficients(g) ];
    content := GCD([ Abs(n) : n in nums | n ne 0 ]);
    if content gt 1 then
        g := Parent(f)!(g/content);
    end if;
    return g;
end function;

function ClearDenominatorAsPolynomial(f, R)
    return IntegralPrimitivePolynomial(R!Numerator(f));
end function;

function PulledBackThreeTorsionEquations(family)
    R<t,L,U,v> := PolynomialRing(Q, 4);
    K := FieldOfFractions(R);
    PX<X> := PolynomialRing(K);

    vals := SectionTuple(family, K!t);
    f := X;
    for a in vals do
        f *:= X + a^2;
    end for;

    e1 := Coefficient(f, 4);
    e2 := Coefficient(f, 3);
    e3 := Coefficient(f, 2);
    e4 := Coefficient(f, 1);

    M := L^2;
    A := 2*M*e1 + 6*(U^2+v^2) - (M+3*U)^2;
    F1 := (M+3*U)*A + 8*v^3 - 4*e2*M - 4*U^3 - 24*U*v^2;
    F2 := A^2 + 16*(M+3*U)*v^3 - 16*e3*M
          - 48*(U^2*v^2+v^4);
    F3 := A*v^3 - 2*e4*M - 6*U*v^4;

    return [
        ClearDenominatorAsPolynomial(F1, R),
        ClearDenominatorAsPolynomial(F2, R),
        ClearDenominatorAsPolynomial(F3, R)
    ];
end function;

procedure PrintEquationSummary(family)
    Fs := PulledBackThreeTorsionEquations(family);
    print "FAMILY", family;
    for i in [1..#Fs] do
        print " F", i,
              "total_degree", TotalDegree(Fs[i]),
              "degree_t", Degree(Fs[i], 1),
              "degree_L", Degree(Fs[i], 2),
              "degree_U", Degree(Fs[i], 3),
              "degree_v", Degree(Fs[i], 4),
              "terms", #Terms(Fs[i]);
    end for;
end procedure;

procedure WriteEquations(output_file)
    out := Open(output_file, "w");
    for family in M2248SectionFamilyNames do
        Fs := PulledBackThreeTorsionEquations(family);
        fprintf out, "FAMILY %o\n", family;
        for i in [1..#Fs] do
            fprintf out, "F%o = %o\n", i, Fs[i];
        end for;
        fprintf out, "\n";
    end for;
    delete out;
    print "Wrote equations to", output_file;
end procedure;

function CurvePolynomialFromTupleFinite(vals, F)
    PF<X> := PolynomialRing(F);
    f := X;
    for a in vals do
        f *:= X + a^2;
    end for;
    return f;
end function;

function HasFiniteTripleContactSolution(f)
    F := BaseRing(Parent(f));
    PF<X> := PolynomialRing(F);

    e1 := Coefficient(f, 4);
    e2 := Coefficient(f, 3);
    e3 := Coefficient(f, 2);
    e4 := Coefficient(f, 1);

    for L in F do
        if L eq 0 then
            continue;
        end if;
        M := L^2;

        for U in F do
            for v in F do
                if v eq 0 or U^2 - 4*v^2 eq 0 then
                    continue;
                end if;

                q := X^2 + U*X + v^2;
                if Degree(GCD(q, f)) gt 0 then
                    continue;
                end if;

                A := 2*M*e1 + 6*(U^2+v^2) - (M+3*U)^2;
                F1 := (M+3*U)*A + 8*v^3 - 4*e2*M - 4*U^3 - 24*U*v^2;
                if F1 ne 0 then
                    continue;
                end if;

                F2 := A^2 + 16*(M+3*U)*v^3 - 16*e3*M
                      - 48*(U^2*v^2+v^4);
                if F2 ne 0 then
                    continue;
                end if;

                F3 := A*v^3 - 2*e4*M - 6*U*v^4;
                if F3 ne 0 then
                    continue;
                end if;

                return true, <L,U,v>;
            end for;
        end for;
    end for;

    return false, <F!0,F!0,F!0>;
end function;

procedure LocalTripleContactResidues()
    print "Finite-field triple-contact residues for named section families";
    print "A rational point must reduce either to a bad section/curve residue";
    print "or to one of the listed good triple-contact residues.";

    for family in M2248SectionFamilyNames do
        print "FAMILY", family;
        for p in local_primes do
            F := GF(p);
            good := 0;
            bad := 0;
            allowed_good := [];
            witnesses := [];

            for t0 in F do
                vals := [];
                ok := true;
                try
                    vals := SectionTuple(family, t0);
                catch e
                    ok := false;
                end try;

                if not ok or not IsUsableTuple(vals) then
                    bad +:= 1;
                    continue;
                end if;

                f := CurvePolynomialFromTupleFinite(vals, F);
                if Degree(f) ne 5 or Discriminant(f) eq 0 then
                    bad +:= 1;
                    continue;
                end if;

                good +:= 1;
                has_solution, witness := HasFiniteTripleContactSolution(f);
                if has_solution then
                    Append(~allowed_good, Integers()!t0);
                    if #witnesses lt 5 then
                        Append(~witnesses, <Integers()!t0, witness>);
                    end if;
                end if;
            end for;

            print " p", p, "good", good, "bad", bad,
                  "good_triple_contact", #allowed_good,
                  "residues", allowed_good,
                  "sample_witnesses", witnesses;
        end for;
    end for;
end procedure;

if mode eq "summary" then
    print "Pulled-back M(2,2,2,8) 3-torsion equation summaries";
    for family in M2248SectionFamilyNames do
        PrintEquationSummary(family);
    end for;
    if assigned equations_output_file then
        WriteEquations(equations_output_file);
    end if;
elif mode eq "local" then
    LocalTripleContactResidues();
else
    error "Unknown mode. Use mode:=\"summary\" or mode:=\"local\".";
end if;

quit;
