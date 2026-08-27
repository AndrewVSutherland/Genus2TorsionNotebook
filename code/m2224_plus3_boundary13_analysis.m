//////////////////////////////////////////////////////////////////////
//  Component-wise p=13 boundary analysis for M(2,2,2,4)+3.
//
//  Base family:
//      C: y^2 = x (x+a^2) (x+b^2) (x+c^2) (x+d^2).
//
//  The good affine M(2,2,2,4) chart has no cubic-contact 3-torsion
//  points over F_13, so rational examples must reduce to the boundary:
//      Zi:   a_i = 0
//      Eij+: a_i =  a_j
//      Eij-: a_i = -a_j
//
//  This script solves the raw monic cubic-contact equations over F_13,
//  allowing bad cover reduction (v=0, repeated q, q meeting f), and
//  attributes every solution to the boundary components containing its
//  base point.
//
//  Typical run:
//      magma code/m2224_plus3_boundary13_analysis.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned max_print then
    max_print := 30;
elif Type(max_print) eq MonStgElt then
    max_print := StringToInteger(max_print);
end if;

p := 13;
F := GF(p);
P<x> := PolynomialRing(F);
Z := Integers();

names := ["a", "b", "c", "d"];

function JoinLabels(labels)
    if #labels eq 0 then
        return "";
    end if;
    out := labels[1];
    for i in [2..#labels] do
        out cat:= "+" cat labels[i];
    end for;
    return out;
end function;

function LabelKey(labels)
    return JoinLabels(Sort(labels));
end function;

function BaseIndex(vals)
    return Z!vals[1] + p*Z!vals[2] + p^2*Z!vals[3] + p^3*Z!vals[4];
end function;

function BoundaryLabels(vals)
    labels := [];
    for i in [1..4] do
        if vals[i] eq 0 then
            Append(~labels, "Z" cat IntegerToString(i));
        end if;
    end for;
    for i in [1..3] do
        for j in [i+1..4] do
            tag := IntegerToString(i) cat IntegerToString(j);
            if vals[i] - vals[j] eq 0 then
                Append(~labels, "E" cat tag cat "+");
            end if;
            if vals[i] + vals[j] eq 0 then
                Append(~labels, "E" cat tag cat "-");
            end if;
        end for;
    end for;
    if #labels eq 0 then
        Append(~labels, "open");
    end if;
    return labels;
end function;

function CurveKeyFromVals(vals)
    A := vals[1]^2;
    B := vals[2]^2;
    C := vals[3]^2;
    D := vals[4]^2;
    e1 := A+B+C+D;
    e2 := A*B + A*C + A*D + B*C + B*D + C*D;
    e3 := A*B*C + A*B*D + A*C*D + B*C*D;
    e4 := A*B*C*D;
    return <Z!e1, Z!e2, Z!e3, Z!e4>;
end function;

function PolynomialFromKey(key)
    return x^5 + (F!key[1])*x^4 + (F!key[2])*x^3 + (F!key[3])*x^2 + (F!key[4])*x;
end function;

function RawContactWitnessesByKey()
    by_key := AssociativeArray();
    raw_count := 0;

    for L in F do
        if L eq 0 then
            continue;
        end if;
        M := L^2;
        for U in F do
            for v in F do
                q := x^2 + U*x + v^2;
                for e1 in F do
                    A := 2*M*e1 + 6*(U^2+v^2) - (M+3*U)^2;
                    e2 := ((M+3*U)*A + 8*v^3 - 4*U^3 - 24*U*v^2)/(4*M);
                    e3 := (A^2 + 16*(M+3*U)*v^3
                           - 48*(U^2*v^2+v^4))/(16*M);
                    e4 := (A*v^3 - 6*U*v^4)/(2*M);
                    key := <Z!e1, Z!e2, Z!e3, Z!e4>;
                    witness := <Z!L, Z!U, Z!v>;
                    if IsDefined(by_key, key) then
                        seq := by_key[key];
                        Append(~seq, witness);
                        by_key[key] := seq;
                    else
                        by_key[key] := [witness];
                    end if;
                    raw_count +:= 1;
                end for;
            end for;
        end for;
    end for;

    return by_key, raw_count;
end function;

function CoverFlags(key, witness)
    f := PolynomialFromKey(key);
    U := F!witness[2];
    v := F!witness[3];
    q := x^2 + U*x + v^2;
    flags := [];
    if v eq 0 then
        Append(~flags, "v=0");
    end if;
    if U^2 - 4*v^2 eq 0 then
        Append(~flags, "discq=0");
    end if;
    if Degree(GCD(q, f)) gt 0 then
        Append(~flags, "gcd(q,f)");
    end if;
    if #flags eq 0 then
        Append(~flags, "cover_open");
    end if;
    return flags;
end function;

procedure Increment(~A, key, amount)
    if IsDefined(A, key) then
        A[key] +:= amount;
    else
        A[key] := amount;
    end if;
end procedure;

procedure IncludeBase(~A, key, base_key)
    if IsDefined(A, key) then
        S := A[key];
        Include(~S, base_key);
        A[key] := S;
    else
        A[key] := {Z | base_key};
    end if;
end procedure;

by_key, raw_count := RawContactWitnessesByKey();

print "M(2,2,2,4)+3 component-wise p=13 boundary analysis";
print "raw_contact_param", raw_count, "contact_coefficient_keys", #Keys(by_key);

base_component_count := AssociativeArray();
base_signature_count := AssociativeArray();
component_contact_count := AssociativeArray();
component_cover_open_count := AssociativeArray();
component_bases_with_contact := AssociativeArray();
component_bases_with_cover_open := AssociativeArray();
signature_contact_count := AssociativeArray();
cover_signature_count := AssociativeArray();
flag_count := AssociativeArray();

base_total := 0;
open_base := 0;
boundary_base := 0;
contact_solutions := 0;
open_contact_solutions := 0;
boundary_contact_solutions := 0;
cover_open_solutions := 0;
open_bases_with_contact := {Z |};
boundary_bases_with_contact := {Z |};
open_bases_with_cover_open := {Z |};
boundary_bases_with_cover_open := {Z |};
samples := [];
cover_open_samples := [];

for a in F do
    for b in F do
        for c in F do
            for d in F do
                vals := [a,b,c,d];
                labels := BoundaryLabels(vals);
                sig := LabelKey(labels);
                base_key := BaseIndex(vals);
                base_total +:= 1;
                Increment(~base_signature_count, sig, 1);
                for lab in labels do
                    Increment(~base_component_count, lab, 1);
                end for;
                if #labels eq 1 and labels[1] eq "open" then
                    open_base +:= 1;
                else
                    boundary_base +:= 1;
                end if;

                key := CurveKeyFromVals(vals);
                if not IsDefined(by_key, key) then
                    continue;
                end if;

                base_has_contact := false;
                base_has_cover_open := false;
                for witness in by_key[key] do
                    flags := CoverFlags(key, witness);
                    flag_key := LabelKey(flags);
                    is_cover_open := (#flags eq 1 and flags[1] eq "cover_open");
                    contact_solutions +:= 1;
                    base_has_contact := true;
                    Increment(~signature_contact_count, sig, 1);
                    Increment(~cover_signature_count, sig cat " | " cat flag_key, 1);
                    Increment(~flag_count, flag_key, 1);
                    for lab in labels do
                        Increment(~component_contact_count, lab, 1);
                    end for;

                    if #labels eq 1 and labels[1] eq "open" then
                        open_contact_solutions +:= 1;
                    else
                        boundary_contact_solutions +:= 1;
                    end if;

                    if is_cover_open then
                        cover_open_solutions +:= 1;
                        base_has_cover_open := true;
                        for lab in labels do
                            Increment(~component_cover_open_count, lab, 1);
                        end for;
                        if #cover_open_samples lt max_print then
                            Append(~cover_open_samples, <vals, sig, witness>);
                        end if;
                    end if;

                    if #samples lt max_print then
                        Append(~samples, <vals, sig, witness, flag_key>);
                    end if;
                end for;

                if base_has_contact then
                    if #labels eq 1 and labels[1] eq "open" then
                        Include(~open_bases_with_contact, base_key);
                    else
                        Include(~boundary_bases_with_contact, base_key);
                    end if;
                    for lab in labels do
                        IncludeBase(~component_bases_with_contact, lab, base_key);
                    end for;
                end if;
                if base_has_cover_open then
                    if #labels eq 1 and labels[1] eq "open" then
                        Include(~open_bases_with_cover_open, base_key);
                    else
                        Include(~boundary_bases_with_cover_open, base_key);
                    end if;
                    for lab in labels do
                        IncludeBase(~component_bases_with_cover_open, lab, base_key);
                    end for;
                end if;
            end for;
        end for;
    end for;
end for;

print "base_total", base_total, "open_base", open_base, "boundary_base", boundary_base;
print "contact_solutions", contact_solutions,
      "open_contact_solutions", open_contact_solutions,
      "boundary_contact_solutions", boundary_contact_solutions,
      "cover_open_solutions", cover_open_solutions;
print "open_bases_with_contact", #open_bases_with_contact,
      "boundary_bases_with_contact", #boundary_bases_with_contact,
      "open_bases_with_cover_open", #open_bases_with_cover_open,
      "boundary_bases_with_cover_open", #boundary_bases_with_cover_open;

print "";
print "component_summary";
for lab in Sort([k : k in Keys(base_component_count)]) do
    ctot := IsDefined(component_contact_count, lab) select component_contact_count[lab] else 0;
    cot := IsDefined(component_cover_open_count, lab) select component_cover_open_count[lab] else 0;
    bct := IsDefined(component_bases_with_contact, lab) select #component_bases_with_contact[lab] else 0;
    bot := IsDefined(component_bases_with_cover_open, lab) select #component_bases_with_cover_open[lab] else 0;
    print lab,
          "base_points", base_component_count[lab],
          "bases_with_contact", bct,
          "contact_solutions", ctot,
          "bases_with_cover_open", bot,
          "cover_open_solutions", cot;
end for;

print "";
print "cover_flag_counts";
for key in Sort([k : k in Keys(flag_count)]) do
    print key, flag_count[key];
end for;

print "";
print "top_boundary_signature_contact_counts";
for key in Sort([k : k in Keys(signature_contact_count)]) do
    if key ne "open" then
        print key, signature_contact_count[key];
    end if;
end for;

print "";
print "cover_signature_counts";
for key in Sort([k : k in Keys(cover_signature_count)]) do
    print key, cover_signature_count[key];
end for;

print "";
print "samples";
print samples;
print "cover_open_samples";
print cover_open_samples;

quit;
