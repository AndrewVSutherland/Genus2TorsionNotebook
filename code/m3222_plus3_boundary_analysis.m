//////////////////////////////////////////////////////////////////////
//  Component-wise boundary analysis for M_1(8,2,2) / [2,2,8]
//  plus possible rational 3-torsion.
//
//  The open necessary condition for rational 3-torsion is empty at
//  p=7,11,13, so any rational example must reduce to boundary at all
//  three primes.  This script solves the raw cubic-contact equations
//  modulo those primes on all finite (u,v) residues, including bad
//  boundary fibers, and classifies the live boundary components.
//
//  We avoid the monic eliminated formulas because the M_1(8,2,2)
//  quintic is not generally monic.  Instead, for
//
//      h = m X^3 + N X^2 + R X + S,
//      q = X^2 + U X + V,
//
//  we impose h^2 - f = m^2 q^3 directly.  As usual, the constant
//  coefficient equation gives V=w^2, S=m w^3; varying w gives both
//  signs of S.
//
//  Typical run:
//      magma code/m3222_plus3_boundary_analysis.m
//////////////////////////////////////////////////////////////////////

if not assigned max_print then
    max_print := 12;
elif Type(max_print) eq MonStgElt then
    max_print := StringToInteger(max_print);
end if;

prime_list := [7,11,13];

function BoundaryLabels(u, v, F)
    s := u + v;
    p := u*v;
    labels := [];
    if p eq 0 then Append(~labels, "p=0"); end if;
    if s + 1 eq 0 then Append(~labels, "s+1=0"); end if;
    if s + 2 eq 0 then Append(~labels, "s+2=0"); end if;
    if p - s + 1 eq 0 then Append(~labels, "p-s+1=0"); end if;
    if s - p + 1 eq 0 then Append(~labels, "s-p+1=0"); end if;
    if s^2 - 4*p eq 0 then Append(~labels, "s^2-4p=0"); end if;
    if 2*s^2 + 3*s + p + 1 eq 0 then Append(~labels, "delta1=0"); end if;
    if s^3 - s^2*p + s^2 - 4*s*p - 4*p eq 0 then Append(~labels, "delta2=0"); end if;
    if #labels eq 0 then
        Append(~labels, "boundary_unknown");
    end if;
    return labels;
end function;

function LabelKey(labels)
    return Join(Sort(labels), "+");
end function;

function CoverFlags(U, w, q, f, F)
    flags := [];
    V := w^2;
    if w eq 0 then Append(~flags, "w=0"); end if;
    if U^2 - 4*V eq 0 then Append(~flags, "discq=0"); end if;
    if Degree(GCD(q, f)) gt 0 then Append(~flags, "gcd(q,f)"); end if;
    if #flags eq 0 then Append(~flags, "cover_open"); end if;
    return flags;
end function;

for p in prime_list do
    F := GF(p);
    P<X> := PolynomialRing(F);

    base_open := 0;
    base_boundary := 0;
    base_open_contact := { Integers() | };
    base_boundary_contact := { Integers() | };
    contact_solutions := 0;
    open_nondeg_contact := 0;
    boundary_nondeg_contact := 0;
    by_base_signature := AssociativeArray();
    by_cover_signature := AssociativeArray();
    samples := [];
    open_samples := [];

    for u in F do
        for v in F do
            qtilde := -X^2
                + (u^2*v - u^2 + u*v^2 - u - v^2 - v - 2)*X
                - (u^2 + u*v + v^2 + u + v + 1);
            f := ((1-u)*X + 1)*((1-v)*X + 1)*((u+v+2)*X + 1)*qtilde;
            yq := u*v*(u+v+1);
            open := true;
            if u eq 0 or v eq 0 or u eq v then open := false; end if;
            if u eq 1 or v eq 1 or u+v+1 eq 0 or u+v+2 eq 0 then open := false; end if;
            if Degree(f) ne 5 or Discriminant(f) eq 0 or yq eq 0 then open := false; end if;
            if yq^2 ne Evaluate(f, F!-1) then open := false; end if;
            if open then
                base_open +:= 1;
            else
                base_boundary +:= 1;
            end if;

            a5 := Coefficient(f, 5);
            a4 := Coefficient(f, 4);
            a3 := Coefficient(f, 3);
            a2 := Coefficient(f, 2);
            a1 := Coefficient(f, 1);
            a0 := Coefficient(f, 0);

            base_key_int := Integers()!u + p*(Integers()!v);
            found_on_base := false;
            found_nondeg_on_base := false;

            for m in F do
                if m eq 0 then
                    continue;
                end if;
                for U in F do
                    for w in F do
                        V := w^2;
                        S := m*w^3;
                        N := (a5 + 3*m^2*U)/(2*m);
                        R := (a4 + 3*m^2*(U^2 + V) - N^2)/(2*m);

                        h := m*X^3 + N*X^2 + R*X + S;
                        q := X^2 + U*X + V;
                        lhs := h^2 - f - m^2*q^3;
                        if lhs ne 0 then
                            continue;
                        end if;

                        contact_solutions +:= 1;
                        found_on_base := true;
                        flags := CoverFlags(U, w, q, f, F);
                        cover_open := (#flags eq 1 and flags[1] eq "cover_open");
                        if open and cover_open then
                            open_nondeg_contact +:= 1;
                            found_nondeg_on_base := true;
                            if #open_samples lt max_print then
                                Append(~open_samples, <u,v,m,U,w>);
                            end if;
                        elif (not open) and cover_open then
                            boundary_nondeg_contact +:= 1;
                            found_nondeg_on_base := true;
                        end if;

                        labels := open select ["open"] else BoundaryLabels(u, v, F);
                        bkey := LabelKey(labels);
                        ckey := bkey cat " | " cat LabelKey(flags);
                        if IsDefined(by_base_signature, bkey) then
                            by_base_signature[bkey] +:= 1;
                        else
                            by_base_signature[bkey] := 1;
                        end if;
                        if IsDefined(by_cover_signature, ckey) then
                            by_cover_signature[ckey] +:= 1;
                        else
                            by_cover_signature[ckey] := 1;
                        end if;

                        if #samples lt max_print then
                            Append(~samples, <u,v,m,U,w,bkey,LabelKey(flags)>);
                        end if;
                    end for;
                end for;
            end for;

            if found_on_base then
                if open then
                    Include(~base_open_contact, base_key_int);
                else
                    Include(~base_boundary_contact, base_key_int);
                end if;
            end if;
        end for;
    end for;

    print "M_1(8,2,2)+3 boundary contact analysis p", p;
    print "base_open", base_open, "base_boundary", base_boundary;
    print "contact_solutions", contact_solutions,
          "open_bases_with_contact", #base_open_contact,
          "boundary_bases_with_contact", #base_boundary_contact;
    print "open_nondeg_contact", open_nondeg_contact,
          "boundary_nondeg_contact", boundary_nondeg_contact;
    if #open_samples gt 0 then
        print "OPEN_NONDEG_SAMPLES", open_samples;
    end if;
    print "samples", samples;
    print "base_signature_counts";
    for key in Sort([ k : k in Keys(by_base_signature) ]) do
        print " ", key, by_base_signature[key];
    end for;
    print "cover_signature_counts";
    for key in Sort([ k : k in Keys(by_cover_signature) ]) do
        print " ", key, by_cover_signature[key];
    end for;
end for;

quit;
