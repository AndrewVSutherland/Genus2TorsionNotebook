//////////////////////////////////////////////////////////////////////
//  Projective q-chart check for M_1(8,2,2)+3 at p=7.
//
//  The finite monic q=X^2+UX+V boundary analysis finds no solutions
//  at p=7.  This script checks the projective q=[A:B:C] equation
//
//      h^2 - f = lambda*q^3
//
//  over F_7, including A=0 charts where the monic Mumford coefficients
//  have a pole modulo 7.
//////////////////////////////////////////////////////////////////////

F := GF(7);
P<X> := PolynomialRing(F);
max_print := 20;

function BoundaryLabels(u, v)
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
    if #labels eq 0 then Append(~labels, "boundary_unknown"); end if;
    return Join(Sort(labels), "+");
end function;

function QFlags(A, B, C, q, f)
    flags := [];
    if A eq 0 then Append(~flags, "A=0"); end if;
    if B^2 - 4*A*C eq 0 then Append(~flags, "discq=0"); end if;
    if Degree(GCD(q, f)) gt 0 then Append(~flags, "gcd(q,f)"); end if;
    if #flags eq 0 then Append(~flags, "q_open"); end if;
    return Join(Sort(flags), "+");
end function;

q_list := [];
for A in F do
    for B in F do
        for C in F do
            if A eq 0 and B eq 0 and C eq 0 then
                continue;
            end if;
            // Normalize projective q by first nonzero coordinate.
            if A ne 0 then
                scale := A^-1;
            elif B ne 0 then
                scale := B^-1;
            else
                scale := C^-1;
            end if;
            qkey := <A*scale, B*scale, C*scale>;
            if qkey notin q_list then
                Append(~q_list, qkey);
            end if;
        end for;
    end for;
end for;

base_open := 0;
base_boundary := 0;
solutions := 0;
open_qopen := 0;
base_with_solution := { Integers() | };
by_signature := AssociativeArray();
samples := [];

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
        if open then base_open +:= 1; else base_boundary +:= 1; end if;

        bkey := open select "open" else BoundaryLabels(u, v);
        base_int := Integers()!u + 7*(Integers()!v);

        for qcoords in q_list do
            A := qcoords[1]; B := qcoords[2]; Cc := qcoords[3];
            q := A*X^2 + B*X + Cc;
            q3 := q^3;
            qflag := QFlags(A, B, Cc, q, f);

            for h3 in F do
                for h2 in F do
                    for h1 in F do
                        for h0 in F do
                            h := h3*X^3 + h2*X^2 + h1*X + h0;
                            dpoly := h^2 - f;
                            idx := -1;
                            lam := F!0;
                            for i in [0..6] do
                                if Coefficient(q3, i) ne 0 then
                                    idx := i;
                                    lam := Coefficient(dpoly, i)/Coefficient(q3, i);
                                    break;
                                end if;
                            end for;
                            if idx lt 0 then
                                continue;
                            end if;
                            if dpoly ne lam*q3 then
                                continue;
                            end if;

                            solutions +:= 1;
                            Include(~base_with_solution, base_int);
                            if open and qflag eq "q_open" then
                                open_qopen +:= 1;
                            end if;
                            skey := bkey cat " | " cat qflag;
                            if IsDefined(by_signature, skey) then
                                by_signature[skey] +:= 1;
                            else
                                by_signature[skey] := 1;
                            end if;
                            if #samples lt max_print then
                                Append(~samples, <u,v,A,B,Cc,h3,h2,h1,h0,lam,bkey,qflag>);
                            end if;
                        end for;
                    end for;
                end for;
            end for;
        end for;
    end for;
end for;

print "M_1(8,2,2)+3 projective q analysis p=7";
print "projective_q_points", #q_list;
print "base_open", base_open, "base_boundary", base_boundary;
print "solutions", solutions, "bases_with_solution", #base_with_solution,
      "open_qopen", open_qopen;
print "samples", samples;
print "signature_counts";
for key in Sort([ k : k in Keys(by_signature) ]) do
    print " ", key, by_signature[key];
end for;

quit;
