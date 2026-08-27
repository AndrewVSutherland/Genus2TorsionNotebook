//////////////////////////////////////////////////////////////////////
//  Intrinsic finite-field diagnostic for contact-7 + root + halving.
//
//  This avoids the s-parametrization and enumerates finite residues
//  (a,b,r) directly:
//
//      h = 1 - (7/2)x + a x^2 + b x^3,
//      f = (h^2 + (x - 1)^7)/x^2,
//      f(r) = 0.
//
//  It then imposes the unsolved first-halving equations for the root class
//  (r,0)-infinity and tests whether the resulting order-4 class is divisible
//  by 2 in J(F_p).  This includes finite boundary residues missed by the
//  s,u,z open chart.
//////////////////////////////////////////////////////////////////////

if not assigned p then
    p := 5;
elif Type(p) eq MonStgElt then
    p := StringToInteger(p);
end if;

if not assigned max_print then
    max_print := 40;
elif Type(max_print) eq MonStgElt then
    max_print := StringToInteger(max_print);
end if;

F := GF(p);
P<x> := PolynomialRing(F);

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

function IsSquareInF(t)
    if t eq 0 then
        return true;
    end if;
    return t^((#F - 1) div 2) eq 1;
end function;

function RootChartTag(r, hval)
    // In the rational-root chart, r=1-s^2 and h(r)=eps*s^7.
    if not IsSquareInF(1 - r) then
        return "not_s_chart";
    end if;
    if r eq 1 then
        return "s=0";
    end if;
    if r eq 0 then
        return "r=0";
    end if;
    return "open_s";
end function;

print "contact-7 intrinsic root + first-halving finite diagnostic";
print "p", p;

curve_checked := 0;
curve_good := 0;
root_good := 0;
surface := 0;
exact_h4 := 0;
h4_div2 := 0;
tag_counts := AssociativeArray();
tag_surface := AssociativeArray();
tag_exact := AssociativeArray();
tag_div2 := AssociativeArray();
samples := [];
div2_samples := [];

procedure Inc(~A, key)
    if not IsDefined(A, key) then
        A[key] := 0;
    end if;
    A[key] +:= 1;
end procedure;

for a in F do
    for b in F do
        curve_checked +:= 1;
        h := 1 - (F!7/F!2)*x + a*x^2 + b*x^3;
        f := ExactQuotient(h^2 + (x - 1)^7, x^2);
        if Degree(f) ne 5 or Discriminant(f) eq 0 then
            continue;
        end if;
        if Evaluate(h, F!1) eq 0 then
            continue;
        end if;
        curve_good +:= 1;

        C := HyperellipticCurve(f);
        J := Jacobian(C);
        D7 := J![x - 1, Evaluate(h, F!1)];

        for r in F do
            if Evaluate(f, r) ne 0 then
                continue;
            end if;
            root_good +:= 1;
            hval := Evaluate(h, r);
            tag := RootChartTag(r, hval);
            Inc(~tag_counts, tag);

            g := ExactQuotient(Evaluate(f, x + r), x);
            c3 := Coefficient(g, 3);
            c2 := Coefficient(g, 2);
            c1 := Coefficient(g, 1);
            c0 := Coefficient(g, 0);
            D2 := J![x - r, F!0];
            X := x - r;

            for u in F do
                for v in F do
                    if c3 ne u^2 - 2*v then
                        continue;
                    end if;
                    for w in F do
                        for z in F do
                            if c2 ne v^2 - 2*u*w + 2*z then
                                continue;
                            end if;
                            if c1 ne w^2 - 2*v*z then
                                continue;
                            end if;
                            if c0 ne z^2 then
                                continue;
                            end if;

                            surface +:= 1;
                            Inc(~tag_surface, tag);
                            if #samples lt max_print then
                                Append(~samples, <tag,a,b,r,u,v,w,z,#J,Order(D7)>);
                            end if;

                            Qpoly := X^2 - v*X + z;
                            alpha := (u*v - w)*X - u*z;
                            H4 := J![Qpoly, alpha];
                            if not (Order(D2) eq 2 and Order(H4) eq 4 and 2*H4 eq D2) then
                                continue;
                            end if;
                            exact_h4 +:= 1;
                            Inc(~tag_exact, tag);

                            if IsDivisibleBy2Finite(J, H4) then
                                h4_div2 +:= 1;
                                Inc(~tag_div2, tag);
                                if #div2_samples lt max_print then
                                    Append(~div2_samples, <tag,a,b,r,u,v,w,z,#J,Order(D7)>);
                                end if;
                            end if;
                        end for;
                    end for;
                end for;
            end for;
        end for;
    end for;
end for;

print "curve_checked", curve_checked;
print "curve_good", curve_good;
print "root_good", root_good;
print "surface", surface;
print "exact_h4", exact_h4;
print "h4_divisible_by_2", h4_div2;

print "TAG_ROOTS";
for key in Sort([ k : k in Keys(tag_counts) ]) do
    print " ", key, tag_counts[key];
end for;
print "TAG_SURFACE";
for key in Sort([ k : k in Keys(tag_surface) ]) do
    print " ", key, tag_surface[key];
end for;
print "TAG_EXACT_H4";
for key in Sort([ k : k in Keys(tag_exact) ]) do
    print " ", key, tag_exact[key];
end for;
print "TAG_DIV2";
for key in Sort([ k : k in Keys(tag_div2) ]) do
    print " ", key, tag_div2[key];
end for;

print "samples", samples;
print "div2_samples", div2_samples;

quit;
