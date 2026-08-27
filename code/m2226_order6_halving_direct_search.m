//////////////////////////////////////////////////////////////////////
//  Direct search for the suggested route: an order-6 class in the
//  M(2,2,2,6) family divisible by 2 in the Jacobian.
//
//  This uses Magma's exact IsDivisibleBy on specialized Jacobians rather
//  than symbolic descent squareclasses.  For each primitive [s:m:n], it
//  forms the odd quintic model from code/m2226_order6_doubling.m and
//  tests all 16 classes g + T with T in J[2].
//
//  Typical run from torsion_jac:
//      magma -b height:=20 code/m2226_order6_halving_direct_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 20;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;
if not assigned progress_interval then
    progress_interval := 500;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

Q := Rationals();
P<x> := PolynomialRing(Q);

function BranchData(vals)
    s := vals[1];
    m := vals[2];
    n := vals[3];
    A := [ Q!1, Q!1, Q!1, Q!2, Q!2 ];
    B := [
        2*s^2 - s*n,
        2*s^2 + s*m - 2*s*n - m*n,
        2*s^2 + s*m - s*n - m*n,
        -m*n,
        4*s^2 - 4*s*n - m*n
    ];
    return A, B;
end function;

function IsPrimitiveNormalizedTriple(a, b, c)
    if a eq 0 and b eq 0 and c eq 0 then
        return false;
    end if;
    if GCD([ Abs(a), Abs(b), Abs(c) ]) ne 1 then
        return false;
    end if;
    for z in [a,b,c] do
        if z ne 0 then
            return z gt 0;
        end if;
    end for;
    return false;
end function;

function HasDistinctEntries(vals)
    return #(Set(vals)) eq #vals;
end function;

function CurveData(vals)
    A, B := BranchData(vals);
    if &or [ b eq 0 : b in B ] then
        return false, _, _, _, _;
    end if;

    betas := [ -A[i]/B[i] : i in [1..5] ];
    if not HasDistinctEntries(betas) then
        return false, _, _, _, _;
    end if;

    f := &*[ A[i] + B[i]*x : i in [1..5] ];
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false, _, _, _, _;
    end if;

    return true, f, A, B, betas;
end function;

function TwoTorsionClasses(J, betas)
    out := [ <"0", J!0> ];
    for i in [1..5] do
        Append(~out, <Sprintf("W%o-inf", i), J![x-betas[i], Q!0]>);
    end for;
    for i in [1..4] do
        for j in [i+1..5] do
            Append(~out, <Sprintf("W%o+W%o", i, j), J![(x-betas[i])*(x-betas[j]), Q!0]>);
        end for;
    end for;
    return out;
end function;

print "Direct order-6 halving search";
print "height", height;
print "max_hits", max_hits;

hits := [];
checked_curves := 0;
checked_classes := 0;

for ss in [-height..height] do
    for mm in [-height..height] do
        for nn in [-height..height] do
            if #hits ge max_hits then
                break ss;
            end if;
            if not IsPrimitiveNormalizedTriple(ss, mm, nn) then
                continue;
            end if;
            vals := [ Q!ss, Q!mm, Q!nn ];
            ok, f, A, B, betas := CurveData(vals);
            if not ok then
                continue;
            end if;

            C := HyperellipticCurve(f);
            J := Jacobian(C);
            g := J![x, Q!2];
            checked_curves +:= 1;

            for data in TwoTorsionClasses(J, betas) do
                label := data[1];
                T := data[2];
                D := g + T;
                checked_classes +:= 1;

                ord := Order(D);
                if ord ne 6 then
                    continue;
                end if;

                if IsDivisibleBy(D, 2) then
                    Append(~hits, <vals, label, f>);
                    print "HIT vals", vals, "class", label, "f", f;
                end if;
            end for;

            if progress_interval gt 0 and checked_curves mod progress_interval eq 0 then
                print "checked curves", checked_curves, "classes", checked_classes, "hits", #hits;
            end if;
        end for;
    end for;
end for;

print "Done";
print "checked curves", checked_curves;
print "checked classes", checked_classes;
print "hits", #hits;
for H in hits do
    print H;
end for;
quit;
