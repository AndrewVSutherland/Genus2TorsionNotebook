//////////////////////////////////////////////////////////////////////
//  Search the first route for M(2,2,2,6): divide a rational 2-torsion
//  class by 2.
//
//  We use the odd quintic model from code/m2226_order6_doubling.m.
//  For each ordered pair of branch points (infinity, zero), normalize
//  the four remaining branch points by cross-ratios.  The 2-torsion
//  class zero-infinity is divisible by 2 iff the four normalized values
//  are all squares, equivalently their ratios to one reference value are
//  squares.
//
//  Typical run from torsion_jac:
//      magma -b height:=50 code/m2226_order2_halving_search.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 30;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_hits then
    max_hits := 100;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;

Q := Rationals();
R<s,m,n> := PolynomialRing(Q, 3);
PX<X> := PolynomialRing(R);

A := [ R!1, R!1, R!1, R!2, R!2 ];
B := [
    2*s^2 - s*n,
    2*s^2 + s*m - 2*s*n - m*n,
    2*s^2 + s*m - s*n - m*n,
    -m*n,
    4*s^2 - 4*s*n - m*n
];
L := [ PX!A[i] + PX!B[i]*X : i in [1..5] ];
f := &*L;
f_disc := Discriminant(f);

// Branch points on P^1.  Index 1 is the old infinity of the odd model;
// indices 2..6 are roots of L_1..L_5, represented as (-A_i : B_i).
branches := [ <R!1, R!0> ] cat [ < -A[i], B[i] > : i in [1..5] ];

function DetP(P, Qp)
    return P[1]*Qp[2] - Qp[1]*P[2];
end function;

function EvalR(g, vals)
    return Q!Evaluate(g, vals);
end function;

function EvalRatio(num, den, vals)
    denq := EvalR(den, vals);
    if denq eq 0 then
        return false, Q!0;
    end if;
    return true, EvalR(num, vals)/denq;
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

function HalvingWitnessForOrderedPair(vals, i_inf, i_zero)
    remaining := [ k : k in [1..6] | k ne i_inf and k ne i_zero ];
    if #remaining ne 4 then
        error "Bad remaining branch count.";
    end if;

    // lambda_k = det(branch_k, zero)/det(branch_k, infinity).
    // We only need lambda_k/lambda_ref.
    ref := remaining[1];
    for idx in [2..#remaining] do
        k := remaining[idx];
        num := DetP(branches[k], branches[i_zero]) * DetP(branches[ref], branches[i_inf]);
        den := DetP(branches[k], branches[i_inf]) * DetP(branches[ref], branches[i_zero]);
        ok, ratio := EvalRatio(num, den, vals);
        if not ok or ratio eq 0 then
            return false, [];
        end if;
        ok_square, root := IsSquare(ratio);
        if not ok_square then
            return false, [];
        end if;
    end for;

    return true, remaining;
end function;

print "Searching primitive [s:m:n] with height <=", height;
print "Max hits:", max_hits;

hits := [];
seen := {};
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
            if EvalR(f_disc, vals) eq 0 then
                continue;
            end if;
            if &or [ EvalR(Bi, vals) eq 0 : Bi in B ] then
                continue;
            end if;

            for i_inf in [1..6] do
                for i_zero in [1..6] do
                    if i_inf eq i_zero then
                        continue;
                    end if;
                    // Avoid reporting the same 2-torsion class twice.
                    unordered := Sort([i_inf, i_zero]);
                    key := Sprint(vals) cat ":" cat Sprint(unordered);
                    if key in seen then
                        continue;
                    end if;
                    ok, remaining := HalvingWitnessForOrderedPair(vals, i_inf, i_zero);
                    if ok then
                        Include(~seen, key);
                        Append(~hits, <vals, unordered, [i_inf, i_zero], remaining>);
                        print "hit vals", vals, "class", unordered, "ordered", [i_inf, i_zero], "remaining", remaining;
                    end if;
                end for;
            end for;
        end for;
    end for;
end for;

print "Total hits", #hits;
for H in hits do
    print H;
end for;
quit;
