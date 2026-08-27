//////////////////////////////////////////////////////////////////////
//  Search rational points on the M(2,2,2,6) theta-doubling curves.
//
//  This is a companion to code/m2226_order6_doubling.m.  It searches
//  primitive triples [s:m:n] and keeps only points where the odd quintic
//  is nonsingular and the corresponding coset element has order 6.
//
//  Typical run from torsion_jac:
//      magma -b height:=50 code/m2226_search_doubling_curves.m
//////////////////////////////////////////////////////////////////////

if not assigned height then
    height := 30;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned max_hits_per_pair then
    max_hits_per_pair := 20;
elif Type(max_hits_per_pair) eq MonStgElt then
    max_hits_per_pair := StringToInteger(max_hits_per_pair);
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

function ProductExcept(seq, omitted)
    out := Universe(seq)!1;
    for i in [1..#seq] do
        if i notin omitted then
            out *:= seq[i];
        end if;
    end for;
    return out;
end function;

function PairRawU(i, j)
    G := ProductExcept(L, {i,j});
    numerator := G - (Q!4)/(A[i]^2*A[j]^2)*L[i]*L[j];
    ok, Bij := IsDivisibleBy(numerator, X);
    assert ok;
    return Bij;
end function;

function PairDiscriminant(i, j)
    return Discriminant(PairRawU(i, j));
end function;

function EvalR(g, vals)
    return Q!Evaluate(g, vals);
end function;

function EvalPX(poly, vals, xvar)
    Qx := Parent(xvar);
    out := Qx!0;
    for k in [0..Degree(poly)] do
        out +:= EvalR(Coefficient(poly, k), vals)*xvar^k;
    end for;
    return out;
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

function CandidateIsValid(vals, i, j)
    Qx<x> := PolynomialRing(Q);
    Lq := [ EvalPX(Li, vals, x) : Li in L ];
    fq := &*Lq;
    if Discriminant(fq) eq 0 then
        return false, "singular", 0;
    end if;

    Uraw := EvalPX(PairRawU(i, j), vals, x);
    if Degree(Uraw) ne 2 or LeadingCoefficient(Uraw) eq 0 then
        return false, "bad_U", 0;
    end if;
    if Discriminant(Uraw) ne 0 then
        return false, "disc_nonzero", 0;
    end if;

    C := HyperellipticCurve(fq);
    J := Jacobian(C);
    g := J![x, Q!2];

    Aq := [ EvalR(Ai, vals) : Ai in A ];
    Bq := [ EvalR(Bi, vals) : Bi in B ];
    if Bq[i] eq 0 or Bq[j] eq 0 then
        return false, "root_at_infinity", 0;
    end if;
    beta_i := -Aq[i]/Bq[i];
    beta_j := -Aq[j]/Bq[j];
    Tij := J![(x-beta_i)*(x-beta_j), Q!0];
    elt := g + Tij;
    ord := Order(elt);
    if ord ne 6 then
        return false, "order_" cat IntegerToString(ord), ord;
    end if;

    return true, "ok", ord;
end function;

pairs := [];
for i in [1..4] do
    for j in [i+1..5] do
        if not (i eq 3 and j eq 5) then
            Append(~pairs, <i,j>);
        end if;
    end for;
end for;
hits := [* *];
for pair in pairs do
    Append(~hits, []);
end for;

print "Searching primitive [s:m:n] with height <=", height;
print "Max hits per pair:", max_hits_per_pair;

for ss in [-height..height] do
    for mm in [-height..height] do
        for nn in [-height..height] do
            if not IsPrimitiveNormalizedTriple(ss, mm, nn) then
                continue;
            end if;
            vals := [ Q!ss, Q!mm, Q!nn ];
            if EvalR(f_disc, vals) eq 0 then
                continue;
            end if;

            for pair_idx in [1..#pairs] do
                pair := pairs[pair_idx];
                if #hits[pair_idx] ge max_hits_per_pair then
                    continue;
                end if;
                i := pair[1];
                j := pair[2];
                disc := PairDiscriminant(i, j);
                if EvalR(disc, vals) ne 0 then
                    continue;
                end if;
                ok, reason, ord := CandidateIsValid(vals, i, j);
                if ok then
                    Append(~hits[pair_idx], vals);
                    print "hit pair", [i,j], "vals", vals;
                end if;
            end for;
        end for;
    end for;
end for;

print "Summary";
for pair_idx in [1..#pairs] do
    pair := pairs[pair_idx];
    print [pair[1], pair[2]], #hits[pair_idx], hits[pair_idx];
end for;
