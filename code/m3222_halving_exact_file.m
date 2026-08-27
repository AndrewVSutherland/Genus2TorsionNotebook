//////////////////////////////////////////////////////////////////////
//  Exact-test candidate (u,v) pairs for halving the distinguished
//  order-8 class in M_1(8,2,2).
//////////////////////////////////////////////////////////////////////

if not assigned candidate_file then
    candidate_file := "data/m3222_halving_candidates_h50_p43.txt";
end if;
if not assigned max_exact then
    max_exact := 1000000;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;

Qq := Rationals();
P<X> := PolynomialRing(Qq);

function ParseQ(s)
    if "/" in s then
        parts := Split(s, "/");
        return Qq!StringToInteger(parts[1]) / Qq!StringToInteger(parts[2]);
    end if;
    return Qq!StringToInteger(s);
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L, Denominator(Coefficient(f, i)));
    end for;
    return P!(L^2*f), L;
end function;

function OddPolynomial(u, v)
    qtilde := -X^2
        + (u^2*v - u^2 + u*v^2 - u - v^2 - v - 2)*X
        - (u^2 + u*v + v^2 + u + v + 1);
    return ((1-u)*X + 1)*((1-v)*X + 1)*((u+v+2)*X + 1)*qtilde;
end function;

function HasBasicDegeneracy(u, v)
    if u eq 0 or v eq 0 or u eq v then return true; end if;
    if u eq 1 or v eq 1 or u+v+1 eq 0 or u+v+2 eq 0 then return true; end if;
    return false;
end function;

function CurveData(u, v)
    if HasBasicDegeneracy(u, v) then
        return false, _, _;
    end if;
    f := OddPolynomial(u, v);
    yq := u*v*(u+v+1);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        return false, _, _;
    end if;
    if yq eq 0 or yq^2 ne Evaluate(f, Qq!-1) then
        return false, _, _;
    end if;
    return true, f, yq;
end function;

function TorsionInvariants(J)
    G, phi := TorsionSubgroup(J);
    return Invariants(G);
end function;

lines := Split(Read(candidate_file), "\n");
print "Exact M_1(8,2,2) halving candidate-file test";
print "candidate_file", candidate_file, "max_exact", max_exact;

rows := 0;
curve_ok := 0;
exact := 0;
order8 := 0;
divisible := 0;
hits := [];

for raw in lines do
    if #hits ge max_hits or exact ge max_exact then
        break;
    end if;
    s := raw;
    if #s gt 0 and s[#s] eq "\r" then
        s := s[1..#s-1];
    end if;
    if #s eq 0 or s[1] eq "#" then
        continue;
    end if;
    parts := [part : part in Split(s, " ") | #part gt 0];
    if #parts lt 2 then
        continue;
    end if;
    rows +:= 1;
    u := ParseQ(parts[1]);
    v := ParseQ(parts[2]);

    ok, f, yq := CurveData(u, v);
    if not ok then
        continue;
    end if;
    curve_ok +:= 1;
    fI, L := IntegralModelPolynomial(f);
    if Degree(fI) ne 5 or Discriminant(fI) eq 0 then
        continue;
    end if;

    C := HyperellipticCurve(fI);
    J := Jacobian(C);
    DQ := J![X + 1, Qq!(L*yq)];
    exact +:= 1;
    ordQ := Order(DQ);
    if ordQ eq 8 then
        order8 +:= 1;
    end if;

    isdiv := IsDivisibleBy(DQ, 2);
    if isdiv then
        _, half := IsDivisibleBy(DQ, 2);
        divisible +:= 1;
        invs := TorsionInvariants(J);
        Append(~hits, <u, v, ordQ, Order(half), invs, fI, L>);
        print "HIT", "u", u, "v", v,
              "order_Q", ordQ, "half_order", Order(half), "torsion", invs;
        print "  fI", fI;
    end if;
end for;

print "DONE";
print "rows", rows;
print "curve_ok", curve_ok;
print "exact", exact;
print "order8", order8;
print "divisible", divisible;
print "hits", #hits;
for H in hits do
    print H;
end for;

quit;
