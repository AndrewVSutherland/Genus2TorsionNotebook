//////////////////////////////////////////////////////////////////////
//  claude_ov_2216_exact3.m
//
//  Exact test of candidate (u,v) parameters on the odd M_1(8,2,2)
//  model.  For each parameter we report:
//    * the exact TorsionSubgroup invariants  (this decides ALL 16
//      twists at once: a rational half of ANY order-8 class is an
//      order-16 rational torsion point, so [2,2,16] <= torsion iff
//      some order-8 class halves)
//    * the order of the marked class Q = (-1, u*v*(u+v+1))
//    * whether Q itself is divisible by 2
//    * the number of order-8 classes in J(Q)_tors that are divisible
//      by 2 (the "twist" count)
//    * the 2-rank
//
//  Usage:
//    magma -b candidate_file:=data/xxx.txt code/claude_ov_2216_exact3.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned candidate_file then
    candidate_file := "data/m3222_halving_candidates_h80_p73_open_py.txt";
end if;
if not assigned max_exact then
    max_exact := 1000000;
elif Type(max_exact) eq MonStgElt then
    max_exact := StringToInteger(max_exact);
end if;
if not assigned MemGB then
    MemGB := 8;
elif Type(MemGB) eq MonStgElt then
    MemGB := StringToInteger(MemGB);
end if;
SetMemoryLimit(MemGB*10^9);

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

// odd (quintic) model of M_1(8,2,2); marked point Q = (-1, u*v*(u+v+1))
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

// 2-rank from the factorization type of the (integral) quintic
function TwoRankQuintic(fI)
    fac := Factorization(fI);
    nlin := #[t : t in fac | Degree(t[1]) eq 1 and t[2] eq 1];
    // quintic: Weierstrass pts = roots + one at infinity.
    // rational 2-torsion rank = (#rational Weierstrass points) - 1
    return nlin;   // (#roots) + inf - 1 = #roots
end function;

lines := Split(Read(candidate_file), "\n");
printf "claude_ov_2216_exact3  candidate_file=%o  max_exact=%o\n",
       candidate_file, max_exact;

rows := 0; curve_ok := 0; exact := 0; hits := 0;

for raw in lines do
    if exact ge max_exact then break; end if;
    s := raw;
    if #s gt 0 and s[#s] eq "\r" then s := s[1..#s-1]; end if;
    if #s eq 0 or s[1] eq "#" then continue; end if;
    parts := [part : part in Split(s, " ") | #part gt 0];
    if #parts lt 2 then continue; end if;
    rows +:= 1;
    u := ParseQ(parts[1]);
    v := ParseQ(parts[2]);

    if HasBasicDegeneracy(u, v) then
        printf "SKIP degenerate u=%o v=%o\n", u, v;
        continue;
    end if;
    f := OddPolynomial(u, v);
    yq := u*v*(u+v+1);
    if Degree(f) ne 5 or Discriminant(f) eq 0 then
        printf "SKIP singular u=%o v=%o\n", u, v;
        continue;
    end if;
    if yq eq 0 or yq^2 ne Evaluate(f, Qq!-1) then
        printf "SKIP markedpoint u=%o v=%o\n", u, v;
        continue;
    end if;
    curve_ok +:= 1;

    fI, L := IntegralModelPolynomial(f);
    if Degree(fI) ne 5 or Discriminant(fI) eq 0 then
        printf "SKIP intmodel u=%o v=%o\n", u, v;
        continue;
    end if;

    C := HyperellipticCurve(fI);
    J := Jacobian(C);
    DQ := J![X + 1, Qq!(L*yq)];
    exact +:= 1;

    ordQ := Order(DQ);
    isdiv := IsDivisibleBy(DQ, 2);
    tr := TwoRankQuintic(fI);

    G, phi := TorsionSubgroup(J);
    invs := Invariants(G);

    // count order-8 classes that are 2-divisible
    n8 := 0; n8div := 0;
    for g in G do
        Dg := phi(g);
        if Order(Dg) eq 8 then
            n8 +:= 1;
            if IsDivisibleBy(Dg, 2) then n8div +:= 1; end if;
        end if;
    end for;

    printf "RESULT u=%o v=%o torsion=%o ordQ=%o Qdiv2=%o tworank=%o n_order8=%o n_order8_div2=%o\n",
           u, v, invs, ordQ, isdiv, tr, n8, n8div;
    printf "   fI = %o\n", fI;
    printf "   L  = %o\n", L;
    if n8div gt 0 then
        hits +:= 1;
        printf "HIT u=%o v=%o torsion=%o\n", u, v, invs;
    end if;
end for;

printf "rows=%o curve_ok=%o exact=%o hits=%o\n", rows, curve_ok, exact, hits;
print "SEARCH_DONE";
quit;
