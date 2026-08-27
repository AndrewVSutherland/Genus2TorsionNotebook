//////////////////////////////////////////////////////////////////////
//  claude_ov_2216_delta_validate.m
//
//  Validate the x-T module code/claude_ov_2216_delta.m against Magma's
//  IsDivisibleBy(D,2) on the odd M_1(8,2,2) family.
//
//  Theorem under test (x-T injective on J(Q)/2J(Q), odd degree):
//        delta(D) trivial   <==>   D in 2 J(Q).
//
//  Every nonzero element of the full torsion subgroup is tested, so the
//  test has genuine POSITIVES (2Q, 4Q, 6Q are always 2-divisible).
//
//  Usage: magma -b nsample:=NN height:=HH code/claude_ov_2216_delta_validate.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
load "code/claude_ov_2216_delta.m";

if not assigned height then height := 12;
elif Type(height) eq MonStgElt then height := StringToInteger(height); end if;
if not assigned nsample then nsample := 8;
elif Type(nsample) eq MonStgElt then nsample := StringToInteger(nsample); end if;
if not assigned MemGB then MemGB := 4;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

P<X> := PP;   // polynomial ring from the module

function OddPolynomial(u, v)
    qt := -X^2
        + (u^2*v - u^2 + u*v^2 - u - v^2 - v - 2)*X
        - (u^2 + u*v + v^2 + u + v + 1);
    return ((1-u)*X + 1)*((1-v)*X + 1)*((u+v+2)*X + 1)*qt, qt;
end function;

hs := [];
for a in [-height..height] do
    for b in [1..height] do
        if GCD(AbsoluteValue(a), b) eq 1 then Append(~hs, Qq!a/b); end if;
    end for;
end for;
hs := Sort(Setseq(Seqset(hs)));

printf "claude_ov_2216_delta_validate height=%o nsample=%o #rationals=%o\n",
        height, nsample, #hs;

tested := 0; totalD := 0; agree := 0; disagree := 0; pos := 0;
homoerr := 0; skipped := 0;
SetSeed(20260725);

for iter in [1..500000] do
    if tested ge nsample then break; end if;
    u := Random(hs); v := Random(hs);
    if u eq 0 or v eq 0 or u eq v or u eq 1 or v eq 1 then continue; end if;
    if u+v+1 eq 0 or u+v+2 eq 0 then continue; end if;
    f0, qt := OddPolynomial(u, v);
    if Degree(f0) ne 5 or Discriminant(f0) eq 0 then continue; end if;
    if not IsIrreducible(qt) then continue; end if;
    y0 := u*v*(u+v+1);
    if y0 eq 0 or y0^2 ne Evaluate(f0, Qq!-1) then continue; end if;

    Lden := 1;
    for i in [0..Degree(f0)] do Lden := LCM(Lden, Denominator(Coefficient(f0,i))); end for;
    f := P!(Lden^2*f0);
    C := HyperellipticCurve(f);
    J := Jacobian(C);
    ok := true;
    try DQ := J![X+1, Qq!(Lden*y0)]; catch e ok := false; end try;
    if not ok then continue; end if;
    if Order(DQ) ne 8 then continue; end if;

    roots := [ Qq!rt[1] : rt in Roots(f) ];
    if #roots ne 3 then continue; end if;
    qtm := qt / LeadingCoefficient(qt);
    KK<th> := NumberField(qtm);
    tested +:= 1;

    // the eight rational 2-torsion classes with their delta'
    lin := [ X - r : r in roots ];
    T2 := [* <J!0, <Qq!1,Qq!1,Qq!1,KK!1>> *];
    subsets := [ qtm ];
    for i in [1..3] do for j in [i+1..3] do
        Append(~subsets, lin[i]*lin[j]);
    end for; end for;
    for i in [1..3] do
        Append(~subsets, qtm * &*[ lin[j] : j in [1..3] | j ne i ]);
    end for;
    for uS in subsets do
        gotit, dp := CO_DeltaPrimeTwo(uS, f, roots, KK, th);
        if not gotit then continue; end if;
        DS := J![uS, P!0];
        Append(~T2, <DS, dp>);
    end for;
    if #T2 ne 8 then printf "WARN u=%o v=%o #T2=%o\n", u, v, #T2; end if;

    // internal homomorphism check on J[2] (mod squares and mod Q*)
    for i in [1..#T2] do for j in [i..#T2] do
        S := T2[i][1] + T2[j][1];
        dpred := CO_MulPrime(T2[i][2], T2[j][2]);
        found := false;
        for k in [1..#T2] do
            if T2[k][1] eq S then
                ratio := CO_MulPrime(dpred, T2[k][2]);
                s1 := CO_SqClass(ratio[1]);
                okh := (CO_SqClass(ratio[2]) eq s1) and (CO_SqClass(ratio[3]) eq s1)
                       and IsSquare((KK!s1)*ratio[4]);
                if not okh then homoerr +:= 1;
                    printf "  HOMO-FAIL u=%o v=%o i=%o j=%o\n", u, v, i, j; end if;
                found := true; break;
            end if;
        end for;
        if not found then homoerr +:= 1; printf "  SUM-NOT-FOUND u=%o v=%o\n", u, v; end if;
    end for; end for;

    G, phi := TorsionSubgroup(J);
    invs := Invariants(G);

    nD := 0; nag := 0; ndis := 0; npos := 0; nsk := 0;
    for g in G do
        D := phi(g);
        if D eq J!0 then continue; end if;
        uD := D[1];
        got := false; dp := <Qq!1,Qq!1,Qq!1,KK!1>;
        if GCD(uD, f) eq 1 then
            got, dp := CO_DeltaPrimeCoprime(uD, roots, KK, th);
        else
            for k in [1..#T2] do
                D2 := D + T2[k][1];
                if D2 eq J!0 then got := true; dp := T2[k][2]; break; end if;
                if GCD(D2[1], f) eq 1 then
                    g2, d2 := CO_DeltaPrimeCoprime(D2[1], roots, KK, th);
                    if g2 then got := true; dp := CO_MulPrime(d2, T2[k][2]); break; end if;
                end if;
            end for;
        end if;
        if not got then nsk +:= 1; continue; end if;
        nD +:= 1;
        pred := CO_IsTrivial(dp);
        truth := IsDivisibleBy(D, 2);
        if pred eq truth then nag +:= 1;
        else ndis +:= 1;
            printf "  DISAGREE u=%o v=%o ord=%o uD=%o pred=%o truth=%o\n",
                   u, v, Order(D), uD, pred, truth;
        end if;
        if truth then npos +:= 1; end if;
    end for;
    totalD +:= nD; agree +:= nag; disagree +:= ndis; pos +:= npos; skipped +:= nsk;
    printf "SAMPLE u=%o v=%o torsion=%o classes=%o agree=%o disagree=%o divisible=%o skipped=%o\n",
           u, v, invs, nD, nag, ndis, npos, nsk;
end for;

printf "TOTAL samples=%o classes=%o agree=%o disagree=%o positives=%o skipped=%o homoerr=%o\n",
       tested, totalD, agree, disagree, pos, skipped, homoerr;
print "SEARCH_DONE";
quit;
