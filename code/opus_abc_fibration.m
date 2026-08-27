//////////////////////////////////////////////////////////////////////
// opus_abc_fibration.m  (2026-07-19)
//
// DEEP search on the (2,2,2,12) moduli surface via its elliptic
// fibration over the (A,B)-plane.  For fixed (A,B) the FIRST condition
//
//     y^2 = F1 = (A^2-B^2)(B^2-C^2)(A^2-B^2+C^2)
//
// is a quartic in C with the rational roots C = +-B, so it is an
// elliptic curve E_{A,B}.  Substituting C = B + 1/t clears it to a
// cubic:
//
//     v^2 = -D*(2Bt+1)*(A^2 t^2 + 2Bt + 1),   D = A^2-B^2,  v = y t^2.
//
// We compute the Mordell-Weil group of E_{A,B}, enumerate lattice
// points, recover C = B + 1/t for each, and test the SECOND condition
//     W*Q = square,  W = A^2-B^2+C^2,  Q = C^4 + D C^2 + A^2 D
// exactly.  Because C is allowed to be any RATIONAL number, this
// searches C of astronomically large height for each (A,B) -- far
// beyond any box scan (the projective point is (A*den : B*den : num)).
//
// Validation: (A,B)=(143,437) must return C=408 and C=1015, and
// (A,B)=(120,241) must return C=143 and C=266 (known curve orbits).
//
// Run: magma -b N:=30 K:=4 code/opus_abc_fibration.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(3*10^9);
SetClassGroupBounds("GRH");

if not assigned N then N := 30;
elif Type(N) eq MonStgElt then N := StringToInteger(N); end if;
if not assigned K then K := 4;
elif Type(K) eq MonStgElt then K := StringToInteger(K); end if;
if not assigned validate then validate := 1;
elif Type(validate) eq MonStgElt then validate := StringToInteger(validate); end if;

QQ := Rationals();
PR<t> := PolynomialRing(QQ);

// Jen's discriminant condition (nonzero <=> genus-2 curve is smooth)
function NonDegenerate(a, b, c)
    if a eq 0 or b eq 0 or c eq 0 then return false; end if;
    a2 := a^2; b2 := b^2; c2 := c^2;
    fs := [ b-c, b+c, a-b, a+b, a-c, a+c,
            -a2+b2-2*c2, -a2+b2-c2, -a2+b2-a*c-c2, -a2+b2+a*c-c2,
            -a2*b2+b2^2+2*a2*c2-2*b2*c2+c2^2,
            -a2^2+a2*b2-a2*c2+b2*c2-c2^2,
            a2^2-2*a2*b2+b2^2+2*a2*c2-b2*c2 ];
    for f in fs do if f eq 0 then return false; end if; end for;
    return true;
end function;

// second condition, exactly, for rational C
function SecondCondition(a, b, c)
    D := a^2 - b^2;
    W := D + c^2;
    if W eq 0 then return false; end if;
    Q := c^4 + D*c^2 + a^2*D;
    return IsSquare(W*Q);
end function;

// the fibre elliptic curve and the C-recovery
function FibreCurve(a, b)
    D := a^2 - b^2;
    g := -D*(2*b*t+1)*(a^2*t^2 + 2*b*t + 1);
    if Degree(g) ne 3 or Discriminant(g) eq 0 then return false, 0, 0; end if;
    Cg := HyperellipticCurve(g);
    ok := true;
    try
        pt := Cg![-1/(2*b), 0];
        E, mp := EllipticCurve(Cg, pt);
    catch e ok := false;
    end try;
    if not ok then return false, 0, 0; end if;
    return true, E, mp;
end function;

procedure ScanPair(a, b, K, ~nhit, ~nfail, ~report)
    okE, E, mp := FibreCurve(a, b);
    if not okE then nfail +:= 1; return; end if;
    okG := true;
    try
        G, gm := MordellWeilGroup(E);
    catch e okG := false; end try;
    if not okG then nfail +:= 1; return; end if;

    ngen := Ngens(G);
    if ngen eq 0 then return; end if;
    // enumerate coefficient vectors: torsion factors fully, free part in [-K,K]
    ranges := [];
    for i in [1..ngen] do
        oi := Order(G.i);
        if oi eq 0 then Append(~ranges, [-K..K]);
        else Append(~ranges, [0..oi-1]); end if;
    end for;
    // guard against absurd enumeration sizes
    tot := 1;
    for r in ranges do tot *:= #r; end for;
    if tot gt 40000 then
        for i in [1..ngen] do
            if Order(G.i) eq 0 then ranges[i] := [-2..2]; end if;
        end for;
    end if;

    minv := 0;
    try
        minv := Inverse(mp);
    catch e
        nfail +:= 1; return;
    end try;

    for vec in CartesianProduct(ranges) do
        cf := [vec[i] : i in [1..ngen]];
        if &and[c eq 0 : c in cf] then continue; end if;
        P := &+[ cf[i]*G.i : i in [1..ngen] ];
        pt := gm(P);
        if pt eq Identity(E) then continue; end if;
        okp := true;
        try
            qt := minv(pt);
            tv := qt[1]/qt[3];
        catch e okp := false; end try;
        if not okp or tv eq 0 then continue; end if;
        c := b + 1/tv;
        if c eq 0 then continue; end if;
        if not NonDegenerate(a, b, c) then continue; end if;
        if SecondCondition(a, b, c) then
            nhit +:= 1;
            num := Numerator(c); den := Denominator(c);
            Append(~report, <a, b, c, a*den, b*den, num>);
            printf "SURFACEPOINT A=%o B=%o C=%o  => projective (%o : %o : %o)\n",
                a, b, c, a*den, b*den, num;
        end if;
    end for;
end procedure;

// ---------- validation on the two known fibres ----------
if validate eq 1 then
    printf "=== VALIDATION ===\n";
    for pr in [<143,437>, <120,241>] do
        nh := 0; nf := 0; rep := [];
        ScanPair(pr[1], pr[2], 3, ~nh, ~nf, ~rep);
        cs := Sort([r[3] : r in rep]);
        printf "VALIDATE (A,B)=(%o,%o): %o surface points, C in %o\n",
            pr[1], pr[2], nh, cs;
    end for;
end if;

// ---------- production scan ----------
printf "=== SCAN N=%o K=%o ===\n", N, K;
total := 0; fails := 0; allrep := [];
for a in [1..N] do
    for b in [1..N] do
        if a eq b then continue; end if;
        if Gcd(a,b) ne 1 then continue; end if;
        nh := 0; nf := 0; rep := [];
        ScanPair(a, b, K, ~nh, ~nf, ~rep);
        total +:= nh; fails +:= nf;
        allrep cat:= rep;
    end for;
    printf "PROGRESS a=%o cumulative_points=%o fails=%o\n", a, total, fails;
end for;

printf "SCAN_DONE pairs_up_to=%o points=%o fails=%o\n", N, total, fails;
printf "DISTINCT (A,B,C) up to scaling:\n";
seen := {};
for r in allrep do
    g := Gcd([r[4], r[5], r[6]]);
    key := <r[4] div g, r[5] div g, r[6] div g>;
    if key notin seen then
        Include(~seen, key);
        printf "  (%o : %o : %o)\n", key[1], key[2], key[3];
    end if;
end for;
quit;
