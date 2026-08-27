//////////////////////////////////////////////////////////////////////
//  claude_ov_2216_isoclass.m   -- lane OV/2216
//
//  How many DISTINCT curves do the norm-surface hits actually represent?
//
//  hits2uv.py deduplicates only by the S_3 that permutes the root triple
//  {u,v,w}.  But the M_1(8,2,2) chart also depends on WHICH of the 16
//  order-8 classes is marked, so one curve can occur at many (u,v).  This
//  script reduces the candidate list to Q-isomorphism classes and prints,
//  for each class, its members, its Igusa-Clebsch invariants, the square
//  class of the discriminant DK of the quadratic Weierstrass pair, and a
//  Frobenius fingerprint (a_p for the first good primes) as a cross-check.
//
//  Usage:
//    code/claude_magma_slot.sh -b candidate_file:=data/xxx.txt \
//        code/claude_ov_2216_isoclass.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned candidate_file then
    candidate_file := "data/claude_ov_2216_all15_uv.txt";
end if;
if not assigned MemGB then MemGB := 6;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Qq := Rationals();
P<X> := PolynomialRing(Qq);

function ParseQ(s)
    if "/" in s then
        pr := Split(s, "/");
        return Qq!StringToInteger(pr[1]) / Qq!StringToInteger(pr[2]);
    end if;
    return Qq!StringToInteger(s);
end function;

function SqClass(a)
    n := Numerator(a); d := Denominator(a);
    m := n*d; s := Sign(m); m := AbsoluteValue(m); r := 1;
    for pe in Factorization(m) do
        if IsOdd(pe[2]) then r *:= pe[1]; end if;
    end for;
    return s*r;
end function;

lines := Split(Read(candidate_file), "\n");
printf "claude_ov_2216_isoclass  candidate_file=%o\n", candidate_file;

curves := []; labels := []; dks := [];
for raw in lines do
    s := raw;
    if #s gt 0 and s[#s] eq "\r" then s := s[1..#s-1]; end if;
    if #s eq 0 or s[1] eq "#" then continue; end if;
    prt := [t : t in Split(s, " ") | #t gt 0];
    if #prt lt 2 then continue; end if;
    u := ParseQ(prt[1]);  v := ParseQ(prt[2]);
    Bc := u^2*v - u^2 + u*v^2 - u - v^2 - v - 2;
    Cc := u^2 + u*v + v^2 + u + v + 1;
    qtm := X^2 - Bc*X + Cc;
    f := ((1-u)*X + 1)*((1-v)*X + 1)*((u+v+2)*X + 1)*(-qtm);
    L := 1;
    for i in [0..5] do L := LCM(L, Denominator(Coefficient(f, i))); end for;
    fI := P!(L^2*f);
    Append(~curves, HyperellipticCurve(fI));
    Append(~labels, Sprintf("(%o, %o)", u, v));
    Append(~dks, SqClass(Bc^2 - 4*Cc));
end for;
printf "candidates read: %o\n\n", #curves;

reps := []; members := []; repdk := [];
for i in [1..#curves] do
    hit := 0;
    for j in [1..#reps] do
        if IsIsomorphic(curves[i], curves[reps[j]]) then hit := j; break; end if;
    end for;
    if hit eq 0 then
        Append(~reps, i); Append(~members, [i]); Append(~repdk, dks[i]);
    else
        Append(~members[hit], i);
    end if;
end for;

printf "DISTINCT_Q_ISOMORPHISM_CLASSES %o   (from %o candidate (u,v) points)\n\n",
       #reps, #curves;

// good primes for ALL representatives simultaneously
goodp := [];
p := 3;
while #goodp lt 12 and p lt 400 do
    ok := true;
    for i in [1..#reps] do
        fp := HyperellipticPolynomials(curves[reps[i]]);
        num := LCM([ Denominator(cc) : cc in Coefficients(fp) ]);
        fz := PolynomialRing(Integers())!(num*fp);
        if Valuation(Integers()!LeadingCoefficient(fz), p) gt 0
           or Valuation(Integers()!Discriminant(fz), p) gt 0
           or Valuation(num, p) gt 0 or p eq 2 then
            ok := false; break;
        end if;
    end for;
    if ok then Append(~goodp, p); end if;
    p := NextPrime(p);
end while;

for j in [1..#reps] do
    i := reps[j];
    C := curves[i];
    Cmin := ReducedModel(SimplifiedModel(C));
    IC := IgusaClebschInvariants(C);
    wt := [1,2,3,5];
    // normalise the weighted-projective IC tuple to a canonical rational tuple
    nz := 0;
    for t in [1..4] do if IC[t] ne 0 then nz := t; break; end if; end for;
    normIC := [ IC[t]^(wt[nz]) / IC[nz]^(wt[t]) : t in [1..4] ];
    aps := []; v2s := [];
    for pp in goodp do
        Np := #Jacobian(ChangeRing(C, GF(pp)));
        Append(~aps, Np);
        Append(~v2s, Valuation(Np, 2));
    end for;
    minv2 := Min(v2s);
    kills := [ goodp[t] : t in [1..#goodp] | v2s[t] lt 6 ];
    printf "CLASS %o  DK_sqclass=%o\n", j, dks[i];
    printf "   members (%o): %o\n", #members[j], [ labels[m] : m in members[j] ];
    printf "   reduced model y^2 = %o\n", HyperellipticPolynomials(Cmin);
    printf "   normalised Igusa-Clebsch = %o\n", normIC;
    printf "   #J(F_p) for p in %o : %o\n", goodp, aps;
    printf "   v_2(#J(F_p))            : %o   min=%o\n", v2s, minv2;
    // [2,2,16] has order 64; if v_2(#J(F_p)) < 6 at a good prime it is impossible
    printf "   FROBENIUS_KILLS_2216 %o   witnesses p = %o\n", #kills gt 0, kills;
end for;

print "SEARCH_DONE";
quit;
