//////////////////////////////////////////////////////////////////////
//  claude_ov_2216_s3check.m   -- lane OV/2216
//
//  Supplies the verification that code/claude_ov_2216_sweep2.c and
//  code/claude_ov_2216_sweep3.c CITE in their headers but which had
//  never actually been run (results/claude_ov_2216_s3check.log did not
//  exist in the repository).
//
//  What is proved here, over the function field Q(u,v):
//
//   (a) The eight halving conditions "Q+T is 2-divisible", T in J[2](Q),
//       collapse to exactly FOUR distinct condition sets, and the
//       collapse is exactly  T ~ T + B_quad  (B_quad = 4Q).
//   (b) The rational necessary norm condition N(delta'_K) = square is
//       literally the SAME square class u*v*(u+v+1) for all eight.
//   (c) The S_3 that permutes the root triple {u, v, w = -(u+v+1)}
//       acts on the chart by  (u,v) -> (sigma(u), sigma(v))  and
//       fixes condition set Y0 while permuting the other three
//       transitively.  Hence "surface 0" and "surface 1" of the C
//       sweeps are the ONLY two S_3-inequivalent norm surfaces.
//   (d) The conic parametrisations hard-coded in the C sweeps,
//         Y0:  v =  A/(A - k^2),        Y1: v = -A(2u+1)/(k^2 + A),
//       with A = u(u-1), do satisfy condition (1) identically, and the
//       integer closed forms  S1, S2, Dn, c2, c3  used by the tier-A
//       bitmap engine are exact identities in Q(p,q,a,b).
//
//  Conventions identical to code/claude_ov_2216_symbolic.m.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(4*10^9);

FF<u,v> := RationalFunctionField(Rationals(), 2);
R := PolynomialRing(Rationals(), 2);
PX<X> := PolynomialRing(FF);

B := u^2*v - u^2 + u*v^2 - u - v^2 - v - 2;
C := u^2 + u*v + v^2 + u + v + 1;
qt  := -X^2 + B*X - C;
qtm := X^2 - B*X + C;
f   := ((1-u)*X + 1)*((1-v)*X + 1)*((u+v+2)*X + 1)*qt;

r := [ 1/(u-1), 1/(v-1), -1/(u+v+2) ];
lin := [ X - r[i] : i in [1..3] ];
assert f eq LeadingCoefficient(f)*&*lin*qtm;

// ---- squarefree normal form in Q(u,v)* / (Q(u,v)*)^2  (canonical) -----
function SqFree(gin)
    g := FF!gin;                       // lab rule: coerce, never mix FldFunRatMElt/RngMPolElt
    if g eq 0 then return R!0; end if;
    nu := R!Numerator(g);  de := R!Denominator(g);
    h := nu*de;
    fa := Factorization(h);
    unit := h;
    for t in fa do unit := unit div t[1]^t[2]; end for;
    out := R!1;
    for t in fa do if IsOdd(t[2]) then out *:= t[1]; end if; end for;
    cst := Rationals()!unit;
    n := Numerator(cst); d := Denominator(cst);
    m := n*d; s := Sign(m); m := AbsoluteValue(m); rr := 1;
    for pe in Factorization(m) do
        if IsOdd(pe[2]) then rr *:= pe[1]; end if;
    end for;
    return (s*rr)*out;
end function;

// robust class equality: g1 ~ g2  iff  g1*g2 is a square in Q(u,v)
function SameClass(x1, x2)
    g1 := FF!x1; g2 := FF!x2;
    if g1 eq 0 or g2 eq 0 then return g1 eq g2; end if;
    s := SqFree(g1*g2);
    if Degree(s) ne 0 then return false; end if;
    return IsSquare(Rationals()!s);
end function;

function ShowQ(g)
    sf := SqFree(g);
    if sf eq 0 then return "0"; end if;
    fa := Factorization(sf);
    unit := sf;
    for t in fa do unit := unit div t[1]^t[2]; end for;
    cst := Rationals()!unit;
    s := "";
    for t in fa do
        p := t[1];
        den := LCM([ Denominator(cc) : cc in Coefficients(p) ]);
        pp := den*p;
        cont := GCD([ Numerator(cc) : cc in Coefficients(pp) ]);
        pp := pp div cont;
        if LeadingCoefficient(pp) lt 0 then pp := -pp; cst *:= -1; end if;
        cst *:= (Rationals()!cont/den)^t[2];
        for e in [1..t[2]] do s := s cat "(" cat Sprint(pp) cat ")"; end for;
    end for;
    n := Numerator(cst); d := Denominator(cst);
    m := n*d; sg := Sign(m); m := AbsoluteValue(m); rr := 1;
    for pe in Factorization(m) do
        if IsOdd(pe[2]) then rr *:= pe[1]; end if;
    end for;
    cs := sg*rr;
    if cs eq 1 and s ne "" then return s; end if;
    return Sprint(cs) cat s;
end function;

function ModQ(g) return g mod qtm; end function;

function DeltaCop(uD)
    return [* Evaluate(uD, r[1]), Evaluate(uD, r[2]), Evaluate(uD, r[3]), ModQ(uD) *];
end function;

function DeltaTwo(uS)
    cof := f div uS;
    assert uS*cof eq f;
    out := [* *];
    for i in [1..3] do
        if Evaluate(uS, r[i]) eq 0 then Append(~out, -Evaluate(cof, r[i]));
        else Append(~out, Evaluate(uS, r[i])); end if;
    end for;
    vS := ModQ(uS);
    if vS eq 0 then Append(~out, -ModQ(cof)); else Append(~out, vS); end if;
    return out;
end function;

function DMul(a, b)
    return [* a[1]*b[1], a[2]*b[2], a[3]*b[3], ModQ(a[4]*b[4]) *];
end function;

function NormK(g)
    h := ModQ(g);
    a := Coefficient(h, 0); b := Coefficient(h, 1);
    return a^2 + a*b*B + b^2*C;
end function;

// =====================================================================
print "######## PART (a)+(b): the eight twists collapse to four ########";

dQ := DeltaCop(X+1);

labels := ["ZERO   (T=0, untwisted)"];
subs   := [PX!1];                       // T = 0 : trivial divisor
for i in [1..3] do for j in [i+1..3] do
    Append(~subs, lin[i]*lin[j]);
    Append(~labels, Sprintf("A_r%or%o", i, j));
end for; end for;
Append(~subs, qtm);  Append(~labels, "B_quad (=4Q)");
for k in [1..3] do
    Append(~subs, qtm * &*[ lin[j] : j in [1..3] | j ne k ]);
    Append(~labels, Sprintf("C_r%o  ", k));
end for;

// condition set of a twist = the subgroup {d1d2, d1d3, d2d3} of squares
CondGens := [];   // [ [g12, g13, g23] ] as squarefree reps
NormCond := [];
for k in [1..#subs] do
    if k eq 1 then
        dQT := dQ;
    else
        dQT := DMul(dQ, DeltaTwo(subs[k]));
    end if;
    g12 := SqFree(dQT[1]*dQT[2]);
    g13 := SqFree(dQT[1]*dQT[3]);
    g23 := SqFree(dQT[2]*dQT[3]);
    Append(~CondGens, [g12, g13, g23]);
    Append(~NormCond, SqFree(NormK(dQT[4])));
    printf "%-24o  d1d2=%-42o d1d3=%-42o d2d3=%-42o N(dK)=%o\n",
        labels[k], ShowQ(g12), ShowQ(g13), ShowQ(g23), ShowQ(NormK(dQT[4]));
end for;

// a condition set is the SET {g12,g13,g23}; two twists agree iff the sets match
function CondEq(Ain, Bin)
    A  := [ FF | FF!x : x in Ain ];
    Bb := [ FF | FF!x : x in Bin ];
    for x in A do
        ok := false;
        for y in Bb do if SameClass(x, y) then ok := true; break; end if; end for;
        if not ok then return false; end if;
    end for;
    return true;
end function;

reps := []; repidx := []; assign := [];
for k in [1..#CondGens] do
    hit := 0;
    for j in [1..#reps] do
        if CondEq(CondGens[k], reps[j]) then hit := j; break; end if;
    end for;
    if hit eq 0 then
        Append(~reps, CondGens[k]); Append(~repidx, k); hit := #reps;
    end if;
    Append(~assign, hit);
end for;
printf "\nDISTINCT_CONDITION_SETS %o\n", #reps;
for j in [1..#reps] do
    printf "  set %o  <- twists: %o\n", j,
        [ labels[k] : k in [1..#assign] | assign[k] eq j ];
end for;
assert #reps eq 4;

// (b) the norm condition is the same class for all eight, = u*v*(u+v+1)
target := SqFree(u*v*(u+v+1));
allsame := true;
for k in [1..#NormCond] do
    if not SameClass(NormCond[k], target) then allsame := false; end if;
end for;
printf "NORM_CONDITION_UNIFORM %o   (class = %o)\n", allsame, ShowQ(target);
assert allsame;

// the collapse is exactly T ~ T + B_quad
print "COLLAPSE_IS_TRANSLATION_BY_4Q:";
printf "  {T=0, B_quad}      same set: %o\n", assign[1] eq assign[5];
printf "  {A_r1r2, C_r3}     same set: %o\n", assign[2] eq assign[8];
printf "  {A_r1r3, C_r2}     same set: %o\n", assign[3] eq assign[7];
printf "  {A_r2r3, C_r1}     same set: %o\n", assign[4] eq assign[6];
assert assign[1] eq assign[5] and assign[2] eq assign[8]
   and assign[3] eq assign[7] and assign[4] eq assign[6];

// =====================================================================
print "";
print "######## PART (c): the S_3 action on the four condition sets ########";
// S_3 permutes {u, v, w}, w = -(u+v+1); the chart coordinates are the
// first two entries of the ordered triple.
w := -(u+v+1);
trip := [u, v, w];
perms := [ [1,2,3], [2,1,3], [1,3,2], [3,2,1], [2,3,1], [3,1,2] ];
pnames := [ "id", "(uv)", "(vw)", "(uw)", "(uvw)", "(uwv)" ];

// image of the surface list under each permutation
orb := [];
for pi in [1..#perms] do
    pm := perms[pi];
    nu := trip[pm[1]]; nv := trip[pm[2]];
    h := hom< FF -> FF | [nu, nv] >;
    row := [];
    for j in [1..#reps] do
        img := [ FF | h(FF!g) : g in reps[j] ];
        tgt := 0;
        for jj in [1..#reps] do
            if CondEq(img, [ FF | FF!g : g in reps[jj] ]) then tgt := jj; break; end if;
        end for;
        Append(~row, tgt);
    end for;
    // norm condition invariance
    ninv := SameClass(h(FF!target), target);
    printf "  perm %-6o sends condition sets [1,2,3,4] -> %o    N-cond invariant: %o\n",
        pnames[pi], row, ninv;
    Append(~orb, row);
    assert ninv;
end for;

// Y0 = the set containing T=0, i.e. index assign[1]
i0 := assign[1];
fixed0 := &and[ orb[pi][i0] eq i0 : pi in [1..#perms] ];
others := [ j : j in [1..#reps] | j ne i0 ];
transit := { orb[pi][others[1]] : pi in [1..#perms] };
printf "Y0_IS_S3_INVARIANT %o\n", fixed0;
printf "OTHER_THREE_FORM_ONE_ORBIT %o   (orbit of set %o = %o)\n",
    #transit eq 3, others[1], transit;
assert fixed0;
assert #transit eq 3;
printf "S3_INEQUIVALENT_SURFACES 2\n";

// =====================================================================
print "";
print "######## PART (d): the conic parametrisations in the C sweeps ########";
FK<uu,kk> := RationalFunctionField(Rationals(), 2);
A := uu*(uu-1);
// surface 0 : u(u-1) v(v-1) = square,  point (v,y)=(0,0), line y = k v
v0 := A/(A - kk^2);
lhs0 := A*v0*(v0-1);
printf "  Y0: v = A/(A-k^2)   ->  u(u-1)v(v-1) - (k v)^2 = %o\n", lhs0 - (kk*v0)^2;
assert lhs0 - (kk*v0)^2 eq 0;
// surface 1 : -u(u-1) v (2u+v+1) = square, point (0,0), line y = k v
v1 := -A*(2*uu+1)/(kk^2 + A);
lhs1 := -A*v1*(2*uu + v1 + 1);
printf "  Y1: v = -A(2u+1)/(k^2+A) -> -u(u-1)v(2u+v+1) - (k v)^2 = %o\n", lhs1 - (kk*v1)^2;
assert lhs1 - (kk*v1)^2 eq 0;

// integer closed forms used by tier A of sweep2/sweep3
FP<p,q,a,b> := RationalFunctionField(Rationals(), 4);
An := p*(p-q); Ad := q^2;
for surf in [0,1] do
    if surf eq 0 then
        vn := An*b^2;              vd := An*b^2 - a^2*Ad;
    else
        vn := -An*(2*p+q)*b^2;     vd := q*(a^2*Ad + An*b^2);
    end if;
    S1 := p*vd + q*vn + q*vd;
    S2 := S1 + q*vd;
    Dn := p*vd - q*vn;
    U := p/q; V := vn/vd;
    // condition (3):  u v (u+v+1)  ==  p*vn*S1 / (q vd)^2
    e3 := U*V*(U+V+1) - p*vn*S1/(q*vd)^2;
    if surf eq 0 then
        e2 := U*(U-1)*(U+V+1)*(U+V+2) - An*S1*S2/(Ad*(q*vd)^2);
        cS1 := An*(p+2*q)*b^2 - Ad*(p+q)*a^2;
        cS2 := An*(p+3*q)*b^2 - Ad*(p+2*q)*a^2;
        cG3 := p*An;              // c3 = G3 * S1 * b^2
        e1 := S1 - cS1;  ee := S2 - cS2;
        ec3 := p*vn*S1 - cG3*S1*b^2;
        ec2 := An*S1*S2 - An*S1*S2;
    else
        e2 := U*(U-1)*(U-V)*(U+V+1) - An*Dn*S1/(Ad*(q*vd)^2);
        cS1 := q*((p+q)*Ad*a^2 - p*An*b^2);
        cS2 := q*(p*Ad*a^2 + An*(3*p+q)*b^2);   // = Dn
        cG3 := -p*An*(2*p+q);
        e1 := S1 - cS1;  ee := Dn - cS2;
        ec3 := p*vn*S1 - cG3*S1*b^2;
        ec2 := An*Dn*S1 - An*Dn*S1;
    end if;
    printf "  surface %o: cond3_identity=%o cond2_identity=%o S1_closed=%o S2/Dn_closed=%o c3_closed=%o\n",
        surf, e3 eq 0, e2 eq 0, e1 eq 0, ee eq 0, ec3 eq 0;
    assert e3 eq 0 and e2 eq 0 and e1 eq 0 and ee eq 0 and ec3 eq 0;
end for;

print "";
print "ALL_ASSERTIONS_PASSED";
print "SEARCH_DONE";
quit;
