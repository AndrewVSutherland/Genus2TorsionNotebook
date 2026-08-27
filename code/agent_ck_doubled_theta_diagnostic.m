//////////////////////////////////////////////////////////////////////
// Doubled-theta diagnostic on Elkies's Clebsch--Klein family.
//
// The split source is
//
//   y^2 = x * Product_i (x-r_i^2),
//   Sum_i r_i = Sum_i r_i^3 = 0.
//
// Sending x=0 to infinity by X=1/x gives the odd model
//
//   Z^2 = g(X) = Product_i (1-r_i^2*X).
//
// Elkies's marked atypical order-5 class is then
//
//   T = (0,1)-infinity = J![X,1].
//
// Let S01 be the 2-class {infinity,1/r_1^2} and S12 the
// 2-class {1/r_1^2,1/r_2^2}.  Up to sign and the S_5 action, the
// nontrivial order-10 classes which can still be doubled points are
//
//   T + S12,       2*T + S01,       2*T + S12.
//
// (For the omitted T+S01, the unique effective degree-two divisor is
// P+W_i with distinct support, so it cannot equal 2*R.)  For D=[u,v],
// the condition
//
//   D = 2*P-K = 2*(P-infinity)
//
// is exactly Degree(u)=2 and Discriminant(u)=0.  Every hit is
// independently reconstructed as Q=P-infinity and checked by 2*Q=D.
//
// Finite-field mode compares the three loci with the exact Stoll--
// Zarhin halving criteria already used in elkies22210_halving_covers.m.
// Rational mode completely enumerates primitive integral projective CK
// tuples of bounded coordinate height, removes square-set duplicates,
// and tests all five orbit-01 and all ten orbit-12 2-classes.
//
// Typical bounded runs:
//
//   magma -b mode:="finite" primes:="11,19,23,29,31,37" \
//       code/agent_ck_doubled_theta_diagnostic.m
//
//   magma -b mode:="rational" height:=30 \
//       code/agent_ck_doubled_theta_diagnostic.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(4*10^9);

Z := Integers();
Q := Rationals();

if not assigned mode then mode := "both"; end if;
require mode in {"finite","rational","both"}: "unknown mode";

if not assigned height then
    height := 30;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;

if not assigned max_print then
    max_print := 12;
elif Type(max_print) eq MonStgElt then
    max_print := StringToInteger(max_print);
end if;

if assigned primes then
    if Type(primes) eq MonStgElt then
        prime_list := [StringToInteger(a) : a in Split(primes, ",") | #a gt 0];
    else
        prime_list := primes;
    end if;
else
    prime_list := [11,19,23,29,31,37];
end if;
prime_list := [Z!p : p in prime_list];

function AbsGCD(vals)
    nz := [Abs(Z!v) : v in vals | v ne 0];
    return #nz eq 0 select Z!0 else GCD(nz);
end function;

function IsOpenTuple(rs)
    return (&and[r ne 0 : r in rs]) and
           (#Seqset(Setseq({r^2 : r in rs})) eq 5);
end function;

function IsCanonicalIntegralTuple(rs)
    if AbsGCD(rs) ne 1 or not IsOpenTuple(rs) then
        return false;
    end if;
    // Fix only global sign.  Permutations are removed by SquareKey.
    for r in rs do
        if r ne 0 then
            return r gt 0;
        end if;
    end for;
    return false;
end function;

function SquareKey(rs)
    return Join([IntegerToString(a) : a in Sort([Z!(r^2) : r in rs])], ",");
end function;

function OddCKPolynomial(K, PX, rs)
    X := PX.1;
    return &*[PX | 1-(K!r)^2*X : r in rs];
end function;

// Exact Stoll--Zarhin criterion for the 2-class represented by the
// unordered pair of entries i,j in branches.  This is copied in
// mathematical form, not loaded, so this file is independently runnable.
function PairHalves(branches, i, j)
    rem := [k : k in [1..6] | k ne i and k ne j];
    for k in rem do
        z := branches[k]-branches[j];
        for ell in rem do
            if ell ne k then
                z *:= branches[i]-branches[ell];
            end if;
        end for;
        if z eq 0 or not IsSquare(z) then
            return false;
        end if;
    end for;
    return true;
end function;

// If D lies on doubled theta, reconstruct the rational point Q=P-infinity.
// The returned alpha,beta satisfy P=(alpha,beta), 2*Q=D.
function DoubledThetaData(J, PX, g, D)
    uv := Eltseq(D);
    u := PX!uv[1];
    v := PX!uv[2];
    if Degree(u) ne 2 or Discriminant(u) ne 0 then
        return false, _, _, _;
    end if;
    lc := LeadingCoefficient(u);
    alpha := -Coefficient(u,1)/(2*lc);
    assert u eq lc*(PX.1-alpha)^2;
    beta := Evaluate(v,alpha);
    assert beta^2 eq Evaluate(g,alpha);
    QP := J![PX.1-alpha,beta];
    assert 2*QP eq D;
    return true, alpha, beta, QP;
end function;

procedure FixedMarkedClasses(K, PX, rs, J, ~T, ~S01, ~S12)
    X := PX.1;
    aa := [(K!r)^2 : r in rs];
    T := J![X,K!1];
    S01 := J![X-1/aa[1],K!0];
    S12 := J![(X-1/aa[1])*(X-1/aa[2]),K!0];
    O := J!0;
    assert T ne O and 5*T eq O;
    assert S01 ne O and 2*S01 eq O;
    assert S12 ne O and 2*S12 eq O;
end procedure;

procedure RunFinitePrime(p)
    assert IsPrime(p) and p notin {2,3,5};
    F := GF(p);
    PX<X> := PolynomialRing(F);
    PR<R4> := PolynomialRing(F);

    total_ck := 0;
    open_ck := 0;
    general01 := 0;
    general12 := 0;
    counts := AssociativeArray();
    counts["T_plus_S12"] := 0;
    counts["2T_plus_S01"] := 0;
    counts["2T_plus_S12"] := 0;
    samples := [];

    // Every open projective point has r1 nonzero; normalize it to 1.
    r1 := F!1;
    for r2 in F do
        for r3 in F do
            ss := r1+r2+r3;
            AA := r1^3+r2^3+r3^3;
            if ss eq 0 then
                r4s := [];
            else
                r4s := [z[1] : z in Roots(3*ss*R4^2+3*ss^2*R4+(ss^3-AA))];
            end if;
            for r4 in r4s do
                r5 := -(r1+r2+r3+r4);
                rs := [r1,r2,r3,r4,r5];
                if &+rs ne 0 or &+[r^3 : r in rs] ne 0 then
                    continue;
                end if;
                total_ck +:= 1;
                if not IsOpenTuple(rs) then
                    continue;
                end if;
                open_ck +:= 1;
                branches := [F!0] cat [r^2 : r in rs];
                half01 := PairHalves(branches,1,2);
                half12 := PairHalves(branches,2,3);
                if half01 then general01 +:= 1; end if;
                if half12 then general12 +:= 1; end if;

                g := OddCKPolynomial(F,PX,rs);
                assert Degree(g) eq 5 and Discriminant(g) ne 0 and
                       Evaluate(g,F!0) eq 1;
                J := Jacobian(HyperellipticCurve(g));
                FixedMarkedClasses(F,PX,rs,J,~T,~S01,~S12);
                cases := [
                    <"T_plus_S12",T+S12,S12,half12>,
                    <"2T_plus_S01",2*T+S01,S01,half01>,
                    <"2T_plus_S12",2*T+S12,S12,half12>
                ];
                for rec in cases do
                    label,D,S,cover_ok := Explode(rec);
                    assert Order(D) eq 10;
                    hit,alpha,beta,QP := DoubledThetaData(J,PX,g,D);
                    if not hit then
                        continue;
                    end if;
                    // A doubled-theta order-10 class forces the associated
                    // 2-class into 2J: 5D=S=2(5Q).
                    assert cover_ok and 5*D eq S and 2*(5*QP) eq S;
                    counts[label] +:= 1;
                    if #samples lt max_print then
                        Append(~samples,<label,[Z!r : r in rs],Z!alpha,Z!beta>);
                    end if;
                end for;
            end for;
        end for;
    end for;

    print "FINITE", "p",p,"total_ck",total_ck,"open_ck",open_ck,
          "general_halving_01",general01,"general_halving_12",general12,
          "double_T_plus_S12",counts["T_plus_S12"],
          "double_2T_plus_S01",counts["2T_plus_S01"],
          "double_2T_plus_S12",counts["2T_plus_S12"];
    if #samples gt 0 then
        print "FINITE_SAMPLES",p,samples;
    end if;
end procedure;

procedure AuditRationalSource(rs, source_index,
                              ~classes_tested, ~hits, ~records)
    PX<X> := PolynomialRing(Q);
    g := OddCKPolynomial(Q,PX,rs);
    assert Degree(g) eq 5 and Discriminant(g) ne 0 and Evaluate(g,Q!0) eq 1;
    J := Jacobian(HyperellipticCurve(g));
    T := J![X,Q!1];
    assert Order(T) eq 5;
    aa := [(Q!r)^2 : r in rs];

    // Case 2T+S01: all five orbit-01 classes.
    for i in [1..5] do
        S := J![X-1/aa[i],Q!0];
        D := 2*T+S;
        assert Order(S) eq 2 and Order(D) eq 10;
        classes_tested +:= 1;
        hit,alpha,beta,QP := DoubledThetaData(J,PX,g,D);
        if not hit then continue; end if;
        branches := [Q!0] cat aa;
        assert PairHalves(branches,1,i+1);
        assert 5*D eq S and 2*(5*QP) eq S and Order(QP) eq 20;
        invs := Invariants(TorsionSubgroup(J));
        assert #invs ge 4 and invs[#invs] mod 20 eq 0;
        hits +:= 1;
        Append(~records,<source_index,"2T_plus_S01",[0,i],rs,
                         alpha,beta,invs>);
    end for;

    // Cases T+S12 and 2T+S12: all ten orbit-12 classes.
    for i in [1..4] do
        for j in [i+1..5] do
            S := J![(X-1/aa[i])*(X-1/aa[j]),Q!0];
            assert Order(S) eq 2;
            branches := [Q!0] cat aa;
            cover_ok := PairHalves(branches,i+1,j+1);
            for rec in [<"T_plus_S12",T+S>,<"2T_plus_S12",2*T+S>] do
                label,D := Explode(rec);
                assert Order(D) eq 10;
                classes_tested +:= 1;
                hit,alpha,beta,QP := DoubledThetaData(J,PX,g,D);
                if not hit then continue; end if;
                assert cover_ok and 5*D eq S and 2*(5*QP) eq S and
                       Order(QP) eq 20;
                invs := Invariants(TorsionSubgroup(J));
                assert #invs ge 4 and invs[#invs] mod 20 eq 0;
                hits +:= 1;
                Append(~records,<source_index,label,[i,j],rs,
                                 alpha,beta,invs>);
            end for;
        end for;
    end for;
end procedure;

procedure RunRational(H)
    assert H ge 1;
    seen := {};
    triple_checked := 0;
    tuple_checked := 0;
    ck_open_canonical := 0;
    source_count := 0;
    classes_tested := 0;
    hits := 0;
    records := [];

    for r1 in [-H..H] do
        for r2 in [-H..H] do
            for r3 in [-H..H] do
                triple_checked +:= 1;
                ss := r1+r2+r3;
                if ss eq 0 then continue; end if;
                AA := r1^3+r2^3+r3^3;
                disc := 3*ss*(4*AA-ss^3);
                is_square,root := IsSquare(disc);
                if not is_square then continue; end if;

                for root_signed in Sort(Setseq({-root,root})) do
                    q4 := Q!(-3*ss^2+root_signed)/(6*ss);
                    if Denominator(q4) ne 1 then continue; end if;
                    r4 := Z!q4;
                    r5 := -(r1+r2+r3+r4);
                    if Abs(r4) gt H or Abs(r5) gt H then continue; end if;
                    rs := [r1,r2,r3,r4,r5];
                    tuple_checked +:= 1;
                    if not IsCanonicalIntegralTuple(rs) then continue; end if;
                    ck_open_canonical +:= 1;
                    key := SquareKey(rs);
                    if key in seen then continue; end if;
                    Include(~seen,key);
                    source_count +:= 1;
                    AuditRationalSource(rs,source_count,
                                        ~classes_tested,~hits,~records);
                end for;
            end for;
        end for;
    end for;

    print "RATIONAL", "height",H,"triple_checked",triple_checked,
          "tuple_checked",tuple_checked,
          "open_canonical_tuples",ck_open_canonical,
          "unique_square_sources",source_count,
          "doubled_theta_classes_tested",classes_tested,"hits",hits;
    if #records gt 0 then
        print "RATIONAL_HITS",records;
    end if;
end procedure;

print "AGENT_CK_DOUBLED_THETA_DIAGNOSTIC";
print "mode",mode,"memory_limit_bytes",4*10^9;
print "cases",["T_plus_S12","2T_plus_S01","2T_plus_S12"];

if mode in {"finite","both"} then
    print "FINITE_PRIMES",prime_list;
    for p in prime_list do RunFinitePrime(p); end for;
end if;

if mode in {"rational","both"} then
    RunRational(height);
end if;

print "AGENT_CK_DOUBLED_THETA_DONE";
quit;
