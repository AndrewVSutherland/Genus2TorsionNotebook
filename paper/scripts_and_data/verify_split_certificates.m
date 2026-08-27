// verify_split_certificates.m -- geometric splitness certificates for ALL
// 77 rows of Table 2 (tab:splitcensus) of the paper.  The certified
// curves are read from table2.txt (the machine-readable table data from
// which table2.tex is typeset); run from this directory:
//   magma -b verify_split_certificates.m
//
// Routes (established once by an exploration over the (2,2)-isogeny graph):
//  R1: RichelotIsogenousSurfaces(Jac) over Q contains a DEGENERATE codomain
//      -- either a product E x E' of elliptic curves over Q (SetCart; the
//      conductors are printed) or an elliptic curve over a quadratic field
//      K (CrvEll; the codomain is then E x E^sigma over K).  Either way Jac
//      admits an explicit (2,2)-isogeny to a product over a number field,
//      so it is geometrically split.  59 of the 77 rows certify this way
//      (56 over Q, one to E x E^sigma over Q(sqrt 11) [the [11] row],
//      2 more over Q(sqrt 2) at depth <= 2: route Q2).
//  COVER: ALL 18 odd-glued curves carry an EXPLICIT nonconstant
//      map phi = (p(x)/q(x), y h(x)/q(x)^2) from the genus-2 curve to a
//      genus-1 curve E, over Q (13 rows, E given by its Cremona label) or
//      over an explicitly given number field K of degree 2 or 3 (4 rows,
//      E : v^2 = 4u^3 - g2 u - g3 with g2, g3 in K).  The verification is a
//      single polynomial identity N(p,q)*q = g*h^2 (plus nonvanishing of
//      the Wronskian p q' - p' q and of disc(E)), checked exactly below; a
//      nonconstant morphism C -> E over any field induces a surjection
//      Jac(C) ->> E, so Jac(C) is geometrically split.  The maps were FOUND
//      numerically (period matrices, endomorphism idempotents, Weierstrass-P
//      evaluation, exact rational/algebraic reconstruction), but the
//      certificate is independent of the search: only the identity matters.
//      All 13 rational covers have degree 3; the K-covers have degrees
//      2, 2, 2, 3, 5 -- in particular the [70] curve of Howe maps to 858k1
//      and [19] = J_1(13) is bielliptic over an explicit sextic field.
//  SIG: (retained for completeness; no row currently uses it): no chain of
//      (2,2)-isogenies can reach a product (the products in their isogeny
//      class are odd-isogenous), so no Richelot-style certificate exists.
//      For these we verify the split SIGNATURE: at EVERY good prime
//      p < 200 the Frobenius polynomial fails the root-power simplicity
//      test (chi_p is reducible, or some power pi^n with n <= 12 generates
//      a proper subfield).  A single strict prime would PROVE geometric
//      simplicity (Part A of verify_simple_certificates.m), so for a
//      genuinely split Jacobian this failure occurs at every good prime;
//      note some of these curves (e.g. [19] = J_1(13)) are Q-simple, so
//      chi_p itself is often irreducible -- the power-drop is the geometric
//      fingerprint.  The splitness statement itself is cited: the certified
//      LMFDB/database endomorphism computation (Costa-Mascot-Sijsling-
//      Voight) for database rows, and Howe 2015 (Thm 2.4) for the [70]
//      curve.  (The signature is strong evidence, not a proof: finitely
//      many Frobenius polynomials cannot prove an endomorphism exists.)
SetColumns(0);
SetMemoryLimit(8*10^9);
Q := Rationals(); Z := Integers();
R<x> := PolynomialRing(Q);
Tp<T> := PolynomialRing(Q);
PxQ<X> := PolynomialRing(Q);

// root-power strictness (as in verify_simple_certificates.m)
function IsStrict(chi)
    if Degree(chi) ne 4 or not IsIrreducible(chi) then return false; end if;
    K := NumberField(chi); pi := K.1;
    for n in [2..12] do
        if Degree(MinimalPolynomial(pi^n)) lt 4 then return false; end if;
    end for;
    return true;
end function;

function SplitNode(L)
    for s in L do
        if Type(s) eq SetCart or Type(s) eq CrvEll then return true, s; end if;
    end for;
    return false, 0;
end function;

// The curves AND their certificate routes are read from table2.txt (the
// machine-readable table data from which table2.tex is typeset by
// make_tables.py), so certificate and table cannot drift.  The route
// column is R1, Q2, or COVER:n, with n indexing the covers/coverK data
// below; the tag printed with each certificate is the curve itself,
// [[f],[h]] coefficient lists.
rows := [* *];
for line in Split(Read("table2.txt"), "\n") do
    if #line eq 0 or line[1] eq "#" then continue; end if;
    parts := Split(line, "|");
    assert #parts eq 7;   // group | [[f],[h]] | label | display | source | route | comment
    fh := eval parts[2];
    rt := Split(parts[6], ":");
    ck := #rt eq 2 select StringToInteger(rt[2]) else 0;
    Append(~rows, <eval parts[1], fh[1], fh[2], rt[1], parts[2], ck>);
end for;
assert #rows eq 77;

covers := AssociativeArray();
coverK := AssociativeArray();
K := NumberField(X^6 + 2925/2*X^5 - 13931515/144*X^4 + 555105005/864*X^3 + 741300755/1152*X^2 - 59482995065/62208*X - 2129491223327/2985984);
PK := PolynomialRing(K);
covers[1] := < "K", "", 2, PK![ K | 1/350499166447916498915753838*(-1047400442066319243264*K.1^5 - 1530244557357870304783296*K.1^4 + 103639941866199892102086384*K.1^3 - 827402670280905916142666448*K.1^2 + 497327212158924622832114979*K.1 + 825849916525810185225936473), 1/175249583223958249457876919*(1181047455487212756480*K.1^5 + 1725717953028064158570240*K.1^4 - 116547640063212079646312544*K.1^3 + 913046320206876174333647256*K.1^2 - 421538903357027186361600372*K.1 - 655021524942587374393280381), 1/1401996665791665995663015352*(15990842261263127814144*K.1^5 + 23367448923609665432475648*K.1^4 - 1575054550234548111640212288*K.1^3 + 12162386878499496157273057680*K.1^2 - 4320129434415791253147032064*K.1 - 8711193399793278371297110063) ], PK![ K | 1/1518829721274304828634933298*(62752690556056115767296*K.1^5 + 91701242410418949541259520*K.1^4 - 6180068925644118438615566976*K.1^3 + 47674259218659743345753384112*K.1^2 - 17057854345563835128788883108*K.1 - 34299957341643443052026541793), 1/1518829721274304828634933298*(16541103711402026827776*K.1^5 + 24166552645798839366764544*K.1^4 - 1636558488650062614742432704*K.1^3 + 13063575221221494738458433744*K.1^2 - 7641246715946420125902252744*K.1 - 11882728243714815255723273527), 1 ], PK![ K | 1/467332221930555331887671784*(-263143200386260062268416*K.1^5 - 384520483678913448260359680*K.1^4 + 25935239216690935277561086464*K.1^3 - 201243384412403003124237441216*K.1^2 + 80517482866691594749291161324*K.1 + 150570853786333685914255361549), 1/155777407310185110629223928*(-12657327167462518499328*K.1^5 - 18495495027938281568656896*K.1^4 + 1247708735528258120286326976*K.1^3 - 9694015543269686983981601520*K.1^2 + 3966652048288737212257203624*K.1 + 7357662702216901913155444543) ], PK!0, PK!0 >;
coverK[1] := < K, K!(K.1), K!(1/4205989997374997986989046056*(-272112803100608279040*K.1^5 - 395664740768467375261056*K.1^4 + 29647589653215861439432320*K.1^3 - 477245639866772050871629032*K.1^2 - 195733771623786479752709040*K.1 - 207843556759852706226835301)) >;
K := NumberField(X^3 - 5175/4*X^2 + 50625/8*X - 421875/64);
PK := PolynomialRing(K);
covers[2] := < "K", "", 2, PK![ K | 1/397500*(64*K.1^2 + 16000*K.1 + 95625), 1/99375*(24*K.1^2 - 15200*K.1 + 73125), 1/397500*(-16*K.1^2 + 64900*K.1 + 50625) ], PK![ K | 1/59625*(16*K.1^2 - 19320*K.1 + 68625), 1/298125*(256*K.1^2 - 328200*K.1 + 978750), 1 ], PK![ K | 1/19875*(-224*K.1^2 + 3360*K.1 - 6750), 1/33125*(-64*K.1^2 - 5400*K.1 + 3750) ], PK!0, PK!0 >;
coverK[2] := < K, K!(K.1), K!(1/795000*(-4192*K.1^2 - 353700*K.1 + 245625)) >;
covers[3] := < "Q", "26b1", 3, PxQ![-144, -36, -144, -45], PxQ![0, 4, 0, 1], PxQ![-864, -432, 0, 432], PxQ!0, PxQ!0 >;
covers[4] := < "Q", "20a2", 3, PxQ![36, 12, 60, 12], PxQ![0, 4, 5, 1], PxQ![216, 648, 216], PxQ!0, PxQ!0 >;
covers[5] := < "Q", "19a1", 3, PxQ![-84, 684, -276, 192], PxQ![-4, 0, -5, 1], PxQ![4104, -4104, 0, 2052], PxQ!0, PxQ!0 >;
covers[6] := < "Q", "35a1", 3, PxQ![-36, -12, 60, 48], PxQ![-6, 5, -4, 1], PxQ![-756, 3780, -2268, 756], PxQ!0, PxQ!0 >;
covers[7] := < "Q", "11a1", 3, PxQ![-156, 2100, -1536, 564], PxQ![-8, -4, 5, 1], PxQ![-26136, 52272, -52272, 13068], PxQ!0, PxQ!0 >;
covers[8] := < "Q", "90c3", 3, PxQ![-1890, 1890, 1647/2, -549], PxQ![-6, 6, -3/2, 1], PxQ![46656, -69984, 31104], PxQ!0, PxQ!0 >;
covers[9] := < "Q", "182a1", 3, PxQ![891, 2853, -891, 891], PxQ![1, -1, -1, 1], PxQ![-11232, 0, -11232, 22464], PxQ!0, PxQ!0 >;
covers[10] := < "Q", "33a2", 3, PxQ![192/11, 156/11, -606/11, -57], PxQ![-4/11, -4/11, 10/11, 1], PxQ![216/121, 1404/121, 216/11, 108/11], PxQ!0, PxQ!0 >;
K := NumberField(X^2 - 32246250*X + 2573441015625);
PK := PolynomialRing(K);
covers[11] := < "K", "", 5, PK![ K | 1/92625*(2*K.1 + 4803750), 1/7125*(K.1 + 264375), 1/18525*(4*K.1 + 499375), 1/85500*(17*K.1 + 6631875), 1/17784*(5*K.1 - 587625), 1/39000*(17*K.1 + 1831875) ], PK![ K | 0, 0, 1/13893750*(K.1 + 2401875), 1/6946875*(K.1 + 7033125), 1/13893750*(K.1 - 2229375), 1 ], PK![ K | 0, 1/55575*(-56*K.1 + 13695000), 1/2223*(-14*K.1 + 1200750), 1/55575*(-808*K.1 + 78510000), 1/55575*(-709*K.1 + 66208125), 1/18525*(23*K.1 - 331875), 1/18525*(601*K.1 - 47735625), 1/325*(11*K.1 - 879375) ], PK!0, PK!0 >;
coverK[11] := < K, K!(K.1), K!(1/52*(113515*K.1 - 9308159375)) >;
K := NumberField(X^2 + 138024*X - 1767690864);
PK := PolynomialRing(K);
covers[12] := < "K", "", 3, PK![ K | 1/10368*(-5*K.1 - 438372), 1/10368*(7*K.1 + 1042956), 1/5184*(13*K.1 - 1342332), 1/648*(-K.1 - 4860) ], PK![ K | 1/93312*(K.1 + 162324), 1/93312*(5*K.1 + 904932), 1/46656*(-K.1 - 69012), 1 ], PK![ K | 1/72*(-K.1 - 162324), 1/72*(-K.1 - 208980), 1/72*(K.1 - 70956), 1/36*(-K.1 + 24300) ], PK!0, PK!0 >;
coverK[12] := < K, K!(K.1), K!(93*K.1 - 1344276) >;
covers[13] := < "Q", "294b2", 3, PxQ![509/2, -148, -367/2, 435], PxQ![5/6, -4/3, -7/6, 1], PxQ![-216, -216, -324, 972], PxQ!0, PxQ!0 >;
covers[14] := < "Q", "462f1", 3, PxQ![4509/2, -654, 2505/2, -501], PxQ![-9/2, 10, -5/2, 1], PxQ![58806, -58806, 26136], PxQ!0, PxQ!0 >;
covers[15] := < "Q", "114a1", 3, PxQ![-1656, 156, -615/2, 147], PxQ![24, 4, -13/2, 1], PxQ![31104, -15552, -648, 648], PxQ!0, PxQ!0 >;
covers[16] := < "Q", "66c1", 3, PxQ![643/2, -4, -209/2, 3], PxQ![-5/6, -4/3, 7/6, 1], PxQ![-432, 432, 972, 324], PxQ!0, PxQ!0 >;
covers[17] := < "Q", "858k1", 3, PxQ![4954371, -5487933/4, 78465, -15693], PxQ![89, -55/4, -5, 1], PxQ![-161951724, 190531440, -80975862, 9526572], PxQ!0, PxQ!0 >;
K := NumberField(X^3 - 53859375/4*X^2 + 116033291015625/8*X - 138877345184326171875/64);
PK := PolynomialRing(K);
covers[18] := < "K", "", 2, PK![ K | 1/58446398437500*(-3376*K.1^2 + 32794375000*K.1 - 31187170107421875), 1/29223199218750*(1312*K.1^2 - 24895000000*K.1 + 17443671416015625), 1/58446398437500*(624*K.1^2 - 12088437500*K.1 - 5401994326171875) ], PK![ K | 1/43834798828125*(16*K.1^2 - 203947500*K.1 + 190810300781250), 1/73057998046875*(-32*K.1^2 + 440450000*K.1 - 352397402343750), 1 ], PK![ K | 1/14611599609375*(64*K.1^2 - 718125000*K.1 + 631736806640625), 1/14611599609375*(-16*K.1^2 + 301612500*K.1 - 322314697265625) ], PK!0, PK!0 >;
coverK[18] := < K, K!(K.1), K!(1/37405695000*(1548256*K.1^2 + 6254107537500*K.1 - 623782079912109375)) >;


t0 := Cputime();
nR := 0; nS := 0;
for row in rows do
    I, fc, hc, route, tag, ck := Explode(row);
    C := SimplifiedModel(HyperellipticCurve(R!fc, R!hc));
    if route eq "R1" then
        ok, s := SplitNode(RichelotIsogenousSurfaces(Jacobian(C)));
        assert ok;
        if Type(s) eq SetCart then
            printf "SPLIT %-14o (2,2)-isogenous /Q to E x E', conductors %o, %o  (%o)\n",
                I, Conductor(s[1]), Conductor(s[2]), tag;
        elif Type(BaseRing(s)) eq FldRat then
            printf "SPLIT %-14o (2,2)-isogenous /Q to E x E, E of conductor %o  (%o)\n",
                I, Conductor(s), tag;
        else
            printf "SPLIT %-14o (2,2)-isogenous to E x E^sigma over %o  (%o)\n",
                I, DefiningPolynomial(BaseRing(s)), tag;
        end if;
        nR +:= 1;
    elif route eq "Q2" then
        K := QuadraticField(2);
        JK := Jacobian(SimplifiedModel(HyperellipticCurve(
            PolynomialRing(K)!(4*(R!fc) + (R!hc)^2))));
        L := RichelotIsogenousSurfaces(JK);
        ok, s := SplitNode(L);
        if not ok then
            for u in L do
                if Type(u) eq JacHyp then
                    ok, s := SplitNode(RichelotIsogenousSurfaces(u));
                    if ok then break; end if;
                end if;
            end for;
        end if;
        assert ok;
        printf "SPLIT %-14o (2,2)-route to a product over Q(sqrt 2)  (%o)\n", I, tag;
        nR +:= 1;
    elif route eq "COVER" then
        cv := covers[ck];
        gg := 4*(R!fc) + (R!hc)^2;
        pv := cv[4]; qv := cv[5]; hv := cv[6];
        assert pv*Derivative(qv) - Derivative(pv)*qv ne 0;   // nonconstant
        if cv[1] eq "Q" then
            Ew := WeierstrassModel(EllipticCurve(CremonaDatabase(), cv[2]));
            WW := HyperellipticPolynomials(Ew);
            NN := pv^3 + Coefficient(WW,2)*pv^2*qv
                + Coefficient(WW,1)*pv*qv^2 + Coefficient(WW,0)*qv^3;
            assert NN*qv eq gg*hv^2;
            printf "SPLIT %-14o degree-%o map to elliptic curve %o over Q  (%o)\n",
                I, cv[3], cv[2], tag;
        else
            dat := coverK[ck];
            Kf := dat[1]; g2K := dat[2]; g3K := dat[3];
            assert g2K^3 - 27*g3K^2 ne 0;                    // E nonsingular
            NN := 4*pv^3 - g2K*pv*qv^2 - g3K*qv^3;
            PKl := Parent(pv);
            assert NN*qv eq (PKl!gg)*hv^2;
            printf "SPLIT %-14o degree-%o map to an elliptic curve over a degree-%o field  (%o)\n",
                I, cv[3], Degree(Kf), tag;
        end if;
        nR +:= 1;
    else // SIG
        g0 := 4*(R!fc) + (R!hc)^2;
        d0 := LCM([Z | Denominator(c) : c in Coefficients(g0)]);
        g := R!(d0^2*g0);
        D := Z!Discriminant(g);
        n := 0;
        for p in PrimesInInterval(3, 200) do
            if D mod p eq 0 or (Z!LeadingCoefficient(g)) mod p eq 0 then
                continue;
            end if;
            gp := PolynomialRing(GF(p))!g;
            if Degree(gp) lt 5 or not IsSquarefree(gp) then continue; end if;
            chi := Tp!Reverse(Coefficients(
                EulerFactor(Jacobian(HyperellipticCurve(gp)))));
            assert not IsStrict(chi);
            n +:= 1;
        end for;
        printf "SIG   %-14o odd gluing: no strict prime among %o good p<200; splitness cited  (%o)\n",
            I, n, tag;
        nS +:= 1;
    end if;
end for;
assert nR + nS eq #rows;
printf "ALL %o SPLITNESS CHECKS PASSED: %o by explicit isogeny or cover, %o signature+citation (%.1o s)\n",
    #rows, nR, nS, Cputime() - t0;
quit;
