
//////////////////////////////////////////////////////////////////////
//  Deep Mordell-Weil enumeration on an ELS fiber of the P_R-halving
//  problem in M_1(8,4).  Proper multi-generator version of the MW part
//  of agent_m18_416_els_fiber_attack.m.
//
//  Usage: magma -b Rnum:=-25 Rden:=4 mex_num:=-25 mex_den:=2 N:=5 \
//             agent_m18_416_els_mw_deep.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
Z := Integers();

if not assigned Rnum then Rnum := -25; end if;
if not assigned Rden then Rden := 4; end if;
if not assigned mex_num then mex_num := -25; end if;
if not assigned mex_den then mex_den := 2; end if;
if not assigned N then N := 5; end if;
for nm in ["Rnum","Rden","mex_num","mex_den","N"] do ; end for;
if Type(Rnum) eq MonStgElt then Rnum := StringToInteger(Rnum); end if;
if Type(Rden) eq MonStgElt then Rden := StringToInteger(Rden); end if;
if Type(mex_num) eq MonStgElt then mex_num := StringToInteger(mex_num); end if;
if Type(mex_den) eq MonStgElt then mex_den := StringToInteger(mex_den); end if;
if Type(N) eq MonStgElt then N := StringToInteger(N); end if;

Rv := Q!Rnum/Rden;
mex := Q!mex_num/mex_den;
Kv := -2*Rv*(Rv^2-1);
HeightCap := 10^60;   // skip m with numerator/denominator beyond this

P<xq> := PolynomialRing(Q);

function IsSquareQ(qv)
    qv := Q!qv;
    if qv le 0 then return false, 0; end if;
    okn, sn := IsSquare(Numerator(qv));
    okd, sd := IsSquare(Denominator(qv));
    if okn and okd then return true, sn/sd; end if;
    return false, 0;
end function;

function FamilyData(Rv, wv)
    tv := (2*Rv^2 + (1-wv^2)*Rv - 2*wv^2)/(4*(wv^2-1));
    c4v := Rv + 2 + 4*tv;
    Av := xq^2 + (Rv^3 + 4*Rv^2*tv + Rv - 8*Rv*tv + 4*tv)*xq + Rv^4;
    Bv := c4v*xq^2 + (Rv^2 + 4*Rv + 1 + 8*tv)*xq + (2*Rv^2 + Rv + 4*tv);
    return xq*Av*Bv, Av, Bv, c4v;
end function;

function SecondStagePass(Rv, wv)
    fq, Av, Bv, c4v := FamilyData(Rv, wv);
    if Degree(fq) ne 5 or Discriminant(fq) eq 0 then return false; end if;
    XR := -c4v*Rv;
    Atq := P![Q!co : co in Coefficients(c4v^2*Evaluate(Av, xq/c4v))];
    Btq := P![Q!co : co in Coefficients(c4v*Evaluate(Bv, xq/c4v))];
    for gq in [Atq, Btq] do
        dsc := Discriminant(gq);
        if dsc eq 0 then return false; end if;
        okd, _ := IsSquareQ(dsc);
        if okd then
            for rt in Roots(gq) do
                val := XR - rt[1];
                if val eq 0 then return false; end if;
                okv, _ := IsSquareQ(val);
                if not okv then return false; end if;
            end for;
        else
            Kf<th> := NumberField(gq);
            if not IsSquare(XR - th) then return false; end if;
        end if;
    end for;
    return true;
end function;

printf "DEEP MW: fiber R=%o, N=%o\n", Rv, N;
Pm<mm> := PolynomialRing(Q);
Gam := 2*(Rv^2-1)*(Rv*(2*Rv+1)*(mm^2-Kv)^2 - (Rv+2)*(mm^2+Kv)^2);
C := HyperellipticCurve(Gam);
okg, g0 := IsSquareQ(Evaluate(Gam, mex));
error if not okg, "bad example point";
pt0 := C![mex, g0];
E, phi := EllipticCurve(C, pt0);
Emin, mpmin := MinimalModel(E);
printf "E = %o\n", aInvariants(Emin);
MW, mMW := MordellWeilGroup(Emin);
invs := Invariants(MW);
printf "MW group invariants: %o\n", invs;
freeidx := [i : i in [1..#invs] | invs[i] eq 0];
torsidx := [i : i in [1..#invs] | invs[i] ne 0];
r := #freeidx;
printf "rank %o, torsion part %o\n", r, [invs[i] : i in torsidx];

phiInv := Inverse(phi);
mpminInv := Inverse(mpmin);

// enumerate torsion combos
torscombos := [[]];
for i in torsidx do
    newc := [];
    for c in torscombos do
        for v in [0..invs[i]-1] do
            Append(~newc, c cat [v]);
        end for;
    end for;
    torscombos := newc;
end for;

// enumerate free combos with sup-norm <= N
count := 0; tested := 0; passes := 0; skippedHeight := 0;
boxes := [ [n : n in [-N..N]] : i in [1..r] ];
idx := [ -N : i in [1..r] ];
done := false;
seenm := {};
while not done do
    // build free part vector
    for tc in torscombos do
        coeffs := [0 : i in [1..#invs]];
        for j in [1..r] do coeffs[freeidx[j]] := idx[j]; end for;
        for j in [1..#torsidx] do coeffs[torsidx[j]] := tc[j]; end for;
        g := MW!coeffs;
        Ept := mMW(g);
        if Ept eq Emin!0 then continue; end if;
        count +:= 1;
        try
            Cpt := phiInv(mpminInv(Ept));
        catch e continue; end try;
        if Cpt[3] eq 0 then continue; end if;
        mv := Cpt[1]/Cpt[3];
        if mv in seenm then continue; end if;
        Include(~seenm, mv);
        if Abs(Numerator(mv)) gt HeightCap or Abs(Denominator(mv)) gt HeightCap then
            skippedHeight +:= 1; continue;
        end if;
        den0 := mv^2 - Kv;
        if den0 eq 0 then continue; end if;
        wv := (mv^2 + Kv)/den0;
        if wv in {Q!0, Q!1, Q!-1} then continue; end if;
        tested +:= 1;
        if SecondStagePass(Rv, wv) then
            passes +:= 1;
            printf "!!!! FULL DESCENT PASS R=%o m=%o\n", Rv, mv;
            // certify
            fq, Av, Bv, c4v := FamilyData(Rv, wv);
            L := 1;
            for i in [0..Degree(fq)] do L := LCM(L, Denominator(Coefficient(fq, i))); end for;
            fI := P!(L^2*fq);
            J := Jacobian(HyperellipticCurve(fI));
            Qf := Rv^2 - (Q!1/2)*Rv*wv^2 + (Q!1/2)*Rv - wv^2;
            YRv := -2*Rv*(Rv-1)^2*Qf/(wv^2-1);
            PR := J![xq + Rv, P!(L*YRv)];
            ok, half := IsDivisibleBy(PR, 2);
            printf "IsDivisibleBy = %o", ok;
            if ok then
                printf ", half order %o", Order(half);
                printf ", torsion %o", Invariants(TorsionSubgroup(J));
            end if;
            printf "\nf = %o\n", fI;
        end if;
    end for;
    // increment idx
    j := 1;
    while j le r do
        idx[j] +:= 1;
        if idx[j] le N then break; end if;
        idx[j] := -N; j +:= 1;
    end while;
    if j gt r then done := true; end if;
    if r eq 0 then done := true; end if;
end while;
printf "MW points visited=%o distinct m tested=%o skipped(height)=%o passes=%o\n",
    count, tested, skippedHeight, passes;
print "DONE";
quit;
