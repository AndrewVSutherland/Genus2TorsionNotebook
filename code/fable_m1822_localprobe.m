//////////////////////////////////////////////////////////////////////
// fable_m1822_localprobe.m  (2026-07-18, Fable session)
//
// Finite-field viability probes on the odd M_1(8,2,2) family
// (f_{u,v} copied VERBATIM from notes/m3222_plus3.md), which carries
// [2,2,8] generically:
//
//   - contains [2,2,24]  (the +3 route, target of m3222_plus3)
//   - contains [2,2,16]  (the halving route, m3222_halving_boundary;
//                         known 7/11-adic boundary walls on the cover —
//                         here we ask the cruder chart-level question,
//                         including at q=13,17 beyond those walls)
//   - canary: contains [2,2,8]
//
// Run:  magma -b code/fable_m1822_localprobe.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(2*10^9);

function ContainsInv(iv, ta)
    k := #iv; m := #ta;
    if k lt m then return false; end if;
    return &and[ IsDivisibleBy(iv[k-m+j], ta[j]) : j in [1..m] ];
end function;

for q in [5,7,11,13,17] do
    K := GF(q);
    Pol<X> := PolynomialRing(K);
    tested := 0; good := 0; c228 := 0; c2224 := 0; c2216 := 0;
    ex24 := []; ex16 := [];
    for uv0, vv0 in K do
        u := uv0; v := vv0;
        f := ((1-u)*X+1)*((1-v)*X+1)*((u+v+2)*X+1)
             *(-X^2 + (u^2*v-u^2+u*v^2-u-v^2-v-2)*X
                     -(u^2+u*v+v^2+u+v+1));
        tested +:= 1;
        if Degree(f) ne 5 or not IsSeparable(f) then continue; end if;
        okc := true;
        try
            C := HyperellipticCurve(f);
            G := AbelianGroup(Jacobian(C));
        catch ee okc := false; end try;
        if not okc then continue; end if;
        good +:= 1;
        iv := Invariants(G);
        if ContainsInv(iv, [2,2,8])  then c228  +:= 1; end if;
        if ContainsInv(iv, [2,2,24]) then
            c2224 +:= 1;
            if #ex24 lt 4 then Append(~ex24, <u,v,iv>); end if;
        end if;
        if ContainsInv(iv, [2,2,16]) then
            c2216 +:= 1;
            if #ex16 lt 4 then Append(~ex16, <u,v,iv>); end if;
        end if;
    end for;
    printf "M1822 q=%o tested=%o good=%o canary228=%o contains2224=%o contains2216=%o\n",
        q, tested, good, c228, c2224, c2216;
    for e in ex24 do printf "M1822HIT24 q=%o u=%o v=%o inv=%o\n", q, e[1], e[2], e[3]; end for;
    for e in ex16 do printf "M1822HIT16 q=%o u=%o v=%o inv=%o\n", q, e[1], e[2], e[3]; end for;
end for;

printf "PROBE_DONE\n";
quit;
