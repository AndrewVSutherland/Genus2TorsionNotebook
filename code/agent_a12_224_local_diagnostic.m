
//////////////////////////////////////////////////////////////////////
//  Z/2 x Z/24 local (mod-p) diagnostic on the A(12) chart.
//
//  Is the [2,24] target locally OBSTRUCTED (like [4,16] over F_7) or
//  merely HEIGHT-LIMITED?  For each prime q, count chart points (p,z,r)
//  over F_q where:
//    (H)  P12 (order 12) or a 2-torsion translate is divisible by 2
//         (=> an order-24 point exists over F_q),
//    (R)  the sextic f = R*F has 2-rank >= 2 (factor type supporting an
//         independent rational 2-torsion class),
//    (HR) both -- the necessary local condition for a rational [2,24].
//  Also reports, among (HR) points, how many have the finite torsion
//  group actually containing Z/2 x Z/24.
//
//  Nonempty (HR) at all small q  =>  no local obstruction; push height.
//  Empty (HR) at some q          =>  boundary/obstruction analysis, as
//                                    for [4,16].
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB := 2; elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

// 2-rank of a squarefree deg-6 poly over F_q from its factor type
function TwoRank(fp)
    degs := [Degree(g[1]) : g in Factorization(fp)];
    // add a root at infinity if odd degree (deg 5 model) -- here deg 6
    // dim = (#even-Galois-stable-subsets)/(∅,all): use orbit parity formula
    k := #degs;
    // count subsets of orbits with even total size
    even := 0;
    for mask in [0..2^k - 1] do
        s := 0;
        for i in [1..k] do
            if (mask div 2^(i-1)) mod 2 eq 1 then s +:= degs[i]; end if;
        end for;
        if s mod 2 eq 0 then even +:= 1; end if;
    end for;
    // even subsets form F_2-space of dim log2(even); mod (∅ ~ all) => -1
    return Ilog2(even) - 1;
end function;

for q in [11, 13, 17, 19, 23] do
    Fq := GF(q);
    PF<xf> := PolynomialRing(Fq);
    nchart := 0; nH := 0; nR := 0; nHR := 0; ntors := 0;
    for pv in Fq do
        if pv eq 0 then continue; end if;
        for zv in Fq do
            if zv eq 0 then continue; end if;
            sv := (zv^2 - 4*pv^2 + 1)/(2*zv);
            if sv^2 eq 1 then continue; end if;
            tv := (zv^2 + 4*pv^2 - 1)^2/(8*pv^2*zv);
            for rv in Fq do
                muv := ((sv^2 - 1)*(2*pv*rv + 1) - pv^2*(2*sv*tv - 4))/(4*pv^3);
                den := sv^2 - 1;
                lav := (4 - muv^2)*pv^2/den;
                if lav eq 0 then continue; end if;
                T1v := pv*xf + rv;
                Rv := (T1v^2 + xf - 1)/lav;
                ellv := sv*xf + tv;
                Qv := 2*T1v + muv*Rv;
                Fpol := Rv*xf^2 + 4*(Rv + xf - 1)*(Rv - 1);
                fv := Rv*Fpol;
                if Degree(fv) ne 6 or not IsSquarefree(fv) then continue; end if;
                nchart +:= 1;
                rank2 := TwoRank(fv) ge 2;
                if rank2 then nR +:= 1; end if;
                // order-12 point + halving
                u4 := Qv/LeadingCoefficient(Qv);
                u6 := (Rv + xf - 1)/LeadingCoefficient(Rv + xf - 1);
                J := Jacobian(HyperellipticCurve(fv));
                O := J!0;
                halvable := false;
                try
                    P4 := J![u4, (Rv*ellv) mod u4];
                    P6 := J![u6, (xf*Rv) mod u6];
                    P12 := P4 + P6;
                    if 12*P12 eq O and &and[n*P12 ne O : n in [1..11]] then
                        G, phi := AbelianGroup(J);
                        // divisible-by-2 test for P12 and its 2-torsion translates
                        Tors2 := [g : g in G | 2*g eq G!0];
                        aP := P12 @@ phi;
                        invs := Invariants(G);
                        for tg in Tors2 do
                            cc := Eltseq(aP + tg);
                            ok := true;
                            for i in [1..#cc] do
                                if cc[i] mod GCD(2, invs[i]) ne 0 then ok := false; break; end if;
                            end for;
                            if ok then halvable := true; break; end if;
                        end for;
                    end if;
                catch e ;
                end try;
                if halvable then nH +:= 1; end if;
                if halvable and rank2 then
                    nHR +:= 1;
                    // does J(F_q) contain Z/2 x Z/24?
                    G2 := AbelianGroup(J);
                    inv := Invariants(G2);
                    has := #[n : n in inv | n mod 24 eq 0] ge 1
                           and #[n : n in inv | n mod 2 eq 0] ge 2;
                    if has then ntors +:= 1; end if;
                end if;
            end for;
        end for;
    end for;
    printf "q=%o: chart_pts=%o  halvable(H)=%o  rank>=2(R)=%o  BOTH(HR)=%o  J(F_q) has [2,24]-type=%o\n",
        q, nchart, nH, nR, nHR, ntors;
end for;
print "DONE";
quit;
