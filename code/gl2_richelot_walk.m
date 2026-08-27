// 37-hunt: EXACT (2,2)-isogeny walk from the two exact cubic-field members.
// Every vertex of the Richelot graph is a PP member of the isogeny class;
// a Q-rational member shows up as a vertex whose G2 invariants have
// degree-1 minimal polynomials.  Pure exact arithmetic, no precision.
//
// magma -b code/gl2_richelot_walk.m > results/gl2_richelot_walk.log

SetColumns(0);
SetSeed(1);
SetMemoryLimit(24*10^9);
Q := Rationals();
Px<X> := PolynomialRing(Q);

// curve 1 over K3a = Q[t]/(t^3+3t-54)  (from gl2_exact_descend.m)
K3a<t> := NumberField(X^3 + 3*X - 54);
Pa<xa> := PolynomialRing(K3a);
f1 := 1/70761*(304360*t^2 + 194256*t - 5052516)*xa^6
    + 1/23587*(-413320*t^2 - 1162368*t + 4982718)*xa^5
    + 1/94348*(315736*t^2 - 397440*t - 17223893)*xa^4
    + 1/141522*(-1137241*t^2 + 801168*t - 5469876)*xa^3
    + 1/94348*(786865*t^2 + 2043498*t + 17488258)*xa^2
    + 1/23587*(-50136*t^2 - 165408*t - 1491258)*xa
    + 1/94348*(-12331*t^2 + 184258*t + 411795);

// curve 2 over K3b = Q[t]/(t^3 - t^2 - 51t - 129)  (v17 log, member n1)
K3b<u> := NumberField(X^3 - X^2 - 51*X - 129);
Pb<xb> := PolynomialRing(K3b);
f2 := 1/188696*(-79338395*u^2 + 458748990*u + 3265747859)*xb^6
    + 1/47174*(75123387*u^2 - 422011134*u - 2913676335)*xb^5
    + 1/188696*(-461732739*u^2 + 2531262798*u + 16983850529)*xb^4
    + 1/94348*(185430019*u^2 - 995756958*u - 6499444593)*xb^3
    + 1/94348*(-82561695*u^2 + 436329990*u + 2768085958)*xb^2
    + 1/23587*(4866660*u^2 - 25476120*u - 156581085)*xb
    + 1/94348*(-1921905*u^2 + 10041210*u + 59563845);

procedure Walk(f, tagbase, depth)
    K := BaseRing(Parent(f));
    C0 := HyperellipticCurve(f);
    seen := {};
    frontier := [ C0 ];
    g20 := G2Invariants(C0);
    Include(~seen, g20);
    for lev in [1..depth] do
        newf := [];
        for C in frontier do
            J := Jacobian(C);
            surfs := [];
            try
                surfs := RichelotIsogenousSurfaces(J);
            catch e
                printf "%o L%o RICHELOT_ERR %o\n", tagbase, lev,
                    Substring(Sprint(e`Object),1,80);
            end try;
            for s in surfs do
                // s is either Jac(curve) or a product of elliptic curves
                if Type(s) eq SetCart then continue; end if;   // E1 x E2 product
                Cs := 0;
                try Cs := Curve(s); catch e ; end try;
                if Type(Cs) eq RngIntElt then continue; end if;
                g2 := G2Invariants(Cs);
                if g2 in seen then continue; end if;
                Include(~seen, g2);
                degs := [ Degree(MinimalPolynomial(g2[k], Q)) : k in [1..3] ];
                printf "%o L%o VERTEX degs=%o\n", tagbase, lev, degs;
                if Max(degs) eq 1 then
                    g2Q := [ Q!(-Coefficient(MinimalPolynomial(g2[k],Q),0)
                              / Coefficient(MinimalPolynomial(g2[k],Q),1)) : k in [1..3] ];
                    printf "%o L%o *** Q_VERTEX G2=%o ***\n", tagbase, lev, g2Q;
                end if;
                Append(~newf, Cs);
            end for;
        end for;
        printf "%o L%o frontier %o total %o\n", tagbase, lev, #newf, #seen;
        frontier := newf;
        if #frontier eq 0 then break; end if;
    end for;
end procedure;

Walk(f1, "K3a", 3);
Walk(f2, "K3b", 3);
printf "RICHELOT_WALK_DONE\n";
quit;
