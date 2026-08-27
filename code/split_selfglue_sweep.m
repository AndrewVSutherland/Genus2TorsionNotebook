// (A) Is our 11a1x11a3 glued [5,5] curve = X_0(22)?
// (B) Cor-7 SELF-gluing sweep: E_t with N-torsion (N=5,7,9) and SQUARE
//     disc of the 2-division cubic (A3-cubic): rotation psi gives
//     C with J ~ (2,2)-isog E x E, torsion >= [N,N].  Target: [9,9].
SetColumns(0);
SetSeed(1);
SetMemoryLimit(4*10^9);
Q := Rationals();
P<x> := PolynomialRing(Q);

// (A)
C1 := HyperellipticCurve(P![-2,4,2,5,2,1], P![1,1,1,1]);  // y^2+(x^3+x^2+x+1)y = x^5+2x^4+5x^3+2x^2+4x-2
X022 := HyperellipticCurve(-2*x^6 - 10*x^4 + 26*x^2 + 242);
printf "OUR55_EQ_X022 %o\n", G2Invariants(SimplifiedModel(C1)) eq G2Invariants(X022);
printf "X022_TORSION %o\n", Invariants(TorsionSubgroup(Jacobian(X022)));

// (B)
function TateE(b, c) return EllipticCurve([1-c, -b, -b, 0, 0]); end function;
function KubertBC(N, d)
    if N eq 5 then return d, d; end if;
    if N eq 7 then return d^3-d^2, d^2-d; end if;
    c := d^2*(d-1); return c*(d*(d-1)+1), c;
end function;
function MonicCubic(E)
    W := WeierstrassModel(E); a := aInvariants(W);
    return x^3 + a[4]*x + a[5];
end function;
function GlueP4(f, g, K)
    rf := [ r[1] : r in Roots(PolynomialRing(K)!f) ];
    rg0 := [ r[1] : r in Roots(PolynomialRing(K)!g) ];
    Df := Discriminant(f); Dg := Discriminant(g);
    out := [];
    for s in SymmetricGroup(3) do
        rg := [ rg0[1^s], rg0[2^s], rg0[3^s] ];
        al := rf; be := rg;
        d_al := [ al[3]-al[2], al[1]-al[3], al[2]-al[1] ];
        d_be := [ be[3]-be[2], be[1]-be[3], be[2]-be[1] ];
        if 0 in d_al or 0 in d_be then continue; end if;
        a1 := d_al[1]^2/d_be[1] + d_al[3]^2/d_be[3] + d_al[2]^2/d_be[2];
        b1 := d_be[1]^2/d_al[1] + d_be[3]^2/d_al[3] + d_be[2]^2/d_al[2];
        a2 := al[1]*d_be[1] + al[2]*d_be[2] + al[3]*d_be[3];
        b2 := be[1]*d_al[1] + be[2]*d_al[2] + be[3]*d_al[3];
        if a1 eq 0 or a2 eq 0 or b1 eq 0 or b2 eq 0 then continue; end if;
        A := Dg*a1/a2; B := Df*b1/b2;
        PK<xx> := PolynomialRing(K);
        h := -( A*d_al[3]*d_al[2]*xx^2 + B*d_be[3]*d_be[2] )
            *( A*d_al[1]*d_al[3]*xx^2 + B*d_be[1]*d_be[3] )
            *( A*d_al[2]*d_al[1]*xx^2 + B*d_be[2]*d_be[1] );
        cf := Coefficients(h);
        ok := true; cq := [];
        for co in cf do
            fl, cr := IsCoercible(Q, co);
            if not fl then ok := false; break; end if;
            Append(~cq, cr);
        end for;
        if not ok then continue; end if;
        hQ := P!cq;
        if Degree(hQ) eq 6 and IsSquarefree(hQ) then Append(~out, hQ); end if;
    end for;
    return out;
end function;

hts := [];
for den in [1..16] do
    for num in [-24..24] do
        if num eq 0 or GCD(Abs(num), den) ne 1 then continue; end if;
        Append(~hts, Q!num/den);
    end for;
end for;
seen := {};
for N in [9, 7, 5] do
    nsq := 0;
    for tv in hts do
        bN, cN := KubertBC(N, tv);
        ok := true;
        try
            E := TateE(bN, cN);
        catch e ok := false; end try;
        if not ok then continue; end if;
        if Order(E![0,0]) ne N then continue; end if;
        f := MonicCubic(E);
        if not IsIrreducible(f) then continue; end if;
        if not IsSquare(Discriminant(f)) then continue; end if;  // A3 cubic
        nsq +:= 1;
        K := SplittingField(f);   // cyclic cubic
        for hQ in GlueP4(f, f, K) do
            C := HyperellipticCurve(hQ);
            ii := G2Invariants(C);
            if ii in seen then continue; end if;
            Include(~seen, ii);
            Cm := C;
            try Cm := ReducedMinimalWeierstrassModel(C); catch e ; end try;
            T := TorsionSubgroup(Jacobian(SimplifiedModel(Cm)));
            fc, hc := HyperellipticPolynomials(Cm);
            printf "SELFGLUE N=%o t=%o TORSION %o : y^2+(%o)y = %o\n", N, tv, Invariants(T), hc, fc;
        end for;
    end for;
    printf "N=%o square-disc curves found: %o\n", N, nsq;
end for;
printf "SELFGLUE_DONE\n";
quit;
