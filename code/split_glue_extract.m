// Rebuild and print minimal models for the [5,5] and [35] gluing hits
SetColumns(0);
Q := Rationals();
P<x> := PolynomialRing(Q);
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

pairs := [ <5,5,Q!-7,Q!7/5>, <5,5,Q!-2,Q!8>, <5,5,Q!-1,Q!11>, <5,5,Q!10,Q!11/2>,
           <5,7,Q!10,Q!-10>, <5,7,Q!14,Q!8>, <5,7,Q!11/2,Q!-10>, <5,7,Q!13/2,Q!-1>, <5,7,Q!5/4,Q!-4> ];
for pr in pairs do
    M := pr[1]; N := pr[2]; tv := pr[3]; uv := pr[4];
    bM, cM := KubertBC(M, tv); E := TateE(bM, cM);
    bN, cN := KubertBC(N, uv); F := TateE(bN, cN);
    f := MonicCubic(E); g := MonicCubic(F);
    K := SplittingField(f);
    for hQ in GlueP4(f, g, K) do
        C := HyperellipticCurve(hQ);
        Cm := C;
        try Cm := ReducedMinimalWeierstrassModel(C); catch e ; end try;
        T := TorsionSubgroup(Jacobian(SimplifiedModel(Cm)));
        fc, hc := HyperellipticPolynomials(Cm);
        printf "%ox%o (t,u)=(%o,%o) torsion %o : y^2 + (%o)y = %o  [E=%o F=%o]\n",
            M, N, tv, uv, Invariants(T), hc, fc,
            CremonaReference(E), CremonaReference(F);
    end for;
end for;
printf "EXTRACT_DONE\n";
quit;
