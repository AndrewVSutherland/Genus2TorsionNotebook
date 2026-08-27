// Minimal models for the new-group and upgraded-row witnesses
SetColumns(0);
Q := Rationals();
P<x> := PolynomialRing(Q);
pairs := [
  <"[2,2,12]", "15.a6", "90.c3">, <"[2,2,2,8]", "21.a5", "210.e6">,
  <"[2,2,4,4]", "210.c5", "2310.o4">, <"[40]", "48.a6", "66.c4">,
  <"[60]", "150.c3", "90.c7">, <"[5,10]", "150.c3", "570.l4">,
  <"[2,24]", "48.a3", "30.a6">, <"[6,12]", "8190.bx1", "155610.fa1"> ];
lines := Split(Read("data/split_ec_torsion_pool.txt"), "\n");
ai := AssociativeArray();
for L in lines do
    if #L eq 0 then continue; end if;
    parts := Split(L, "|");
    ai[parts[1]] := eval parts[2];
end for;
function MonicCubic(lab, ai)
    E := EllipticCurve([Q!a : a in ai[lab]]);
    W := WeierstrassModel(E);
    av := aInvariants(W);
    return x^3 + av[4]*x + av[5];
end function;
function GlueP4All(f, g)
    // splitting field: compositum handling both reducible cases
    K := SplittingField(f*g);
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
for pr in pairs do
    tgt := pr[1];
    f := MonicCubic(pr[2], ai); g := MonicCubic(pr[3], ai);
    for hQ in GlueP4All(f, g) do
        C := HyperellipticCurve(hQ);
        Cm := C;
        try Cm := ReducedMinimalWeierstrassModel(C); catch e ; end try;
        T := TorsionSubgroup(Jacobian(SimplifiedModel(Cm)));
        istr := Sprintf("[%o]", &cat[ IntegerToString(v) cat "," : v in Invariants(T)]);
        if Sprint(Invariants(T)) eq Sprint(eval Sprintf("%o", tgt)) or true then
            fc, hc := HyperellipticPolynomials(Cm);
            printf "WITNESS %o glue %o x %o : torsion %o : y^2+(%o)y = %o\n",
                tgt, pr[2], pr[3], Invariants(T), hc, fc;
        end if;
    end for;
end for;
printf "EXTRACT2_DONE\n";
quit;
