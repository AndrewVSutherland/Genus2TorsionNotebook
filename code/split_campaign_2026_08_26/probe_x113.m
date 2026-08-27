SetColumns(0);
SetMemoryLimit(12*10^9);
Q := Rationals();
Qrs<r,s> := FunctionField(Q, 2);
cS := s*(r-1); bS := r*cS;
Egen := EllipticCurve([1-cS, -bS, -bS, 0, 0]);
psiN := DivisionPolynomial(Egen, 13);
F := Numerator(Evaluate(psiN, 0));
P2<R,S> := PolynomialRing(Q, 2);
hh := hom< Parent(F) -> P2 | [R, S] >;
FP := hh(F);
for fc in Factorization(FP) do
    printf "factor degR=%o degS=%o mult=%o small=%o\n", Degree(fc[1],R), Degree(fc[1],S), fc[2],
        TotalDegree(fc[1]) le 8 select fc[1] else P2!0;
end for;
// the X1(13) factor: try the largest; check genus + hyperelliptic
big := [ fc[1] : fc in Factorization(FP) | Degree(fc[1],R) ge 2 and Degree(fc[1],S) ge 2 ];
for FB in big do
    A2<uu,vv> := AffineSpace(Q,2);
    hcr := hom< P2 -> CoordinateRing(A2) | [uu,vv] >;
    C := ProjectiveClosure(Curve(A2, hcr(FB)));
    printf "curve from degR=%o degS=%o: genus %o\n", Degree(FB,R), Degree(FB,S), Genus(C);
    if Genus(C) eq 2 then
        okh, H, mp := IsHyperelliptic(C);
        printf "hyperelliptic: %o\n", okh;
        if okh then
            printf "model: %o\n", HyperellipticPolynomials(H);
        end if;
    end if;
end for;
quit;
