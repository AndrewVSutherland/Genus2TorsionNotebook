SetColumns(0);
SetMemoryLimit(4*10^9);
QQ := Rationals();
R<x> := PolynomialRing(QQ);

function StrictPrime(C : lo := 5, hi := 200)
    for p in PrimesInInterval(lo, hi) do
        ok := true;
        try
            Cp := ChangeRing(C, GF(p));
            L := LPolynomial(Cp);
        catch e
            ok := false;
        end try;
        if not ok then continue; end if;
        if Degree(L) ne 4 then continue; end if;
        chi := Polynomial(Reverse(Coefficients(L)));
        if not IsIrreducible(chi) then continue; end if;
        K := NumberField(chi);
        pi := K.1;
        good := true;
        for n in [1..12] do
            if Degree(MinimalPolynomial(pi^n)) ne 4 then good := false; break; end if;
        end for;
        if good then return p; end if;
    end for;
    return 0;
end function;

procedure TestFiber(name, C0)
    try
        C := IntegralModel(SimplifiedModel(C0));
        J := Jacobian(C);
        T, mp := TorsionSubgroup(J);
        sp := StrictPrime(C);
        printf "FIBER|%o|inv=%o|strictprime=%o\n", name, Invariants(T), sp;
    catch e
        printf "FIBER|%o|ERR\n", name;
    end try;
end procedure;

TestFiber("control66", HyperellipticCurve(1872*x^5-3000*x^4+6969*x^3-1691*x^2+4875*x));

// contact-6 fibers (exact [2,6] earlier)
for ab in [[1,1],[2,-1],[3,1]] do
    a := ab[1]; b := ab[2];
    h6 := 1 + a*x + b*x^2 + x^3;
    f := h6^2 - (x-1)^6;
    if Degree(f) ge 5 and IsSeparable(f) then
        TestFiber(Sprintf("contact6(a=%o,b=%o)", a, b), HyperellipticCurve(f));
    end if;
end for;

// contact-7 fibers (exact [7] earlier)
for ab in [[1,1],[2,1],[1,-1],[3,2]] do
    a := ab[1]; b := ab[2];
    h := 1 - (7/2)*x + a*x^2 + b*x^3;
    g := h^2 + (x-1)^7;
    f := g div x^2;
    if g mod x^2 eq 0 and IsSeparable(f) and Degree(f) ge 5 then
        TestFiber(Sprintf("contact7(a=%o,b=%o)", a, b), HyperellipticCurve(f));
    end if;
end for;

// plain M(12) chart fibers (exact [6] earlier)
for ar in [[1,3],[2,3],[1,-2],[3,2]] do
    a := ar[1]; r := ar[2];
    T12 := a*x^2 - x + r;
    hh := (x-r)*(T12+1);
    ff := a*x^2*T12*(T12+1);
    if IsSeparable(hh^2+4*ff) then
        TestFiber(Sprintf("M12plain(a=%o,r=%o)", a, r), HyperellipticCurve(ff, hh));
    end if;
end for;

// M(12) rational-Weierstrass z-chart: a = (1-z^2)/(4*(r+1)) -> order 12 marked
for rz in [[3,2],[2,3],[5,2],[3,4]] do
    r := rz[1]; z := rz[2];
    a := (1-z^2)/(4*(r+1));
    T12 := a*x^2 - x + r;
    hh := (x-r)*(T12+1);
    ff := a*x^2*T12*(T12+1);
    if IsSeparable(hh^2+4*ff) then
        TestFiber(Sprintf("M12z(r=%o,z=%o)", r, z), HyperellipticCurve(ff, hh));
    end if;
end for;

// [2,12] line fibers (exact [2,12] earlier)
for rv in [3, 5, -3] do
    a := (1-rv)/4;
    T12 := a*x^2 - x + rv;
    hh := (x-rv)*(T12+1);
    ff := a*x^2*T12*(T12+1);
    if IsSeparable(hh^2+4*ff) then
        TestFiber(Sprintf("Z2xZ12(r=%o)", rv), HyperellipticCurve(ff, hh));
    end if;
end for;

print "FIBERS_DONE";
quit;
