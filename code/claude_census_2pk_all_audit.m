// Binary Theta-column audit: does EVERY known example of the group have a
// maximal-order doubled-point class 2P-K?
// Universe audited: all LMFDB geometrically-simple curves with the exact
// torsion group (data/claude_census_2pk_lmfdb_pool.m) + sample fibers of
// every repo-recorded infinite family for the bullet rows.
// Detector (as code/claude_census_2pk_audit.m): reduced Mumford <a,b,d>,
// doubled  <=>  d = 2 and (deg a = 2 with disc(a) = 0, or deg a <= 0).
SetColumns(0);
SetMemoryLimit(6*10^9);
QQ := Rationals();
R<x> := PolynomialRing(QQ);

procedure Audit(grp, src, C0)
    try
        C := IntegralModel(SimplifiedModel(C0));
        J := Jacobian(C);
        T, mp := TorsionSubgroup(J);
        inv := Invariants(T);
        if #inv eq 0 then printf "AUD|%o|%o|inv=[]|SKIP\n", grp, src; return; end if;
        nr := inv[#inv];
        found := false;
        for a in T do
            if IsIdentity(a) or Order(a) ne nr then continue; end if;
            g := mp(a);
            apol := g[1]; d := g[3];
            if d eq 2 and ((Degree(apol) eq 2 and Discriminant(apol) eq 0)
                           or Degree(apol) le 0) then
                found := true; break;
            end if;
        end for;
        printf "AUD|%o|%o|inv=%o|nr=%o|maxdoubled=%o\n", grp, src, inv, nr,
            found select "YES" else "NO";
    catch e
        printf "AUD|%o|%o|ERR\n", grp, src;
    end try;
end procedure;

// ---------- A) LMFDB pool ----------
load "data/claude_census_2pk_lmfdb_pool.m";
for rec in POOL do
    grp := rec[1]; lab := rec[2];
    fc := rec[3]; hc := rec[4];
    f := &+[ QQ!fc[i]*x^(i-1) : i in [1..#fc] ];
    h := #hc eq 0 select R!0 else &+[ QQ!hc[i]*x^(i-1) : i in [1..#hc] ];
    ok, C := IsHyperellipticCurve([f, h]);
    if not ok then printf "AUD|%o|%o|ERR|not_hyp\n", grp, lab; continue; end if;
    Audit(grp, lab, C);
end for;
printf "POOL_DONE\n";

// ---------- B) family fibers ----------
// contact-7 -> [7]
for ab in [[1,1],[2,1],[1,-1],[3,2]] do
    a := ab[1]; b := ab[2];
    h := 1 - (7/2)*x + a*x^2 + b*x^3;
    g := h^2 + (x-1)^7;
    f := g div x^2;
    if g mod x^2 eq 0 and IsSeparable(f) and Degree(f) ge 5 then
        Audit("7", Sprintf("contact7(%o,%o)", a, b), HyperellipticCurve(f));
    end if;
end for;
// contact-9 -> [9]
for a in [1, 2, -1, 3] do
    h := 1 - (9/2)*x + (63/8)*x^2 - (105/16)*x^3 + a*x^4;
    g := h^2 + (x-1)^9;
    f := g div x^4;
    if g mod x^4 eq 0 and IsSeparable(f) and Degree(f) ge 5 then
        Audit("9", Sprintf("contact9(%o)", a), HyperellipticCurve(f));
    end if;
end for;
// Flynn + Daowsud-Schmidt -> [11]
for t in [1, 2, 3, -2] do
    f := x^6 + 2*x^5 + (2*t + 3)*x^4 + 2*x^3 + (t^2 + 1)*x^2 + 2*t*(1 - t)*x + t^2;
    if IsSeparable(f) then Audit("11", Sprintf("Flynn11(%o)", t), HyperellipticCurve(f)); end if;
end for;
for u in [1, 2, 3, -1] do
    f := x^6 - 4*x^5 + 8*(1 + u)*x^4 - (10 + 32*u)*x^3 + 8*(1 + 6*u + 2*u^2)*x^2
         - 4*(1 + 6*u + 16*u^2)*x + 64*u^2 + 1;
    if IsSeparable(f) then Audit("11", Sprintf("DS11(%o)", u), HyperellipticCurve(f)); end if;
end for;
// M(12) z-chart -> [12]
for rz in [[3,2],[2,3],[5,2],[3,4]] do
    r := rz[1]; z := rz[2];
    a := (1-z^2)/(4*(r+1));
    T12 := a*x^2 - x + r;
    hh := (x-r)*(T12+1);
    ff := a*x^2*T12*(T12+1);
    if IsSeparable(hh^2+4*ff) then
        Audit("12", Sprintf("M12z(%o,%o)", r, z), HyperellipticCurve(ff, hh));
    end if;
end for;
// Leprevost 21 (deg-5 model, verbatim polynomials from z21_leprevost_family_verify.m)
K21<t> := FunctionField(QQ); R21<X> := PolynomialRing(K21);
p2 := t^14+4*t^13+19*t^12+32*t^11+113*t^10+188*t^9+379*t^8+448*t^7+379*t^6+188*t^5+113*t^4+32*t^3+19*t^2+4*t+1;
p1 := t^10+4*t^9+17*t^8+24*t^7+54*t^6+56*t^5+54*t^4+24*t^3+17*t^2+4*t+1;
p0 := t^6+4*t^5+15*t^4+16*t^3+15*t^2+4*t+1;
A21 := p2*X^2-2*(t^2+1)^2*p1*X+(t^2+1)^4*p0;
k21 := 64*t^4*(t+1)^2*(t^2+1)^3*(t^4+2*t^3+6*t^2+2*t+1)^3;
f21 := A21^2-k21*X^3*(X-1)^2;
for tv in [1, 2, 3, -3] do
    fv := &+[ QQ!Evaluate(Coefficient(f21, i), tv)*x^i : i in [0..5] ];
    if IsSeparable(fv) and Degree(fv) eq 5 then
        Audit("21", Sprintf("Lep21(%o)", tv), HyperellipticCurve(fv));
    end if;
end for;
// Kuru-Sadek 23
for tv in [2, 3, 4, -2] do
    tq := QQ!tv;
    beta := (tq^2 + 1)^2/(4*tq^2);
    sbeta := (tq^2 + 1)/(2*tq);
    s := (tq^2 - 1)/(2*tq);
    alpha := beta - s^5/(beta*sbeta);
    lambda := (alpha - 1)^4/((alpha - beta)^2*alpha);
    num := x^3*(x - alpha)^2 - (x - 1)*((x - 1)^4 - lambda*(x - beta)^2*x);
    den := 2*(x - alpha)*(x - beta);
    A, rem := Quotrem(num, den);
    if rem ne 0 then continue; end if;
    f := A^2 - lambda*x^4*(x - 1);
    if IsSeparable(f) and Degree(f) ge 5 then
        Audit("23", Sprintf("KS23(%o)", tv), HyperellipticCurve(f));
    end if;
end for;
printf "FIBERS_DONE\n";
quit;
