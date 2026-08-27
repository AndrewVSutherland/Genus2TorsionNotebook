SetColumns(0);
Q := Rationals();
P<z> := PolynomialRing(Q);
Rv := Q!-8;
Kv := -2*Rv*(Rv^2-1);
Gam := 2*(Rv^2-1)*(Rv*(2*Rv+1)*(z^2-Kv)^2 - (Rv+2)*(z^2+Kv)^2);
C := HyperellipticCurve(Gam);
m0 := Q!-28;
okg0, g0 := IsSquare(Evaluate(Gam, m0)); assert okg0;
P0 := C![m0, g0, Q!1];
Eraw, phi := EllipticCurve(C, P0);
Emin, mp := MinimalModel(Eraw);
MW, toE := MordellWeilGroup(Emin);
invs := Invariants(MW);
freeidx := [i : i in [1..#invs] | invs[i] eq 0];
Gmw := (MW.freeidx[1]);
Gpt := toE(Gmw);
Qp := pAdicField(11, 40);
Eminp := BaseChange(Emin, Qp);
Gp := Eminp![Qp!(Gpt[1]/Gpt[3]), Qp!(Gpt[2]/Gpt[3]), Qp!1];
printf "Gpt=%o\n", Gpt;
printf "Gp=%o\n", Gp;
mpInv := Inverse(mp);
phiInv := Inverse(phi);
printf "mpInv equations=%o\n", DefiningEquations(mpInv);
printf "phiInv equations=%o\n", DefiningEquations(phiInv);
try
    psi := mpInv*phiInv;
    printf "composition equations=%o\n", DefiningEquations(psi);
catch e
    printf "composition construction failed: %o\n", e`Object;
end try;
try
    Cpt := phiInv(mpInv(Gpt));
    printf "Q rational pullback=%o\n", Cpt;
catch e
    printf "Q rational pullback failed: %o\n", e`Object;
end try;
try
    Cptp := phiInv(mpInv(Gp));
    printf "Qp pullback through Q maps=%o\n", Cptp;
catch e
    printf "Qp pullback through Q maps failed: %o\n", e`Object;
end try;
try
    mpInvp := BaseChange(mpInv, Qp);
    phiInvp := BaseChange(phiInv, Qp);
    Cptp2 := phiInvp(mpInvp(Gp));
    printf "Qp pullback through basechanged maps=%o\n", Cptp2;
catch e
    printf "Qp pullback through basechanged maps failed: %o\n", e`Object;
end try;
quit;
