//////////////////////////////////////////////////////////////////////
// Probe the elliptic-surface model for A/M(2,4,4).
//
// E2: y^2 = X(X+t^2)(X+4s+t^2)
// E3: y^2 = X(X^2 - 4(t^2+2s)X + 16s^2)
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();

function Curves(s,t)
    E2 := EllipticCurve([Q!0, Q!(2*t^2+4*s), Q!0,
                         Q!(t^2*(t^2+4*s)), Q!0]);
    E3 := EllipticCurve([Q!0, Q!(-4*(t^2+2*s)), Q!0,
                         Q!(16*s^2), Q!0]);
    return E2, E3;
end function;

s := Q!1;
t := Q!1;
if assigned s0 then
    s := Q!StringToRational(s0);
end if;
if assigned t0 then
    t := Q!StringToRational(t0);
end if;

E2, E3 := Curves(s,t);
print "s", s, "t", t;
print "E2", E2;
print "E3", E3;

E2q := IsogenyFromKernel(E2, Polynomial([Q!0,Q!1]));
E3q := IsogenyFromKernel(E3, Polynomial([Q!0,Q!1]));
print "E2/(0,0)", E2q;
print "E3/(0,0)", E3q;
print "E2 quotient isomorphic to E3", IsIsomorphic(E2q, E3);
print "E3 quotient isomorphic to E2", IsIsomorphic(E3q, E2);
print "isomorphism E3quot -> E2", Isomorphism(E3q, E2);

function DualPhiE3ToE2(R, s, E2)
    E3 := Curve(R);
    if R eq E3!0 then
        return E2!0;
    end if;
    X := R[1];
    Y := R[2];
    if X eq 0 then
        return E2!0;
    end if;
    return E2![Y^2/(4*X^2), Y*(16*s^2 - X^2)/(8*X^2)];
end function;

print "E2 small points", Points(E2 : Bound := 20);
print "E3 small points", Points(E3 : Bound := 20);

for R in Points(E3 : Bound := 20) do
    print "dual_phi", R, "->", DualPhiE3ToE2(R, s, E2);
end for;

quit;
