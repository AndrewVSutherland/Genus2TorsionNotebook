// probe_ps.m — compare Points(E : Bound := B) vs PointSearch(E, B) on a
// representative big-coefficient fiber curve, and calibrate PointSearch's
// bound semantics (does it bound H(x) like Points?).
SetColumns(0);
SetClassGroupBounds("GRH");
E := EllipticCurve([0,0,1,-1,0]);   // 37a: known points, cheap
t0 := Realtime();
p1 := Points(E : Bound := 200);
printf "37a Points(200): %o pts %o s\n", #p1, Realtime()-t0;
ok := true;
try
    t0 := Realtime();
    p2 := PointSearch(E, 200);
    printf "37a PointSearch(200): %o pts %o s; sample %o\n", #p2, Realtime()-t0, [ P[1] : P in p2[1..Min(#p2,5)] ];
catch e
    printf "PointSearch failed: %o\n", e`Object;
    ok := false;
end try;
// a heavy fiber curve (si2 u=-22 scale): rebuild quickly
RQ := Rationals();
Qx<x> := PolynomialRing(RQ);
function AV(a, b, xv)
    if a lt b then
        if a eq 1 then
            return b eq 2 select (xv+3)*(xv-5) else 2*(xv-3);
        else
            return -(xv-1)*(xv-9);
        end if;
    else
        return -AV(b, a, xv);
    end if;
end function;
sg := [2,1,3];
u0 := RQ! -16;
d := 2*AV(sg[1], sg[3], u0); cA := AV(sg[1], sg[2], u0);
quart := cA*(x^2+6*d)*(x^2-2*d);
C := HyperellipticCurve(quart);
pts := Points(C : Bound := 200);
aff := [ P : P in pts | P[3] ne 0 ];
E2, mE := EllipticCurve(C, aff[1]);
Em, phi := MinimalModel(E2);
printf "fiber curve Em: %o\n", aInvariants(Em);
t0 := Realtime();
q1 := Points(Em : Bound := 101470);
printf "fiber Points(101470): %o pts %o s\n", #q1, Realtime()-t0;
if ok then
    t0 := Realtime();
    q2 := PointSearch(Em, 101470);
    printf "fiber PointSearch(101470): %o pts %o s\n", #q2, Realtime()-t0;
    // do the searches agree on affine points?
    s1 := { P[1] : P in q1 | P[3] ne 0 };
    s2 := { P[1] : P in q2 | P[3] ne 0 };
    printf "x-sets: Points %o, PointSearch %o, Points-minus-PS %o\n", #s1, #s2, #(s1 diff s2);
    t0 := Realtime();
    q3 := PointSearch(Em, 2000000);
    printf "fiber PointSearch(2e6): %o pts %o s\n", #q3, Realtime()-t0;
end if;
printf "PROBE_PS_DONE\n";
quit;
