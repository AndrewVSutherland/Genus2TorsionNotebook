// 37-hunt endgame for the Bcontrol_d2 member: Mestre model from the exact
// rational G2 invariants, staged reduction (prints between steps), twist
// match to 2190.2.a.v traces, exact torsion.
SetColumns(0);
SetSeed(1);
SetMemoryLimit(24*10^9);
Q := Rationals();
g2Q := [ -162405241570776167323105/241097315486732,
         -7901191237270473674525/304422749567068,
         -4673522090926046225269/3440903425783038 ];
trtargets := [<7,2>, <11,-6>, <13,-2>, <17,6>, <19,10>];
printf "MESTRE start\n";
CQ := HyperellipticCurveFromG2Invariants(g2Q);
printf "MESTRE done over %o\n", BaseRing(CQ);
error if Type(BaseRing(CQ)) ne FldRat, "Mestre obstruction";
fC, hC := HyperellipticPolynomials(CQ);
gC := 4*fC + hC^2;
printf "RAW_COEFF_DIGITS %o\n", Max([#Sprint(Numerator(c)) : c in Coefficients(gC)]);
// clear denominators to integral binary sextic
den := LCM([Denominator(c) : c in Coefficients(gC)]);
gZ := gC * den^2;  // keep square class
// remove square content
cont := GCD([Numerator(c) : c in Coefficients(gZ)]);
sq := 1;
for pf in Factorization(cont) do
    sq *:= pf[1]^(2*(pf[2] div 2));
end for;
gZ := gZ / sq;
printf "INT_COEFF_DIGITS %o\n", Max([#Sprint(Numerator(c)) : c in Coefficients(gZ)]);
// Stoll reduction of the binary form via curve
C0 := HyperellipticCurve(gZ);
printf "REDUCE start\n";
Cred := C0;
try
    Cred := ReducedModel(C0);
    printf "REDUCE done (Stoll)\n";
catch e
    printf "REDUCE unavailable/failed: %o\n", e`Object;
end try;
fR, hR := HyperellipticPolynomials(Cred);
gR := 4*fR + hR^2;
printf "RED_COEFF_DIGITS %o\n", Max([#Sprint(Numerator(c)) : c in Coefficients(gR)]);
printf "RED_MODEL y^2 = %o\n", gR;
// twist scan on the reduced model
for dt in [1,-1,2,-2,3,-3,5,-5,6,-6,10,-10,15,-15,30,-30,73,-73,146,-146,219,-219,365,-365,438,-438,730,-730,1095,-1095,2190,-2190] do
    Cd := HyperellipticCurve(dt*gR);
    match := true;
    nch := 0;
    for tt in trtargets do
        p := tt[1];
        lpok := false;
        try
            Cp := ChangeRing(Cd, GF(p));
            chi := Reverse(Coefficients(LPolynomial(Cp)));
            lpok := true;
        catch e ; end try;
        if not lpok then continue; end if;   // bad reduction at p: skip
        nch +:= 1;
        if -Integers()!chi[2] ne tt[2] then match := false; break; end if;
    end for;
    if match and nch ge 3 then
        printf "TWIST_MATCH d=%o (checked %o primes)\n", dt, nch;
        printf "MINIMIZE start\n";
        Cmin := Cd;
        try
            Cmin := MinimalWeierstrassModel(Cd);
            printf "MINIMIZE done\n";
        catch e printf "MINIMIZE failed, using unminimized\n"; end try;
        Cs := SimplifiedModel(Cmin);
        printf "MATCHED_CURVE %o\n", Cs;
        printf "TORSION start\n";
        Tt := TorsionSubgroup(Jacobian(Cs));
        printf "TORSION %o (order %o)\n", Invariants(Tt), #Tt;
        if #Tt mod 37 eq 0 then
            printf "*** THEOREM: genus-2 Jacobian over Q with a rational point of order 37 ***\n";
            printf "CURVE %o\n", Cmin;
        end if;
        break;
    end if;
end for;
printf "ENDGAME_DONE\n";
quit;
