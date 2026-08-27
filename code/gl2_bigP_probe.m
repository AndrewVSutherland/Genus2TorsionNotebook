// 37-hunt: fast convention probe on the dumped big period matrices.
// Loads results/gl2_bigP_dump.m (12 tuples <combo, Eltseq(bigP) row-major>),
// tries the finite set of block conventions (identity, block swap, sign
// flips) against IsBigPeriodMatrix, and reconstructs with the first valid
// variant per candidate.  Base := true forces Q-recognition.
//
// Run: magma -b code/gl2_bigP_probe.m > results/gl2_bigP_probe.log

SetColumns(0);
SetSeed(1);
SetMemoryLimit(12*10^9);

Attach("~/.claude/jobs/a1db5dd4/tmp/polredabs_shim.m");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/endomorphisms/endomorphisms/magma/spec");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/quartic/magma/spec");
AttachSpec("~/.claude/jobs/a1db5dd4/tmp/curve_reconstruction/magma/spec");

Prec := 140;
FX := RationalsExtra(Prec);
CC := FX`CC;

// preprocessed by python into results/gl2_bigP_clean.m:
// cands := [ <[c1..c4], [re1,im1,...,re8,im8]>, ... ];
load "results/gl2_bigP_clean.m";
dumps := [* *];
RR := RealField(Prec);
for t in cands do
    es := [ CC!(RR!t[2][2*i-1]) + CC.1*(RR!t[2][2*i]) : i in [1..8] ];
    Append(~dumps, < t[1], es >);
end for;
printf "LOADED %o dumped matrices\n", #dumps;

function Var(P, which)
    P1 := Submatrix(P, 1,1, 2,2); P2 := Submatrix(P, 1,3, 2,2);
    if which eq 1 then return HorizontalJoin(P1, P2);
    elif which eq 2 then return HorizontalJoin(P2, P1);
    elif which eq 3 then return HorizontalJoin(-P1, P2);
    elif which eq 4 then return HorizontalJoin(P2, -P1);
    elif which eq 5 then return HorizontalJoin(P1, -P2);
    elif which eq 6 then return HorizontalJoin(-P2, P1);
    else return P;
    end if;
end function;

for d in dumps do
    combo := d[1]; es := d[2];
    P := Matrix(CC, 2, 4, [ CC!x : x in es ]);
    found := false;
    for w in [1..6] do
        Pv := Var(P, w);
        isb := false;
        try
            isb := IsBigPeriodMatrix(Pv);
        catch e ; end try;
        if not isb then continue; end if;
        printf "COMBO %o variant %o : big period matrix OK\n", combo, w;
        found := true;
        try
            Crec := ReconstructCurve(Pv, FX : Base := true);
            printf "QCURVE combo=%o : %o\n", combo, Crec;
        catch e
            msg := Sprint(e`Object);
            printf "RECFAIL combo=%o variant=%o err=%o\n", combo, w,
                   #msg gt 200 select Substring(msg,1,200) else msg;
        end try;
        break;
    end for;
    if not found then printf "NOVARIANT combo=%o\n", combo; end if;
end for;
printf "PROBE_DONE\n";
quit;
