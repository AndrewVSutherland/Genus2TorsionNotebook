// recon.m — Lane B reconnaissance: for each live T5 member and each quartic
// cover (cv=3: y^2=W3, cv=4: y^2=W4), build integral quartic model, take a
// point at infinity as origin, pass to MinimalModel FIRST, then RankBounds.
// Known cover points (hits/nears) are checked to lie on the cover.
SetClassGroupBounds("GRH");
P<u> := PolynomialRing(Rationals());
q := 4*u^2 - 6*u + 3;

// <rn, rd, known-u-list on cover3, known-u-list on cover3&4 (full hits)>
members := [
  <-49, 240, [-97/48], [-97/48]>,
  <289, 240, [133/145], [133/145]>,
  <-1, 143, [17, 13/4], []>,
  <-25, 551, [10, 41/25], []>,
  <-169, 1431, [-4], []>,
  <841, 697, [43/52], []>
];

for mem in members do
    r5 := mem[1]/mem[2];
    printf "=========== member rho' = %o/%o ===========\n", mem[1], mem[2];
    W1 := r5*(q*r5-(2*u-1));
    W2 := q*r5^2-(4*u^2-4*u+2)*r5+(2*u-1);
    W3 := r5*((16*u^4-40*u^3+40*u^2-18*u+3)*r5^3+(-16*u^4+32*u^3-28*u^2+10*u-1)*r5^2
            +(8*u^3-12*u^2+10*u-3)*r5+(-2*u+1));
    W4 := (16*u^4-40*u^3+40*u^2-18*u+3)*r5^2+(-16*u^3+28*u^2-18*u+4)*r5+(4*u^2-4*u+1);
    covs := [<3, W3>, <4, W4>];
    for cvp in covs do
        cv := cvp[1]; Wc := cvp[2];
        dc := LCM([Denominator(co) : co in Coefficients(Wc)]);
        F := dc^2*Wc;
        // clear square content
        cont := GCD([Integers()!co : co in Coefficients(F)]);
        sq := 1;
        for pe in Factorization(cont) do
            sq *:= pe[1]^(2*(pe[2] div 2));
        end for;
        F := F div sq;
        printf "-- cover %o: F = %o\n", cv, F;
        if Degree(F) ne 4 or Discriminant(F) eq 0 then
            printf "   DEGENERATE quartic, skip\n"; continue;
        end if;
        lc := LeadingCoefficient(F);
        oksq, _ := IsSquare(lc);
        printf "   lc = %o, square: %o\n", lc, oksq;
        C := HyperellipticCurve(F);
        ptsI := PointsAtInfinity(C);
        printf "   points at infinity: %o\n", #ptsI;
        smp := Points(C : Bound := 10000);
        printf "   small points (H<=1e4): %o  e.g. %o\n", #smp,
            [smp[i] : i in [1..Min(#smp, 6)]];
        // known cover points sanity
        kn := cv eq 3 select mem[3] else mem[4];
        for u0 in kn do
            okk, _ := IsSquare(Evaluate(F, u0));
            printf "   known u=%o on cover: %o\n", u0, okk;
        end for;
        if #ptsI eq 0 then
            if #smp eq 0 then printf "   no origin found, skip\n"; continue; end if;
            org := Rep(smp);
        else
            org := Rep(ptsI);
        end if;
        E, mp := EllipticCurve(C, org);
        Em, phi := MinimalModel(E);
        printf "   Em = %o cond %o tors %o\n", aInvariants(Em),
            Conductor(Em), Invariants(TorsionSubgroup(Em));
        t0 := Cputime();
        rlo, rhi := RankBounds(Em);
        printf "   RankBounds [%o,%o]  (%.1o s)\n", rlo, rhi, Cputime(t0);
    end for;
end for;
quit;
