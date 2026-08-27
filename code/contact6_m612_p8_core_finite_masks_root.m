//////////////////////////////////////////////////////////////////////
// Finite-field parameter masks for the P8-pulled-back cubic-contact
// cover with L retained.  This is a fast sieve/diagnostic, not a global
// component decomposition.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned primes then primes:=[7,11,13,17,19,23,29,31]; end if;

print "CONTACT6_M612_P8_CORE_FINITE_MASKS_ROOT";
for p in primes do
    k:=GF(p);
    if p in {2,3,5} then continue; end if;
    PR<q>:=PolynomialRing(k);
    tnum:=k!4*(q^2+q-k!6); tden:=q^2+k!6;
    enum:=-(k!25/k!3)*tnum^2*tden^2;
    eden:=tnum^4-k!25*tnum^2*tden^2+(k!1250/k!3)*tden^4;

    masks:=[]; total_points:=0; parameter_poles:=[];
    for q0 in k do
        de:=Evaluate(eden,q0);
        if de eq 0 then Append(~parameter_poles,Integers()!q0); continue; end if;
        ee:=Evaluate(enum,q0)/de;
        if ee eq 0 or k!1+k!2*ee eq 0 then continue; end if;
        count:=0;
        for LL in k do
            if LL eq 0 then continue; end if;
            M:=LL^2;
            for UU in k do for vv in k do
                if vv eq 0 or UU^2 eq k!4*vv^2 then continue; end if;
                N:=(k!3*UU+k!6*M)/k!2;
                Rc:=(k!3*UU^2+k!3*vv^2+(k!2/ee-k!15)*M-N^2)/k!2;
                F3:=k!2*vv^3+k!2*N*Rc-UU^3-k!6*UU*vv^2-k!22*M;
                F2:=Rc^2+k!2*N*vv^3-k!3*UU^2*vv^2-k!3*vv^4
                    -(k!1/ee^2-k!15)*M;
                F1:=k!2*Rc*vv^3-k!3*UU*vv^4-(k!2/ee+k!6)*M;
                if F1 eq 0 and F2 eq 0 and F3 eq 0 then count+:=1; end if;
            end for; end for;
        end for;
        if count gt 0 then Append(~masks,<Integers()!q0,count>); end if;
        total_points+:=count;
    end for;

    // Projective q=infinity: take leading coefficients of e(q).
    einf:=LeadingCoefficient(enum)/LeadingCoefficient(eden);
    infinity_count:=0;
    if einf ne 0 and k!1+k!2*einf ne 0 then
        for LL in k do
            if LL eq 0 then continue; end if;
            M:=LL^2;
            for UU in k do for vv in k do
                if vv eq 0 or UU^2 eq k!4*vv^2 then continue; end if;
                N:=(k!3*UU+k!6*M)/k!2;
                Rc:=(k!3*UU^2+k!3*vv^2+(k!2/einf-k!15)*M-N^2)/k!2;
                F3:=k!2*vv^3+k!2*N*Rc-UU^3-k!6*UU*vv^2-k!22*M;
                F2:=Rc^2+k!2*N*vv^3-k!3*UU^2*vv^2-k!3*vv^4
                    -(k!1/einf^2-k!15)*M;
                F1:=k!2*Rc*vv^3-k!3*UU*vv^4-(k!2/einf+k!6)*M;
                if F1 eq 0 and F2 eq 0 and F3 eq 0 then
                    infinity_count+:=1;
                end if;
            end for; end for;
        end for;
    end if;
    print "PRIME",p,"FINITE_MASK_COUNT",#masks,
          "FINITE_MASKS",masks,"TOTAL_AFFINE_POINTS",total_points,
          "PARAMETER_POLES",parameter_poles,
          "INFINITY_POINTS",infinity_count;
end for;
print "CONTACT6_M612_P8_CORE_FINITE_MASKS_ROOT_DONE";
quit;
