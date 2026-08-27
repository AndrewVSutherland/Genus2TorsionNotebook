
//////////////////////////////////////////////////////////////////////
//  Boundary equation solver for the cleared second-halving systems.
//
//  This works over small finite fields only.  It evaluates the cleared
//  polynomial closures of FIRST_COVER, TARGET_416, and TARGET_88 on base
//  boundary residues and checks whether the corresponding auxiliary
//  variables exist over F_p.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

prime_list := [3,5,7];
if assigned primes then
    if Type(primes) eq MonStgElt then
        prime_list := [StringToInteger(s) : s in Split(primes, ",")];
    else
        prime_list := primes;
    end if;
end if;

max_samples := 8;
Z := Integers();

procedure Increment(~A, key)
    if IsDefined(A, key) then
        A[key] +:= 1;
    else
        A[key] := 1;
    end if;
end procedure;

procedure AddSample(~A, key, pt)
    if not IsDefined(A, key) then
        A[key] := [];
    end if;
    if #A[key] lt max_samples then
        Append(~A[key], pt);
    end if;
end procedure;

function BoundaryLabels(F, R, w)
    labels := [];
    if R eq 0 then Append(~labels, "R"); end if;
    if w eq 0 then Append(~labels, "w"); end if;
    if w - 1 eq 0 then Append(~labels, "w-1"); end if;
    if w + 1 eq 0 then Append(~labels, "w+1"); end if;
    if R - 1 eq 0 then Append(~labels, "R-1"); end if;
    if R + 1 eq 0 then Append(~labels, "R+1"); end if;
    if R - w eq 0 then Append(~labels, "R-w"); end if;
    if R + w eq 0 then Append(~labels, "R+w"); end if;
    if R*w - 3*R + 3*w - 1 eq 0 then Append(~labels, "Lplus"); end if;
    if R*w + 3*R + 3*w + 1 eq 0 then Append(~labels, "Lminus"); end if;
    if 2*R^2 - R*w^2 + R - 2*w^2 eq 0 then Append(~labels, "Q"); end if;
    if R^4 - 2*R^3 + R^2*w^2 - R^2 + 2*R*w^2 - w^2 eq 0 then
        Append(~labels, "Quartic");
    end if;
    return labels;
end function;

function BoundaryKey(labels)
    if #labels eq 0 then
        return "none";
    end if;
    return Join(Sort(labels), "&");
end function;

function BuildEquations(F)
    Rng<R,w,U,V,M,N,a,b,c,d,e,rho,sigma> := PolynomialRing(F, 13, "grevlex");
    K := FieldOfFractions(Rng);
    PX<x> := PolynomialRing(K);
    t := (2*R^2 + (1-w^2)*R - 2*w^2)/(4*(w^2-1));
    A := x^2 + (R^3 + 4*R^2*t + R - 8*R*t + 4*t)*x + R^4;
    B := (R + 2 + 4*t)*x^2 + (R^2 + 4*R + 1 + 8*t)*x
         + (2*R^2 + R + 4*t);
    f := x*A*B;
    h := A*B;
    c4 := Coefficient(h, 4);

    a0 := x^2 + U*x + V;
    F0 := h - x*(M*x+N)^2 - c4*a0^2;
    E0 := [Numerator(Coefficient(F0, i)) : i in [0..4]];

    q416 := x^2 + a*x + b;
    ell416 := c*x^2 + d*x + e;
    F416 := f - ell416^2 - c4*(x+R)*q416^2;
    E416 := [Numerator(Coefficient(F416, i)) : i in [0..4]];

    v0 := (M*U - N)*x + M*V;
    q88 := x^2 + a*x + b;
    ell88 := -v0 + a0*(rho*x + sigma);
    F88 := f - ell88^2 + rho^2*a0*q88^2;
    E88 := [Numerator(Coefficient(F88, i)) : i in [0..5]];

    return Rng, E0, E416, E88;
end function;

function AllZero(eqs, vals)
    for pol in eqs do
        if Evaluate(pol, vals) ne 0 then
            return false;
        end if;
    end for;
    return true;
end function;

function Has416(E416, elems, r0, w0)
    z := elems[1];
    for aa in elems do
        for bb in elems do
            for cc in elems do
                for dd in elems do
                    for ee in elems do
                        vals := [r0,w0,z,z,z,z,aa,bb,cc,dd,ee,z,z];
                        if AllZero(E416, vals) then
                            return true, <aa,bb,cc,dd,ee>;
                        end if;
                    end for;
                end for;
            end for;
        end for;
    end for;
    return false, _;
end function;

function HasFirstAnd88(E0, E88, elems, r0, w0)
    z := elems[1];
    has_first := false;
    first_sample := <>;
    for UU in elems do
        for VV in elems do
            for MM in elems do
                for NN in elems do
                    vals0 := [r0,w0,UU,VV,MM,NN,z,z,z,z,z,z,z];
                    if not AllZero(E0, vals0) then
                        continue;
                    end if;
                    has_first := true;
                    if #first_sample eq 0 then
                        first_sample := <UU,VV,MM,NN>;
                    end if;
                    for aa in elems do
                        for bb in elems do
                            for rr in elems do
                                for ss in elems do
                                    vals88 := [r0,w0,UU,VV,MM,NN,aa,bb,z,z,z,rr,ss];
                                    if AllZero(E88, vals88) then
                                        return true, true, first_sample, <UU,VV,MM,NN,aa,bb,rr,ss>;
                                    end if;
                                end for;
                            end for;
                        end for;
                    end for;
                end for;
            end for;
        end for;
    end for;
    return has_first, false, first_sample, _;
end function;

for p in prime_list do
    if p eq 2 then
        print "p", p, "skipped_char2";
        continue;
    end if;
    F := GF(p);
    elems := [x : x in F];
    Rng, E0, E416, E88 := BuildEquations(F);

    total_boundary := 0;
    first_count := 0;
    target416_count := 0;
    target88_count := 0;
    exact_total := AssociativeArray();
    exact_first := AssociativeArray();
    exact_416 := AssociativeArray();
    exact_88 := AssociativeArray();
    samples_first := AssociativeArray();
    samples_416 := AssociativeArray();
    samples_88 := AssociativeArray();

    print "BOUNDARY_EQUATION_SOLVER p", p;
    for r0 in F do
        for w0 in F do
            labels := BoundaryLabels(F, r0, w0);
            if #labels eq 0 then
                continue;
            end if;
            key := BoundaryKey(labels);
            total_boundary +:= 1;
            Increment(~exact_total, key);

            has_first, has_88, first_sample, sol88 := HasFirstAnd88(E0, E88, elems, r0, w0);
            if has_first then
                first_count +:= 1;
                Increment(~exact_first, key);
                AddSample(~samples_first, key, <Z!r0,Z!w0,first_sample>);
            end if;
            ok416, sol416 := Has416(E416, elems, r0, w0);
            if has_first and ok416 then
                target416_count +:= 1;
                Increment(~exact_416, key);
                AddSample(~samples_416, key, <Z!r0,Z!w0,sol416>);
            end if;
            if has_88 then
                target88_count +:= 1;
                Increment(~exact_88, key);
                AddSample(~samples_88, key, <Z!r0,Z!w0,sol88>);
            end if;
        end for;
    end for;

    print "summary", "boundary", total_boundary,
          "first_closure", first_count,
          "target416_closure", target416_count,
          "target88_closure", target88_count;
    print "strata";
    for key in Sort([k : k in Keys(exact_total)]) do
        fcnt := IsDefined(exact_first, key) select exact_first[key] else 0;
        c416 := IsDefined(exact_416, key) select exact_416[key] else 0;
        c88 := IsDefined(exact_88, key) select exact_88[key] else 0;
        print " ", key, "total", exact_total[key], "first", fcnt,
              "416", c416, "88", c88;
        if fcnt gt 0 then print "    first_samples", samples_first[key]; end if;
        if c416 gt 0 then print "    416_samples", samples_416[key]; end if;
        if c88 gt 0 then print "    88_samples", samples_88[key]; end if;
    end for;
end for;

quit;
