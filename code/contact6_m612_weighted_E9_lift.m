//////////////////////////////////////////////////////////////////////
// The next weighted layer of the exceptional E9 stratum above
// [A:B:T]=[1:0:0] for the contact-6 [6,12] core, at p=5.
//
// On the exact endpoint b=0, a=1/e, use the E9 scaling
//
//   e=5^2*E, L=5^-1*l, U=5^-2*u, v=5^-2*w,
//   N=5^-1*n, R=5^-4*r.
//
// The script derives the five integral normalized equations, enumerates
// the 16 contact-open E9 points over F_5, and solves the first correction
// equations exactly.  Direct evaluation modulo 25 verifies every count.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Z:=Integers();
Q:=Rationals();
p:=5;
k:=GF(p);

A<E,l,u,w,n,r>:=PolynomialRing(Q,6,"grevlex");

// These are respectively 5^2*E5, E*5^4*E4, 5^6*E3,
// E^2*5^8*E2, and E*5^10*E1 after the E9 substitution.  Multiplication
// by E or E^2 only clears a unit denominator.
F5:=2*p*n-3*u-6*l^2;
F4:=E*(p^2*n^2+2*r-3*u^2-3*w^2+15*p^2*l^2)-2*l^2;
F3:=2*w^3+2*p*n*r-u^3-6*u*w^2-22*p^4*l^2;
F2:=E^2*(r^2+2*p*n*w^3-3*u^2*w^2-3*w^4
          +15*p^6*l^2)-p^2*l^2;
F1:=E*(2*r*w^3-3*u*w^4-6*p^8*l^2)-2*p^6*l^2;
Fs:=[F5,F4,F3,F2,F1];

// Verify the normalization against the exact five coefficient equations.
K:=FieldOfFractions(A);
ee:=K!(p^2*E);
LL:=K!l/p;
UU:=K!u/p^2;
WW:=K!w/p^2;
NN:=K!n/p;
RR:=K!r/p^4;
c5:=K!6;
c4:=2/ee-15;
c3:=K!22;
c2:=1/ee^2-15;
c1:=2/ee+6;
O5:=2*NN-3*UU-c5*LL^2;
O4:=NN^2+2*RR-3*UU^2-3*WW^2-c4*LL^2;
O3:=2*WW^3+2*NN*RR-UU^3-6*UU*WW^2-c3*LL^2;
O2:=RR^2+2*NN*WW^3-3*UU^2*WW^2-3*WW^4-c2*LL^2;
O1:=2*RR*WW^3-3*UU*WW^4-c1*LL^2;
assert K!F5 eq p^2*O5;
assert K!F4 eq E*p^4*O4;
assert K!F3 eq p^6*O3;
assert K!F2 eq E^2*p^8*O2;
assert K!F1 eq E*p^10*O1;

function EvalZ(g,pt)
    z:=Evaluate(g,[Q!a:a in pt]);
    assert Denominator(z) eq 1;
    return Z!z;
end function;

print "CONTACT6_M612_WEIGHTED_E9_LIFT";
print "SCALING", "e=25E L=l/5 U=u/25 v=w/25 N=n/5 R=r/625";
print "NORMALIZED_EXACT_EQUATIONS",Fs;
print "MOD25_EQUATIONS";
print " F5=-3u-6l^2+10n";
print " F4=E*(2r-3u^2-3w^2)-2l^2";
print " F3=2w^3-u^3-6uw^2+10nr";
print " F2=E^2*(r^2-3u^2w^2-3w^4+10nw^3)";
print " F1=E*(2rw^3-3uw^4)";

// Recompute the special-fiber points rather than hard-coding the list.
initial_points:=[];
for e0 in [1..p-1] do for l0 in [1..p-1] do
for u0 in [1..p-1] do for w0 in [1..p-1] do
for n0 in [1..p-1] do for r0 in [1..p-1] do
    pt:=[e0,l0,u0,w0,n0,r0];
    if (u0^2-4*w0^2) mod p eq 0 then continue; end if;
    if &and[(EvalZ(f,pt) mod p) eq 0:f in Fs] then
        Append(~initial_points,pt);
    end if;
end for; end for; end for; end for; end for; end for;

assert #initial_points eq 16;
print "INITIAL_POINT_COUNT",#initial_points;
print "INITIAL_POINTS_E_l_u_w_n_r",initial_points;

survivors:=[];
total_lifts:=0;
for pt in initial_points do
    // F(pt+5*d)/5 = F(pt)/5 + Jac(F mod 5)*d (mod 5).
    J:=Matrix(k,5,6,
        &cat[[k!EvalZ(Derivative(Fs[i],j),pt):j in [1..6]]
             :i in [1..5]]);
    rhs:=Vector(k,[k!((-EvalZ(Fs[i],pt) div p) mod p):i in [1..5]]);
    aug:=HorizontalJoin(J,Matrix(k,5,1,[rhs[i]:i in [1..5]]));
    consistent:=Rank(aug) eq Rank(J);

    lift_count:=0;
    correction_samples:=[];
    for dE in [0..p-1] do for dl in [0..p-1] do
    for du in [0..p-1] do for dw in [0..p-1] do
    for dn in [0..p-1] do for dr in [0..p-1] do
        ds:=[dE,dl,du,dw,dn,dr];
        if not &and[
            &+[J[i,j]*k!ds[j]:j in [1..6]] eq rhs[i]
            :i in [1..5]
        ] then
            continue;
        end if;
        vals:=[pt[j]+p*ds[j]:j in [1..6]];
        assert &and[(EvalZ(f,vals) mod p^2) eq 0:f in Fs];
        lift_count+:=1;
        if #correction_samples lt 3 then
            Append(~correction_samples,<ds,vals>);
        end if;
    end for; end for; end for; end for; end for; end for;

    assert consistent eq (lift_count gt 0);
    if consistent then
        assert lift_count eq p^(6-Rank(J));
        Append(~survivors,pt);
        total_lifts+:=lift_count;
    end if;

    // On all 16 rational E9 points the scalar next-layer obstruction is,
    // up to a nonzero residue factor, n-3E.  The rank comparison and the
    // direct mod-25 enumeration independently verify this simplification.
    obstruction:=k!(pt[5]-3*pt[1]);
    assert consistent eq (obstruction eq 0);
    print "POINT",pt,"JAC_RANK",Rank(J),"NEXT_RHS",Eltseq(rhs),
          "n_MINUS_3E",Z!obstruction,"LIFTS_MOD25",lift_count;
    if lift_count gt 0 then
        print " CORRECTION_AND_RESIDUE_SAMPLES",correction_samples;
    end if;
end for;

assert survivors eq [
    [1,2,2,2,3,1], [1,3,2,2,3,1],
    [4,1,3,3,2,1], [4,4,3,3,2,1]
];
assert total_lifts eq 100;

print "NEXT_WEIGHTED_OBSTRUCTION", "n-3E=0 on E9(F_5)";
print "SURVIVING_INITIAL_COUNT",#survivors;
print "SURVIVING_INITIAL_POINTS",survivors;
print "TOTAL_MOD25_LIFTS",total_lifts;
print "VERDICT", "4 of 16 initial points lift through the next weighted layer";
print "CONTACT6_M612_WEIGHTED_E9_LIFT_DONE";
quit;
