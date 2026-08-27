//////////////////////////////////////////////////////////////////////
// E8 Prym-Chabauty, stage 2: pair the Prym differentials against the
// transported generator; solve for the annihilating differential.
//
// Generator of J(C2'_min)(Q) (rank 1): Mumford (x^2+2x+4, 5x+5);
// support x = -1 +- sqrt(-3), y = +-5 sqrt(-3).  To C2'par: t = 2x+1,
// s = y/3.  Bigonal transport: xE4 = -72/(t^2+3), v=-72t/(t^2+3),
// u = 6(xE4^2-6xE4-36)v/(xE4^2(xE4+24)); w,w' roots of T^2-sT+u;
// yE4 = (w^2-A)/B, A = a/d, B = b/d (bigonal step-1 data);
// Tuitman coords: X = e = (xE4/3+2)/(w^2-xE4), Y = e*w.
// The generator class = [D1 - Dinf]; Dinf = C2' infinities transported
// (2 finite points via limits w' -> u/s, plus 2 boundary points over
// e=0).  Pairing(omega) = sum over D1-points of Int_b^P omega
//                       - sum over Dinf-points.
//////////////////////////////////////////////////////////////////////
load "coleman.m";
Q8 := y^8 + (216*x^4+72*x^3-24*x^2)*y^4 + (-1296*x^6-1728*x^5-432*x^4+64*x^3)*y^2
     + (-3888*x^8-2592*x^7+432*x^6+288*x^5-48*x^4);
if not assigned p then p := 37; elif Type(p) eq MonStgElt then p := StringToInteger(p); end if;
N := 20;
printf "RUNNING AT p = %o\n", p;
data := coleman_data(Q8, p, N);
printf "coleman_data ready, g=%o, N=%o\n", data`g, N;

Qp := pAdicField(p, N);
evalQ8 := func<X0, Y0 | &+[ Evaluate(ChangeRing(Coefficient(Q8, i), Qp), X0)*Y0^i : i in [0..8]]>;
// embedding z (zeta_3) -> Q_7: root of T^2+T+1 with z = 2 mod 7
Pq<T> := PolynomialRing(Qp);
rts := [r[1] : r in Roots(T^2+T+1)];
z7 := rts[1];   // at p=37 both deltas are squares under either embedding
lam := 2*z7+1;   // sqrt(-3)
printf "z -> %o (mod 7: 2), lam^2 = %o\n", z7, lam^2;

// generator support over K embedded in Q_7
xm := [-1+lam, -1-lam];
ts := [2*xm[1]+1, 2*xm[2]+1];
ss := [(5*xm[1]+5)/3, (5*xm[2]+5)/3];   // s-values on the SQUAREFREE model

// model-to-true s conversion: S(t)=2A+2u = f2sf*sq^2/den^2 ;
// s_true = s_model * sq(t)/den(t).  Recompute sq, den symbolically.
Kt<tt> := FunctionField(RationalField());
xt := -72/(tt^2+3); vt := -72*tt/(tt^2+3);
ut := 6*(xt^2-6*xt-36)*vt/(xt^2*(xt+24));
At := (xt^4/4 - 15*xt^3/2 + 324*xt + 972)/(xt^3*(xt+24)/24);
St := 2*At + 2*ut;
numS := Numerator(St); denS := Denominator(St);
f2t := numS*denS;
f2sft := &*[ fe[1]^(fe[2] mod 2) : fe in Factorization(f2t) ] * LeadingCoefficient(f2t);
sqt := Sqrt(f2t div (f2t div &*[fe[1]^(2*(fe[2] div 2)) : fe in Factorization(f2t)]))
       where _ := 0;
// simpler: sq^2 = f2t / f2sft
sq2 := f2t / f2sft;
Pu := Parent(numS);
assert Denominator(sq2) eq 1;
sqpoly := Sqrt(Numerator(sq2));
printf "den deg %o, sq deg %o, f2sf deg %o\n", Degree(denS), Degree(sqpoly), Degree(Numerator(f2sft));

// bigonal a,b,d (from step 1)
aP := func<X | X^4/4 - 15*X^3/2 + 324*X + 972>;
bP := func<X | X + 6>;
dP := func<X | X^3*(X+24)/24>;

pts_D1 := [];   // transported generator support: 4 points (X,Y) on Tuitman model
for i in [1,2] do
    t0 := ts[i]; s0m := ss[i];
    // convert model s to true s = w + w'
    sqv := &+[Qp | Qp!Coefficient(sqpoly,k)*t0^k : k in [0..Degree(sqpoly)]];
    denv := &+[Qp | Qp!Coefficient(denS,k)*t0^k : k in [0..Degree(denS)]];
    s0 := s0m * sqv/denv;
    x0 := -72/(t0^2+3);
    v0 := -72*t0/(t0^2+3);
    u0 := 6*(x0^2-6*x0-36)*v0/(x0^2*(x0+24));
    // consistency: s_true^2 = 2A + 2u
    A0 := (x0^4/4 - 15*x0^3/2 + 324*x0 + 972)/(x0^3*(x0+24)/24);
    printf "consistency s^2-(2A+2u) val: %o\n", Valuation(s0^2 - 2*A0 - 2*u0);
    dsc := s0^2 - 4*u0;
    sq := Sqrt(dsc);
    for sgn in [1,-1] do
        w0 := (s0 + sgn*sq)/2;
        e0 := (x0/3 + 2)/(w0^2 - x0);
        X0 := e0; Y0 := e0*w0;
        err := evalQ8(X0, Y0);
        printf "D1 point %o/%o: X val %o, model err val %o\n",
            i, sgn, Valuation(X0), Valuation(err);
        Append(~pts_D1, <X0, Y0>);
    end for;
end for;

// transported infinities of C2': as t -> inf on each sheet
// s ~ sgn*(lam/24) t^3 (lead of f2' = -1/192 = (lam/24)^2 * (-... check):
// (lam/24)^2 = -3/576 = -1/192  OK
// u ~ ? compute u(t) series: u = 6(x^2-6x-36) v/(x^2(x+24)), x=-72/(t^2+3):
// leading behavior: x ~ -72/t^2, v ~ -72/t
// x^2-6x-36 ~ -36; x^2 ~ 5184/t^4; x+24 ~ 24
// u ~ 6*(-36)*(-72/t) / ((5184/t^4)*24) = (6*36*72/(5184*24)) * t^3 = t^3/8
// so w ~ s ~ sgn(lam/24)t^3 -> inf (boundary pt over e=0),
//    w' ~ u/s = (t^3/8)/(sgn lam t^3/24) = 3/(sgn*lam) = sgn*3/lam = -sgn*lam ( 3/lam = -lam since lam^2=-3 )
// finite limit points: w'inf(sgn) = -sgn*lam, over xE4 -> 0:
// e = (x/3+2)/(w'^2 - x) -> 2/w'^2 = 2/(-3) = -2/3
// Y = e*w' = -2/3 * (-sgn lam) = sgn*2lam/3
print "transported infinity finite points: X=-2/3, Y=+-2*lam/3";
pts_Dinf_fin := [ <Qp!(-2)/3, 2*lam/3>, <Qp!(-2)/3, -2*lam/3> ];
for P in pts_Dinf_fin do
    err := evalQ8(P[1], P[2]);
    printf "Dinf finite pt: X=%o model err val %o\n", P[1], Valuation(err);
end for;
// boundary limit points (w -> inf): over e=0; identify with the rational
// boundary points found by Q_points (x=0).  These contribute integrals
// Int_b^{bd} omega; with b also boundary this is between the two x=0 points.
qpts := Q_points(data, 1000);
printf "rational points: %o\n", #qpts;

// set points and integrate basis[3], basis[4] (the odd/anti-invariant ones)
b := qpts[1];
rpoly := data`r;
W0m := data`W0;
allpts := [];
for P in pts_D1 cat pts_Dinf_fin do
    X0 := P[1]; Y0 := P[2];
    rv := &+[Qp | Qp!Coefficient(rpoly,k)*X0^k : k in [0..Degree(rpoly)]];
    if Valuation(rv) eq 0 then
        Append(~allpts, set_point(X0, Y0, data));
    else
        // bad disk: supply the integral-basis vector b_j = sum_k W0[j,k](X0) Y0^(k-1)
        bvec := [];
        for j in [1..8] do
            bj := Qp!0;
            for k in [1..8] do
                w0jk := W0m[j,k];
                if w0jk ne 0 then
                    nn := Numerator(w0jk); dd := Denominator(w0jk);
                    nv := &+[Qp | Qp!Coefficient(nn,h)*X0^h : h in [0..Degree(nn)]];
                    dv := &+[Qp | Qp!Coefficient(dd,h)*X0^h : h in [0..Degree(dd)]];
                    bj +:= (nv/dv) * Y0^(k-1);
                end if;
            end for;
            Append(~bvec, bj);
        end for;
        printf "bad-disk point at X val %o: using set_bad_point\n", Valuation(X0);
        Append(~allpts, set_bad_point(X0, bvec, false, data));
    end if;
end for;
// integrals from b to each point
Ints := [];
for i in [1..#allpts] do
    II, NII := coleman_integrals_on_basis(b, allpts[i], data : e := 100);
    printf "point %o: integrals (basis 3,4) = %o , %o  (prec %o)\n",
        i, II[3], II[4], NII;
    Append(~Ints, II);
end for;
// integral between the two boundary rational points (for Dinf boundary part)
IIb, NIIb := coleman_integrals_on_basis(b, qpts[2], data : e := 300);
printf "b -> other boundary point: basis3,4 = %o, %o\n", IIb[3], IIb[4];

// pairing: <omega_k, Z_* gen> = sum_{D1} - sum_{Dinf}
// Dinf = 2 finite pts + 2 boundary pts (assumed = the two rational pts;
// their contribution: Int_b^{b}=0 and Int_b^{qpts[2]})
for k in [3,4] do
    pair := &+[Ints[i][k] : i in [1..4]] - Ints[5][k] - Ints[6][k] - IIb[k];
    printf "PAIRING basis[%o] = %o\n", k, pair;
end for;

// ================= STAGE 3: the Chabauty sweep =================
// omega_ann = pair4 * basis[3] - pair3 * basis[4]; scale out common 37
pair3 := &+[Ints[i][3] : i in [1..4]] - Ints[5][3] - Ints[6][3] - IIb[3];
pair4 := &+[Ints[i][4] : i in [1..4]] - Ints[5][4] - Ints[6][4] - IIb[4];
mv := Min(Valuation(pair3), Valuation(pair4));
c3 := pair4/37^mv; c4 := -pair3/37^mv;
printf "omega_ann coefficients (basis3, basis4): %o , %o\n", c3, c4;
QpN := pAdicField(p, data`N);
v := [Vector(QpN, [0, 0, QpN!c3, QpN!c4])];

// local coordinate setup for the known rational points (as in effective_chabauty)
ee := 300;
for i:=1 to #qpts do
    _,index:=local_data(qpts[i],data);
    data:=update_minpolys(data,qpts[i]`inf,index);
    if is_bad(qpts[i],data) then
        if is_very_bad(qpts[i],data) then
            xt,bt,index:=local_coord(qpts[i],tadicprec(data,ee),data);
            qpts[i]`xt:=xt; qpts[i]`bt:=bt; qpts[i]`index:=index;
        end if;
    else
        xt,bt,index:=local_coord(qpts[i],tadicprec(data,1),data);
        qpts[i]`xt:=xt; qpts[i]`bt:=bt; qpts[i]`index:=index;
    end if;
end for;
Qppoints, data := Qp_points(data : points := qpts);
printf "residue disks to sweep: %o\n", #Qppoints;
for i:=1 to #Qppoints do
    if is_bad(Qppoints[i],data) then
        xt,bt,index:=local_coord(Qppoints[i],tadicprec(data,ee),data);
    else
        xt,bt,index:=local_coord(Qppoints[i],tadicprec(data,1),data);
    end if;
    Qppoints[i]`xt:=xt; Qppoints[i]`bt:=bt; Qppoints[i]`index:=index;
end for;
zerocount := 0;
zeroxs := [];
for i:=1 to #Qppoints do
    tdisk := Cputime();
    pts := zeros_on_disk(qpts[1], Qppoints[i], v, data : e:=ee);
    printf "disk %o/%o swept (%.1os): %o zero(s)\n", i, #Qppoints, Cputime(tdisk), #pts;
    if #pts gt 0 then
        for P in pts do
            zerocount +:= 1;
            printf "ZERO in disk %o : x = %o , inf = %o\n", i, P`x, P`inf;
            Append(~zeroxs, P`x);
        end for;
    end if;
end for;
printf "TOTAL ZEROS of the Prym-Chabauty function: %o\n", zerocount;

// ============ post-processing: rational reconstruction ============
// Any zero whose x-coordinate is a small rational e0 (with e0 ne 0) is a
// candidate NEW rational point of E8 => live [6,12] candidate: test the
// fiber exactly over Q.
print "=== rational reconstruction of zeros ===";
QQ := Rationals(); Pw := PolynomialRing(QQ);
for i in [1..#zeroxs] do
    x0 := zeroxs[i];
    if Valuation(x0) lt 0 then printf "zero %o : negative valuation, skip\n", i; continue; end if;
    k := Minimum(Precision(Parent(x0)), AbsolutePrecision(x0));
    Zp := IntegerRing(Parent(x0));
    ii := Integers()!(Zp!x0);
    ok, q := RationalReconstruction(IntegerRing(p^k)!ii);
    if ok then
        printf "zero %o reconstructs: e0 = %o", i, q;
        if q ne 0 then
            // exact fiber test: rational roots of Q8(e0, y)
            fib := Pw ! [Evaluate(Coefficient(Q8, j), q) : j in [0..8]];
            rr := Roots(fib);
            if #rr gt 0 then
                printf "   !!! RATIONAL POINT(S) ON E8: y = %o  [6,12] CANDIDATE !!!", [r[1] : r in rr];
            else
                printf "   (no rational fiber point: mock zero)";
            end if;
        else
            printf "   (boundary e=0)";
        end if;
        printf "\n";
    else
        printf "zero %o : x irrational (mock zero at this precision)\n", i;
    end if;
end for;
print "E8_COLEMAN3_DONE";
quit;
