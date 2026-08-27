// claude_census_2pk_audit.m
// For every displayed curve of the realization census (paper/torsion_realizations.tex),
// verify the exact rational torsion group and decide whether a torsion class of order
// equal to the largest invariant factor has the doubled-theta form 2P - K_C
// (P a rational point of C, including doubled rational points at infinity).
// Criterion: in the reduced Mumford representation <a,b,d> of a class g,
//   g = [2P - K]  <=>  d = 2 and ( deg a = 2 with disc(a) = 0, or deg a <= 0 ).
// (deg a <= 0 with d = 2 is a doubled point at infinity; the zero class never
// occurs because we only inspect elements of order n_r >= 2.)
SetColumns(0);
SetMemoryLimit(6*10^9);
QQ := Rationals();
R<x> := PolynomialRing(QQ);

procedure AuditRow(idx, grp, hpol, fpol)
    t0 := Cputime();
    ok, C := IsHyperellipticCurve([fpol, hpol]);
    if not ok then
        printf "ROW|%o|%o|ERR|not_hyperelliptic\n", idx, grp;
        return;
    end if;
    Cs := SimplifiedModel(C);
    Cs := IntegralModel(Cs);
    J := Jacobian(Cs);
    try
        T, mp := TorsionSubgroup(J);
    catch e
        printf "ROW|%o|%o|ERR|torsion_failed\n", idx, grp;
        return;
    end try;
    inv := Invariants(T);
    if #inv eq 0 then
        printf "ROW|%o|%o|inv=[]|nr=1|theta=[]|verdict=NO|t=%o\n", idx, grp, Cputime(t0);
        return;
    end if;
    nr := inv[#inv];
    thetaOrders := {};
    for a in T do
        if IsIdentity(a) then continue; end if;
        g := mp(a);
        apol := g[1]; d := g[3];
        if d eq 2 and ((Degree(apol) eq 2 and Discriminant(apol) eq 0)
                       or Degree(apol) le 0) then
            Include(~thetaOrders, Order(a));
        end if;
    end for;
    verdict := (nr in thetaOrders) select "YES" else "NO";
    printf "ROW|%o|%o|inv=%o|nr=%o|theta=%o|verdict=%o|t=%o\n",
        idx, grp, inv, nr, Sort(SetToSequence(thetaOrders)), verdict, Cputime(t0);
end procedure;


AuditRow(1, "2", R!((x^2+x+1)), R!(x^5-40*x^3+22*x^2+389*x-608));
AuditRow(2, "3", R!(1), R!(x^5-2*x^4+2*x^3-x^2));
AuditRow(3, "4", R!((x^3+1)), R!(x^3-x));
AuditRow(4, "5", R!(1), R!(x^5-9*x^4+14*x^3-19*x^2+11*x-6));
AuditRow(5, "6", R!((x^2+x)), R!(x^5-12*x^4+26*x^3+46*x^2+21*x+3));
AuditRow(6, "7", R!(x^3), R!(x^5-3*x^3+3*x-2));
AuditRow(7, "8", R!((x+1)), R!(-x^6-2*x^5-2*x^4-x^3));
AuditRow(8, "9", R!((x^3+x+1)), R!(-x^4));
AuditRow(9, "10", R!((x^3+x)), R!(x^5-2*x^4-8*x^3+16*x+7));
AuditRow(10, "11", R!((x^3+x+1)), R!(x^2));
AuditRow(11, "12", R!((x^3+x^2+x)), R!(x^2+x+1));
AuditRow(12, "13", R!((x^3+x^2+x+1)), R!(-x^3-x^2));
AuditRow(13, "14", R!((x^3+1)), R!(x^2+x));
AuditRow(14, "15", R!((x^3+x^2+x+1)), R!(-x^2-x));
AuditRow(15, "16", R!((x^3+1)), R!(-x^5+x^4-2*x^2+x+1));
AuditRow(16, "17", R!((x^3+x^2+x)), R!(x^3+x^2+3*x+1));
AuditRow(17, "18", R!((x^3+1)), R!(-2*x^4+4*x^2+2*x));
AuditRow(18, "19", R!((x^2+x+1)), R!(x^6+x^5+x^4-5*x^3-3*x^2+5*x-2));
AuditRow(19, "20", R!((x+1)), R!(-x^5));
AuditRow(20, "21", R!((x^3+x+1)), R!(-x^4+2*x^2+x));
AuditRow(21, "22", R!((x^3+x)), R!(x^3-2*x^2-x+1));
AuditRow(22, "23", R!(x^2), R!(x^6-x^5-2*x^4+x^3+2*x^2-2*x+1));
AuditRow(23, "24", R!((x^3+1)), R!(2*x^4+3*x^3+4*x^2+2*x));
AuditRow(24, "25", R!((x+1)), R!(-9*x^5+25*x^4+75*x^3+41*x^2-23*x+2));
AuditRow(25, "26", R!((x+1)), R!(-12*x^5-3*x^4-3*x^3-3*x^2));
AuditRow(26, "27", R!((x^3+1)), R!(-x^4+x^3+x^2-x));
AuditRow(27, "28", R!((x^3+1)), R!(-x^5+x^3+x^2+3*x+2));
AuditRow(28, "29", R!((x+1)), R!(x^6-2*x^5+2*x^3-x^2));
AuditRow(29, "30", R!((x^3+x^2)), R!(2*x^6+4*x^5-10*x^4-12*x^3+16*x^2-x));
AuditRow(30, "31", R!((-x^2-x)), R!(-839*x^6+2841*x^5-4587*x^4+4300*x^3-2466*x^2+816*x-126));
AuditRow(31, "32", R!((x^2+x)), R!(x^6-x^5+x^4+45*x^3+49*x^2-49*x-47));
AuditRow(32, "33", R!((x^2+x+1)), R!(x^6-7*x^5+9*x^4+7*x^3+21*x^2+13*x+2));
AuditRow(33, "34", R!((x^3+x^2)), R!(-2*x^5+3*x^3+7*x^2+12*x+15));
AuditRow(34, "36", R!(x^2), R!(3*x^5+30*x^4+13*x^3-14*x^2+6*x+9));
AuditRow(35, "39", R!((x^3+1)), R!(x^4+2*x^3+x^2-x));
AuditRow(36, "40", R!((x+1)), R!(3*x^5+4*x^4-2*x^3-2*x^2));
AuditRow(37, "2,2", R!(x), R!(4*x^5+33*x^4+72*x^3+16*x^2+x));
AuditRow(38, "2,4", R!(x), R!(x^5-8*x^4+16*x^3-x));
AuditRow(39, "2,6", R!(1), R!(4*x^5+4*x^4-x^3-2*x^2));
AuditRow(40, "2,8", R!((x+1)), R!(8*x^5+3*x^4-4*x^3-2*x^2));
AuditRow(41, "2,10", R!((x+1)), R!(3*x^5-2*x^4-4*x^3+x^2+x));
AuditRow(42, "2,12", R!((x^2+x)), R!(x^5-8*x^4+14*x^3+2*x^2-x));
AuditRow(43, "2,14", R!((x^3+x)), R!(-2*x^4-x^3+x+1));
AuditRow(44, "2,16", R!((x^3+x^2)), R!(-x^5-x^4+x^3-3*x^2+4*x+8));
AuditRow(45, "2,18", R!((x^2+x)), R!(x^6-x^5-5*x^4+3*x^3+x^2+2*x+7));
AuditRow(46, "2,20", R!(x^2), R!(4*x^6-30*x^5+28*x^4+185*x^3-251*x^2-180*x));
AuditRow(47, "2,22", R!((x^2+x)), R!(x^6-3*x^5+9*x^4-5*x^3+12*x^2-6*x));
AuditRow(48, "2,26", R!((x^2+1)), R!(9*x^6-3*x^5-12*x^4+3*x^3+2*x^2));
AuditRow(49, "2,28", R!((x^2+1)), R!(9*x^6-3*x^5-42*x^4+34*x^3+4*x^2-3*x));
AuditRow(50, "3,3", R!(1), R!(x^6+4*x^4+2*x^3+7*x^2+3*x+3));
AuditRow(51, "3,6", R!((x^3+x^2+x+1)), R!(-x^5-5*x^3+7*x^2-2*x+21));
AuditRow(52, "3,9", R!((x^3+1)), R!(16*x^6+75*x^5+40*x^4-113*x^3-230*x^2-75*x+326));
AuditRow(53, "4,4", R!((x^2+x)), R!(99*x^6+242*x^5+1823*x^4+2601*x^3+10832*x^2+6534*x+20493));
AuditRow(54, "4,8", R!((x^2+x)), R!(-81*x^5+272*x^4-435*x^3+450*x^2+682*x+164));
AuditRow(55, "6,6", R!(0), R!(x*(39*x^2-69*x+125)*(48*x^2+8*x+39)));
AuditRow(56, "2,2,2", R!(x), R!(10*x^5+8*x^4-5*x^3-3*x^2+x));
AuditRow(57, "2,2,4", R!((x^2+x)), R!(x^5+3*x^4-3*x^3-8*x^2+6*x));
AuditRow(58, "2,2,6", R!((x^2+1)), R!(3*x^5-4*x^3-x^2+x));
AuditRow(59, "2,2,8", R!((x^2+x)), R!(9*x^5+20*x^4-3*x^3-6*x^2-x));
AuditRow(60, "2,2,10", R!((x^2+1)), R!(-12*x^5+6*x^4+12*x^3-7*x^2));
AuditRow(61, "2,2,12", R!((x^2+x)), R!(x^6-7*x^5+8*x^4+16*x^3-5*x^2-5*x));
AuditRow(62, "2,2,14", R!(0), R!((x+1)*(x+9)*(2*x-115)*(6*x+55)*(3*x^2-220*x+2652)));
AuditRow(63, "2,2,20", R!(0), R!(-(x-1)*(6*x+1)*(2*x+7)*(6217*x+1008)*(21*x^2-161*x+144)));
AuditRow(64, "2,4,4", R!((x^2+x)), R!(-116*x^5+2204*x^4-1850*x^3-210*x^2+180*x));
AuditRow(65, "2,4,8", R!(0), R!(x*(3*x-5)*(5*x+333)*(16*x^2+23*x+576)));
AuditRow(66, "2,2,2,2", R!((x^2+x)), R!(-27*x^5+23*x^4+8*x^3-4*x^2-x));
AuditRow(67, "2,2,2,4", R!((x^2+x)), R!(-225*x^5+146*x^4-16*x^3-6*x^2+x));
AuditRow(68, "2,2,2,6", R!((x^2+1)), R!(-6*x^5-14*x^4+18*x^3+11*x^2-12*x+2));
AuditRow(69, "2,2,2,8", R!(0), R!(x*(x+1)*(x+55^2)*(x+99^2)*(x+125^2)));
AuditRow(70, "2,2,2,10", R!(0), R!(x*(x+1)*(x-1)*(3*x-7)*(8*x-13)*(24*x+25)));
AuditRow(71, "2,2,2,12", R!(0), R!(x*(x-120^2)*(x-143^2)*(x-266^2)*(x-218^2)*(x-241^2)));
AuditRow(72, "2,2,4,4", R!(0), R!(x*(x+36^2)*(x+57^2)*(x+64^2)*(x+132^2)));
AuditRow(73, "2,2,2,12#fiber1", R!(0), R!(x*(x-408^2)*(x-143^2)*(x-1015^2)*(x-437^2)*(x-1013^2)));
AuditRow(74, "2,2,2,12#fiber3", R!(0), R!(x*(x-16660^2)*(x-78793^2)*(x-21456^2)*(x-78644^2)*(x-27593^2)));
print "AUDIT_DONE";
quit;
