//////////////////////////////////////////////////////////////////////
//  Genuine 7-division cover above the contact-7 family.
//
//  The contact family is
//
//      h = 1 - (7/2)*x + alpha*x^2 + beta*x^3,
//      f = (h^2 + (x-1)^7)/x^2,
//      D7 = [(x-1), h(1)] in Jac(y^2=f).
//
//  A degree-2 Mumford class Q=[U,V], U monic quadratic, satisfies
//  7*Q = +/-D7 on the generic chart precisely when there are A,B with
//
//      A^2 - B^2*f = -(x-1)*U^7,                         (N49)
//
//  deg(A)<=7 and B monic of degree 5, subject to the usual open
//  conditions gcd(B,U)=gcd(U,f)=1.  Then V is recovered as
//
//      V = -A/B mod U,
//
//  and the sign is read at x=1.  The divisor of
//  (A+B*y)/(x-1), or its conjugate, proves the relation.
//
//  Modes:
//      mode:="finite"   exact D7 in 7*J(F_p) screens (default)
//      mode:="symbolic" report the compact generic norm chart
//      mode:="verify"   verify supplied norm parameters
//      mode:="exact"    exact rational checks on five small samples
//
//  Typical low-memory runs:
//      magma -b mode:="finite" primes:="3,5,11,13" \
//          code/contact7_to49_division.m
//      magma -b mode:="symbolic" code/contact7_to49_division.m
//
//  In verify mode assign the nine rational parameters
//      alpha,beta,s,t,c0,c1,c2,c3,r
//  where U=x^2+s*x+t, r=lc(A), and c0..c3 are the low
//  coefficients of B.  Its x^4 coefficient and the other coefficients
//  of A are forced by the high norm equations.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then
    mode := "finite";
end if;
if not assigned primes then
    primes := "3,5,11,13";
end if;
if not assigned max_print then
    max_print := 12;
elif Type(max_print) eq MonStgElt then
    max_print := StringToInteger(max_print);
end if;

function ParseIntegerList(S)
    if Type(S) eq SeqEnum then
        return [Integers()!z : z in S];
    end if;
    return [StringToInteger(t) : t in Split(S, ",")];
end function;

function Contact7Data(k, av, bv)
    P<x> := PolynomialRing(k);
    h := 1 - (k!7/k!2)*x + av*x^2 + bv*x^3;
    num := h^2 + (x-1)^7;
    if Coefficient(num,0) ne 0 or Coefficient(num,1) ne 0 then
        return false, P!0, P!0;
    end if;
    f := ExactQuotient(num,x^2);
    return true,f,h;
end function;

// Solve n*q=a in each invariant-factor coordinate.  Only existence is
// needed for the local obstruction.  This is exact, not merely #J mod 49.
function IsDivisibleInFiniteGroup(D, n)
    J := Parent(D);
    G,phi := AbelianGroup(J);
    g := D @@ phi;
    coords := Eltseq(g);
    invs := Invariants(G);
    for i in [1..#coords] do
        if (Integers()!coords[i]) mod GCD(n,Integers()!invs[i]) ne 0 then
            return false,invs,coords;
        end if;
    end for;
    return true,invs,coords;
end function;

// Generic a7 != 0 chart for (N49).  The x^15 equation is automatic.
// The x^14 equation forces b4, and x^13,...,x^7 successively force
// a6,...,a0.  Seven residual equations, in degrees 0,...,6, remain.
function Norm49Data(k,av,bv,sv,tv,c0v,c1v,c2v,c3v,rv)
    P<x> := PolynomialRing(k);
    if rv eq 0 then
        return false,P!0,P!0,P!0,P!0,[],"r=0 chart";
    end if;
    ok,f,h := Contact7Data(k,av,bv);
    if not ok then
        return false,P!0,P!0,P!0,P!0,[],"bad contact formula";
    end if;

    U := x^2 + sv*x + tv;
    // coeff x^14 of A^2-B^2*f+(x-1)U^7 is
    // r^2-beta^2-2*b4+7*s+6.
    b4v := (rv^2-bv^2+(k!7)*sv+k!6)/(k!2);
    B := x^5+b4v*x^4+c3v*x^3+c2v*x^2+c1v*x+c0v;
    A := rv*x^7;
    E := A^2-B^2*f+(x-1)*U^7;
    if Coefficient(E,15) ne 0 or Coefficient(E,14) ne 0 then
        return false,A,B,f,U,[],"high-equation failure";
    end if;

    for d in [13..7 by -1] do
        ad := -Coefficient(E,d)/(k!2*rv);
        A +:= ad*x^(d-7);
        E := A^2-B^2*f+(x-1)*U^7;
    end for;
    residuals := [Coefficient(E,i) : i in [0..6]];
    return &and[z eq 0 : z in residuals],A,B,f,U,residuals,"ok";
end function;

procedure FiniteScreen()
    print "# contact-7 to cyclic-49 exact finite-group screen";
    print "# DIVISIBLE means the marked D7 lies in 7*J(F_p), not just 49 | #J(F_p).";
    for p in ParseIntegerList(primes) do
        if not IsPrime(p) or p in {2,7} then
            error "primes must be odd and different from 7";
        end if;
        k := GF(p);
        P<x> := PolynomialRing(k);
        checked := 0;
        smooth := 0;
        marked7 := 0;
        order49_possible := 0;
        divisible := 0;
        shown := 0;
        allowed := [];

        for av in k do
            for bv in k do
                checked +:= 1;
                ok,f,h := Contact7Data(k,av,bv);
                if not ok or Degree(f) ne 5 or Discriminant(f) eq 0 then
                    continue;
                end if;
                yp := Evaluate(h,k!1);
                if yp eq 0 then
                    continue;
                end if;
                smooth +:= 1;
                J := Jacobian(HyperellipticCurve(f));
                D7 := J![x-1,yp];
                if Order(D7) ne 7 then
                    continue;
                end if;
                marked7 +:= 1;
                if (#J mod 49) eq 0 then
                    order49_possible +:= 1;
                end if;
                isdiv,invs,coords := IsDivisibleInFiniteGroup(D7,7);
                if not isdiv then
                    continue;
                end if;
                divisible +:= 1;
                Append(~allowed,<Integers()!av,Integers()!bv>);
                if shown lt max_print then
                    printf "DIVISIBLE p=%o alpha=%o beta=%o inv=%o coordsD7=%o #J=%o\n",
                        p,Integers()!av,Integers()!bv,invs,coords,#J;
                    shown +:= 1;
                end if;
            end for;
        end for;
        printf "SUMMARY p=%o checked=%o smooth=%o marked7=%o order49_possible=%o D7_divisible_by_7=%o\n",
            p,checked,smooth,marked7,order49_possible,divisible;
        print "ALLOWED_PAIRS",p,allowed;
    end for;
end procedure;

procedure SymbolicReport()
    Q := Rationals();
    R<aa,bb,ss,tt,a0,a1,a2,a3,a4,a5,a6,a7,b0,b1,b2,b3,b4> :=
        PolynomialRing(Q,17,"grevlex");
    P<x> := PolynomialRing(R);
    h := 1-(R!7/R!2)*x+aa*x^2+bb*x^3;
    f := ExactQuotient(h^2+(x-1)^7,x^2);
    U := x^2+ss*x+tt;
    A := a7*x^7+a6*x^6+a5*x^5+a4*x^4+a3*x^3+a2*x^2+a1*x+a0;
    B := x^5+b4*x^4+b3*x^3+b2*x^2+b1*x+b0;
    E := A^2-B^2*f+(x-1)*U^7;
    assert Degree(E) le 14;
    print "# compact uneliminated contact-7 to 49 norm cover";
    print "# variables (17): alpha,beta,s,t,a0..a7,b0..b4; B is monic";
    print "# equations (15): coeff_x^i(E)=0 for i=0..14";
    print "# E=A^2-B^2*f+(x-1)U^7, U=x^2+s*x+t";
    print "# coefficient 15 cancels identically; expected dimension is 2";
    print "TOP14",Coefficient(E,14);
    for i in [0..14] do
        ei := Coefficient(E,i);
        printf "E%o total_degree=%o terms=%o\n",i,TotalDegree(ei),#Terms(ei);
    end for;
end procedure;

procedure VerifyNormPoint()
    if not (assigned alpha and assigned beta and assigned s and assigned t and
            assigned c0 and assigned c1 and assigned c2 and assigned c3 and
            assigned r) then
        error "verify mode needs alpha,beta,s,t,c0,c1,c2,c3,r";
    end if;
    Q := Rationals();
    vals := [Q!alpha,Q!beta,Q!s,Q!t,Q!c0,Q!c1,Q!c2,Q!c3,Q!r];
    ok,A,B,f,U,res,msg := Norm49Data(Q,vals[1],vals[2],vals[3],vals[4],
        vals[5],vals[6],vals[7],vals[8],vals[9]);
    print "NORM_OK",ok,"message",msg;
    print "f",f;
    print "U",U;
    print "A",A;
    print "B",B;
    print "residuals",res;
    if ok and Discriminant(f) ne 0 and GCD(U,f) eq 1 and
            GCD(U,B) eq 1 then
        V := (-A*InverseMod(B,U)) mod U;
        print "V",V,"V2_minus_f_mod_U",(V^2-f) mod U;
        L := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
        fI := L^2*f;
        J := Jacobian(HyperellipticCurve(fI));
        h := 1-(Q!7/Q!2)*Parent(f).1+vals[1]*Parent(f).1^2+vals[2]*Parent(f).1^3;
        D7 := J![Parent(f).1-1,L*Evaluate(h,Q!1)];
        QQ := J![U,L*V];
        print "orders",Order(D7),Order(QQ),"sevenQ",7*QQ,"D7",D7,
              "relation_plus",7*QQ eq D7,"relation_minus",7*QQ eq -D7;
    end if;
end procedure;

procedure ExactSamples()
    Q := Rationals();
    P<x> := PolynomialRing(Q);
    samples := [<Q!0,Q!0>,<Q!1,Q!0>,<Q!0,Q!1>,<Q!2,Q!-1>,<Q!-3,Q!2>];
    print "# exact IsDivisibleBy(D7,7) sample probe";
    for ab in samples do
        ok,f,h := Contact7Data(Q,ab[1],ab[2]);
        if not ok or Degree(f) ne 5 or Discriminant(f) eq 0 or Evaluate(h,Q!1) eq 0 then
            print "EXACT_SAMPLE",ab,"bad";
            continue;
        end if;
        L := LCM([Denominator(Coefficient(f,i)) : i in [0..Degree(f)]]);
        fI := L^2*f;
        J := Jacobian(HyperellipticCurve(fI));
        D7 := J![x-1,L*Evaluate(h,Q!1)];
        try
            isdiv,Q49 := IsDivisibleBy(D7,7);
            print "EXACT_SAMPLE",ab,"D7order",Order(D7),"div7",isdiv;
            if isdiv then
                print "  Q",Q49,"Qorder",Order(Q49),"check",7*Q49 eq D7;
            end if;
        catch e
            print "EXACT_SAMPLE",ab,"D7order",Order(D7),"ERROR",e`Object;
        end try;
    end for;
end procedure;

if mode eq "finite" then
    FiniteScreen();
elif mode eq "symbolic" then
    SymbolicReport();
elif mode eq "exact" then
    ExactSamples();
elif mode eq "verify" then
    VerifyNormPoint();
else
    error "unknown mode";
end if;

quit;
