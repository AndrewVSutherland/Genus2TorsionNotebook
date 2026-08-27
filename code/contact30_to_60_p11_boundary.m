//////////////////////////////////////////////////////////////////////
//  The p=11 boundary for halving the contact-(5,6) order-30 class.
//
//  It classifies P^1(F_11), prints the factors at which the reduced
//  rational parametrisation ceases to be a good quintic, and tests
//  2-divisibility of D_6 (equivalently D_30) in every good finite
//  Jacobian.  The blow-up charts also certify the local obstructions
//  in every boundary disk.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

F := GF(11);
PR<R> := PolynomialRing(F);
K := FieldOfFractions(PR);
PK<x> := PolynomialRing(K);

function FacString(z)
    if z eq 0 then
        return "0";
    end if;
    return Sprint(Factorization(PR!z));
end function;

function AtInfinity(z)
    if z eq 0 then
        return true, F!0;
    end if;
    n := Numerator(z);
    d := Denominator(z);
    if Degree(n) gt Degree(d) then
        return false, F!0;
    elif Degree(n) lt Degree(d) then
        return true, F!0;
    end if;
    return true, LeadingCoefficient(n)/LeadingCoefficient(d);
end function;

function IsDoubleFinite(J,D)
    G, phi := AbelianGroup(J);
    a := D @@ phi;
    H := sub<G | [2*G.i : i in [1..Ngens(G)]]>;
    return a in H, Invariants(G), Eltseq(a);
end function;

// Coefficient comparison for
//   (T^3+pT^2+qT+r)^2 - (T^2+sT+z) f(T)
// with f(T)=h(T)^2-T^6.  Eliminating z and the upper coefficients
// gives these four equations for a half of [T,h(0)].
function HalfCoverData(A,B,C)
    sols := [];
    smooth := 0;
    ranks := [0,0,0,0,0];
    for p,q,r,s in F do
        H1 := A*q^2 - C*(B+s);
        H2 := -2*A*C + 4*A*p*q - B^2 - 2*C*r + s^2;
        H3 := -2*A*B + 2*A*p^2 + 4*A*q - 2*C + 2*r*s;
        H4 := -A^2 + 4*A*p - 2*B + r^2;
        if H1 eq 0 and H2 eq 0 and H3 eq 0 and H4 eq 0 then
            M := Matrix(F,4,4,[
                0,                 2*A*q,        0,       -C,
                4*A*q,             4*A*p,       -2*C,    2*s,
                4*A*p,             4*A,         2*s,     2*r,
                4*A,               0,           2*r,     0
            ]);
            rk := Rank(M);
            ranks[rk+1] +:= 1;
            if rk eq 4 then smooth +:= 1; end if;
            if #sols lt 8 then Append(~sols,<p,q,r,s,rk>); end if;
        end if;
    end for;
    return &+ranks,smooth,ranks,sols;
end function;

function ValAtLinear(z,r)
    if z eq 0 then return 100; end if;
    n := Numerator(z);
    d := Denominator(z);
    v := 0;
    while Evaluate(n,r) eq 0 do
        n := ExactQuotient(n,R-r);
        v +:= 1;
    end while;
    while Evaluate(d,r) eq 0 do
        d := ExactQuotient(d,R-r);
        v -:= 1;
    end while;
    return v;
end function;

function FamilyK(eps)
    denR := R^2 - 5;
    t := (5*R^2 - 20*R + 19)/denR;
    Y := -2*(5*R^2 - 22*R + 25)/denR;
    u := t^3;
    s := t^5 + t^4 + (F!5/F!2)*t^3 + (F!1/F!2)*t
          + eps*t*(t-F!1/F!2)*(t+1)*Y;
    C := (u^2+1)/(2*u);
    c := (u^2-1)/(2*u);
    denq := u^6 + 6*u^4*s - 2*u^4 + 15*u^3*s
            - u*s^3 + u^2;
    numq := 15*u^5 + 90*u^4 + 20*u^3*s - 6*u^2*s^2
            + 231*u^3 + 2*u^2*s - 15*u*s^2 + 90*u^2
            - 20*u*s + 15*u - 2*s;
    q := numq/denq;
    A := (s+q)/2;
    e := (s-q)/2;
    B := (15-s*q)/2;
    d := (B*C+3)/c;
    h6 := x^3 + A*x^2 + B*x + C;
    h5 := e*x^2 + d*x + c;
    f := h6^2 - (x-1)^6;
    return f,h6,h5,t,Y,u,s,q,C,c,A,B,d,denq,numq;
end function;

for eps in [F!1,-F!1] do
    f,h6,h5,t,Y,u,s,q,C,c,A,B,d,denq,numq := FamilyK(eps);
    printf "\n=== eps=%o rational factors over F_11 ===\n", Integers()!eps;
    printf "t.num=%o\nt.den=%o\n", FacString(Numerator(t)), FacString(Denominator(t));
    printf "u.num=%o\nu.den=%o\n", FacString(Numerator(u)), FacString(Denominator(u));
    printf "c.num=%o\nc.den=%o\n", FacString(Numerator(c)), FacString(Denominator(c));
    printf "q.num=%o\nq.den=%o\n", FacString(Numerator(q)), FacString(Denominator(q));
    printf "A.den=%o\nB.den=%o\nC.den=%o\nd.den=%o\n",
           FacString(Denominator(A)), FacString(Denominator(B)),
           FacString(Denominator(C)), FacString(Denominator(d));
    coeffs := [Coefficient(f,i) : i in [0..5]];
    discf := Discriminant(f);
    for i in [0..5] do
        printf "f%o.num=%o\nf%o.den=%o\n", i,
               FacString(Numerator(coeffs[i+1])), i,
               FacString(Denominator(coeffs[i+1]));
    end for;

    printf "\nprojective residues eps=%o\n", Integers()!eps;
    for rr in [0..10] do
        r := F!rr;
        printf "R=%o valuations f[0..5]=%o disc=%o ", rr,
               [ValAtLinear(z,r) : z in coeffs], ValAtLinear(discf,r);
        vals := [];
        defined := true;
        for z in coeffs cat [Evaluate(h6,F!1)] do
            if Evaluate(Denominator(z),r) eq 0 then
                defined := false;
                break;
            end if;
            Append(~vals,Evaluate(z,r));
        end for;
        if not defined then
            print "BOUNDARY coefficient-pole";
            continue;
        end if;
        Pf<xx> := PolynomialRing(F);
        ff := &+[vals[i+1]*xx^i : i in [0..5]];
        if Degree(ff) ne 5 then
            printf "BOUNDARY degree=%o f=%o\n", Degree(ff), ff;
            Ac := Evaluate(A+3,r);
            Bc := Evaluate(B+2*A+3,r);
            Cc := Evaluate(1+A+B+C,r);
            nsol,nsm,ranks,sols := HalfCoverData(Ac,Bc,Cc);
            printf "  t=%o Y=%o u=%o s=%o q=%o centered_h=(%o,%o,%o) halfcover n=%o smooth=%o ranks=%o sample=%o\n",
                   Evaluate(t,r),Evaluate(Y,r),Evaluate(u,r),Evaluate(s,r),
                   Evaluate(q,r),Ac,Bc,Cc,nsol,nsm,ranks,sols;
            continue;
        elif Discriminant(ff) eq 0 then
            printf "BOUNDARY singular factorization=%o\n", Factorization(ff);
            Ac := Evaluate(A+3,r);
            Bc := Evaluate(B+2*A+3,r);
            Cc := Evaluate(1+A+B+C,r);
            nsol,nsm,ranks,sols := HalfCoverData(Ac,Bc,Cc);
            printf "  t=%o Y=%o u=%o s=%o q=%o centered_h=(%o,%o,%o) halfcover n=%o smooth=%o ranks=%o sample=%o\n",
                   Evaluate(t,r),Evaluate(Y,r),Evaluate(u,r),Evaluate(s,r),
                   Evaluate(q,r),Ac,Bc,Cc,nsol,nsm,ranks,sols;
            continue;
        end if;
        hh6 := &+[Evaluate(Coefficient(h6,i),r)*xx^i : i in [0..3]];
        J := Jacobian(HyperellipticCurve(ff));
        D6 := J![xx-1,Evaluate(hh6,F!1)];
        isdbl,invs,coords := IsDoubleFinite(J,D6);
        Ac := Evaluate(A+3,r);
        Bc := Evaluate(B+2*A+3,r);
        Cc := Evaluate(1+A+B+C,r);
        nsol,nsm,ranks,sols := HalfCoverData(Ac,Bc,Cc);
        printf "GOOD #J=%o inv=%o D6ord=%o D6coords=%o double=%o halfcover=%o\n",
               #J,invs,Order(D6),coords,isdbl,nsol;
    end for;

    valsinf := [];
    definedinf := true;
    for z in coeffs cat [Evaluate(h6,F!1)] do
        ok,v := AtInfinity(z);
        if not ok then
            definedinf := false;
            break;
        end if;
        Append(~valsinf,v);
    end for;
    if not definedinf then
        print "R=inf BOUNDARY coefficient-pole";
    else
        Pf<xx> := PolynomialRing(F);
        ff := &+[valsinf[i+1]*xx^i : i in [0..5]];
        if Degree(ff) ne 5 then
            printf "R=inf BOUNDARY degree=%o f=%o\n", Degree(ff),ff;
        elif Discriminant(ff) eq 0 then
            printf "R=inf BOUNDARY singular factorization=%o\n", Factorization(ff);
        else
            okh,h1 := AtInfinity(Evaluate(h6,F!1));
            hcoeffinf := [];
            for i in [0..3] do
                oki,vi := AtInfinity(Coefficient(h6,i));
                Append(~hcoeffinf,vi);
            end for;
            hh6 := &+[hcoeffinf[i+1]*xx^i : i in [0..3]];
            J := Jacobian(HyperellipticCurve(ff));
            D6 := J![xx-1,h1];
            isdbl,invs,coords := IsDoubleFinite(J,D6);
            printf "R=inf GOOD #J=%o inv=%o D6ord=%o D6coords=%o double=%o\n",
                   #J,invs,Order(D6),coords,isdbl;
        end if;
    end if;
end for;

//////////////////////////////////////////////////////////////////////
// Generic points of the first blow-up R=r+11*z over Q_11.
// For each rational function g(z), print its Gauss valuation and the
// residue of 11^{-v(g)}g on the exceptional P^1_z/F_11.
//////////////////////////////////////////////////////////////////////

QQ := Rationals();
PZ<z> := PolynomialRing(QQ);
L := FieldOfFractions(PZ);
PL<X> := PolynomialRing(L);
FZ<zb> := PolynomialRing(F);
LF := FieldOfFractions(FZ);
LFU<U> := PolynomialRing(LF);

function V11Q(a)
    if a eq 0 then return 1000000; end if;
    return Valuation(Numerator(a),11)-Valuation(Denominator(a),11);
end function;

function ReduceIntegralPolynomial(g)
    ans := FZ!0;
    for i in [0..Degree(g)] do
        a := Coefficient(g,i);
        if a ne 0 then
            ans +:= (F!(Numerator(a) mod 11)/F!(Denominator(a) mod 11))*zb^i;
        end if;
    end for;
    return ans;
end function;

function GaussLead(g)
    if g eq 0 then return 1000000,LF!0; end if;
    n := Numerator(g);
    d0 := Denominator(g);
    vn := Minimum([V11Q(a) : a in Coefficients(n) | a ne 0]);
    vd := Minimum([V11Q(a) : a in Coefficients(d0) | a ne 0]);
    nn := n/(QQ!11)^vn;
    dd := d0/(QQ!11)^vd;
    nb := ReduceIntegralPolynomial(nn);
    db := ReduceIntegralPolynomial(dd);
    return vn-vd,LF!nb/LF!db;
end function;

function FamilyQ(parameter,eps)
    denR := parameter^2 - 5;
    t := (5*parameter^2 - 20*parameter + 19)/denR;
    Y := -2*(5*parameter^2 - 22*parameter + 25)/denR;
    u := t^3;
    s := t^5 + t^4 + (QQ!5/QQ!2)*t^3 + (QQ!1/QQ!2)*t
          + eps*t*(t-QQ!1/QQ!2)*(t+1)*Y;
    C := (u^2+1)/(2*u);
    c := (u^2-1)/(2*u);
    denq := u^6 + 6*u^4*s - 2*u^4 + 15*u^3*s
            - u*s^3 + u^2;
    numq := 15*u^5 + 90*u^4 + 20*u^3*s - 6*u^2*s^2
            + 231*u^3 + 2*u^2*s - 15*u*s^2 + 90*u^2
            - 20*u*s + 15*u - 2*s;
    q := numq/denq;
    A := (s+q)/2;
    B := (15-s*q)/2;
    h6 := X^3 + A*X^2 + B*X + C;
    f := h6^2-(X-1)^6;
    Ac := A+3;
    Bc := B+2*A+3;
    Cc := 1+A+B+C;
    return f,Ac,Bc,Cc,t,Y,u,s,q;
end function;

print "\n=== generic first blow-up charts over Q_11 ===";
print "\n=== exact norm obstructions ===";
for epsi in [1,-1] do
    fq,Acq,Bcq,Ccq,tq,Yq,uq,sq,qq := FamilyQ(L!z,QQ!epsi);
    Nm := Ccq/Acq;
    Np := -Ccq/2;
    printf "eps=%o Nminus=Ccenter/Acenter\n  num=%o\n  den=%o\n",
           epsi,Factorization(Numerator(Nm)),Factorization(Denominator(Nm));
    printf "eps=%o Nplus=-Ccenter/2\n  num=%o\n  den=%o\n",
           epsi,Factorization(Numerator(Np)),Factorization(Denominator(Np));
end for;

print "\n=== generic first blow-up charts over Q_11 ===";
for epsi in [1,-1] do
    printf "\neps=%o\n",epsi;
    for rr in [1,2,3,4,5,6,7,10] do
        parameter := L!(rr+11*z);
        fq,Acq,Bcq,Ccq,tq,Yq,uq,sq,qq := FamilyQ(parameter,QQ!epsi);
        va,la := GaussLead(Acq);
        vb,lb := GaussLead(Bcq);
        vc,lc := GaussLead(Ccq);
        vf := [];
        lf := [];
        for i in [0..5] do
            vi,li := GaussLead(Coefficient(fq,i));
            Append(~vf,vi); Append(~lf,li);
        end for;
        vd,ld := GaussLead(Discriminant(fq));
        vm,lm := GaussLead(Ccq/Acq);
        vp,lp := GaussLead(-Ccq/2);
        printf "R=%o+11z centered v=(%o,%o,%o) lead=(%o,%o,%o)\n",
               rr,va,vb,vc,la,lb,lc;
        printf "  norms Nminus=(v=%o,lead=%o) Nplus=(v=%o,lead=%o)\n",
               vm,lm,vp,lp;
        if rr in {4,5,7,10} then
            qminus_lead := la*U^2+lb*U+lc;
            // The normalized quadratic etale factor always specializes
            // to lambda*(U+1)^2.  Hence its quadratic algebra is split or
            // ramified with residue field F_11.  The image of U is -1,
            // a nonsquare, so U cannot be a square in that factor.
            assert Degree(qminus_lead) eq 2;
            assert Discriminant(qminus_lead) eq 0;
            assert Evaluate(qminus_lead,-LF!1) eq 0;
            printf "  quadratic Kummer lead=%o = lambda*(U+1)^2\n",
                   qminus_lead;
        end if;
        printf "  fvals=%o fleads=%o disc=(v=%o, lead=%o)\n",vf,lf,vd,ld;
    end for;
end for;

print "\n=== second charts inside the four pole disks ===";
pole_data := [<4,4>,<7,6>,<5,5>,<10,5>];
for epsi in [1,-1] do
    printf "eps=%o\n",epsi;
    for rz in pole_data do
        rr := rz[1]; z0 := rz[2];
        parameter := L!(rr+11*z0+121*z);
        fq,Acq,Bcq,Ccq,tq,Yq,uq,sq,qq := FamilyQ(parameter,QQ!epsi);
        vm,lm := GaussLead(Ccq/Acq);
        vp,lp := GaussLead(-Ccq/2);
        va,la := GaussLead(Acq);
        vb,lb := GaussLead(Bcq);
        vc,lc := GaussLead(Ccq);
        qminus_lead := la*U^2+lb*U+lc;
        assert Degree(qminus_lead) eq 2;
        assert Discriminant(qminus_lead) eq 0;
        assert Evaluate(qminus_lead,-LF!1) eq 0;
        printf "R=%o+11*%o+121w Nminus=(v=%o,lead=%o) Nplus=(v=%o,lead=%o)\n",
               rr,z0,vm,lm,vp,lp;
        printf "  quadratic Kummer lead=%o = lambda*(U+1)^2\n",
               qminus_lead;
        vals := [];
        for ww in [0..10] do
            dn := Evaluate(Denominator(lp),F!ww);
            if dn eq 0 then
                Append(~vals,Sprintf("%o:deeper",ww));
            else
                vv := Evaluate(Numerator(lp),F!ww)/dn;
                Append(~vals,Sprintf("%o:%o:%o",ww,Integers()!vv,IsSquare(vv)));
            end if;
        end for;
        printf "  Nplus exceptional-line values=%o\n",vals;
    end for;
end for;

quit;
