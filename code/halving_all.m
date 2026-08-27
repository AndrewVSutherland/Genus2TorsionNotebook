//Written by GPT 5.5 Thinking; last part checks it on all tuples in tor2244.txt;
//all sections pass and in particular Section 3 passed
//this verifies the halving formula on all 26653 tuples from tor2244.txt
//////////////////////////////////////////////////////////////////////
//  1. Generic check on the transformed curve
//
//  Curve:  Y^2 = z(z+1)(z+rho^2)(z+sigma^2)(z+tau^2)
//  Claim:  [ z^2 - s2*z + s4 , (s1*s2 - s3)*z - s1*s4 ]
//          doubles to [z,0], i.e. to (0,0)-infinity.
//////////////////////////////////////////////////////////////////////

Q := Rationals();
K<rho,sigma,tau> := RationalFunctionField(Q, 3);
Pz<z> := PolynomialRing(K);

s1 := 1 + rho + sigma + tau;
s2 := rho + sigma + tau + rho*sigma + rho*tau + sigma*tau;
s3 := rho*sigma + rho*tau + sigma*tau + rho*sigma*tau;
s4 := rho*sigma*tau;

ftr := z*(z+1)*(z+rho^2)*(z+sigma^2)*(z+tau^2);
Ctr := HyperellipticCurve(ftr);
Jtr := Jacobian(Ctr);

Utr := z^2 - s2*z + s4;
Vtr := (s1*s2 - s3)*z - s1*s4;

H := Jtr![Utr, Vtr];
T0 := Jtr![z, 0];

assert 2*H eq T0;
print "Section 1 passed: generic halving formula on transformed curve verified.";


//////////////////////////////////////////////////////////////////////
//  2. Generic symbolic check of the pullback to the original x-model
//
//  We verify:
//
//    (a) the birational change of variables
//    (b) the pulled-back U(x)
//    (c) the explicit expanded V(x)
//
//  This part does NOT use Jacobian arithmetic; it is just algebra.
//////////////////////////////////////////////////////////////////////

K2<A,B,C,D,bv,s1g,s2g,s3g,s4g> := RationalFunctionField(Q, 9);

// ----- (a) Check the curve transformation algebraically -----

L<zz,YY> := RationalFunctionField(K2, 2);

xx := -(A + B*zz)/(1 + zz);
yy := bv*(B - A)/(1 + zz)^3 * YY;

// cleared numerator of y^2 - x(x+A)(x+B)(x+C)(x+D)
num := Numerator(yy^2 - xx*(xx+A)*(xx+B)*(xx+C)*(xx+D));

// expected cleared identity
target :=
-(A - B)^2 *
(
    bv^2*YY^2
    - zz*(1+zz)*(A + B*zz)*(A - C + (B - C)*zz)*(A - D + (B - D)*zz)
);

assert num eq target or num eq -target;
print "Section 2(a) passed: birational change of variables checked.";


// ----- (b) Check the pulled-back U(x) -----

Px<x> := PolynomialRing(K2);

delta := 1 + s2g + s4g;

U_from_pullback :=
((x + A)^2 + s2g*(x + A)*(x + B) + s4g*(x + B)^2) / delta;

U_explicit :=
x^2
+ ((2*A + s2g*(A+B) + 2*B*s4g)/delta)*x
+ (A^2 + A*B*s2g + B^2*s4g)/delta;

assert U_from_pullback eq U_explicit;
print "Section 2(b) passed: explicit U(x) checked.";


// ----- (c) Check the expanded V(x) -----

m := s1g*s2g - s3g;

// pulled-back cubic before reduction mod U
Vpull :=
(bv/(B-A)^2) * (x+B)^2 * ( -m*(x+A) - s1g*s4g*(x+B) );

// divide by U_explicit and take remainder
q, Vrem := Quotrem(Vpull, U_explicit);

Lambda :=
    s1g*s2g - s1g*s4g^2 + 3*s1g*s4g + s2g*s3g*s4g + 3*s3g*s4g - s3g;

M :=
    A*s1g*s2g + 2*A*s1g*s4g + A*s3g*s4g - A*s3g
    - B*s1g*s4g^2 + B*s1g*s4g + B*s2g*s3g*s4g + 2*B*s3g*s4g;

Vexplicit := -(bv/delta^2) * (Lambda*x + M);

assert Vrem eq Vexplicit;
print "Section 2(c) passed: explicit V(x) checked.";


//////////////////////////////////////////////////////////////////////
//  3. Concrete checks on every tuple in tor2244.txt
//
//  The file tor2244.txt should have one tuple per line, e.g.
//      [a,b,c,d]
//  Any non-tuple lines, such as the final status line, are ignored.
//
//  For each tuple we check directly in the Jacobian of
//      y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2)
//  that the proposed Mumford divisor D' satisfies
//      2 D' = (P_a - infinity) + (P_b - infinity).
//////////////////////////////////////////////////////////////////////

Qx := Rationals();
Px0<x0> := PolynomialRing(Qx);

tuple_file := "tor2244.txt";
progress_interval := 100;

function ReadTupleFile(filename)
    S := Read(filename);
    rows := Split(S, "\n");
    tuples := [];

    for rawrow in rows do
        // Handle either Unix or Windows line endings without assigning
        // to the loop variable, which Magma disallows.
        n := #rawrow;
        if n gt 0 then
            if rawrow[n] eq "\r" then
                n := n - 1;
            end if;
        end if;

        // tor2244.txt has bracketed tuple rows and a final non-tuple
        // status row.  We only parse rows of the form [a,b,c,d].
        if n ge 2 then
            if rawrow[1] eq "[" then
                if rawrow[n] eq "]" then
                    body := rawrow[2..n-1];
                    parts := Split(body, ",");
                    tup := [ StringToInteger(part) : part in parts ];
                    if #tup eq 4 then
                        Append(~tuples, tup);
                    end if;
                end if;
            end if;
        end if;
    end for;

    return tuples;
end function;

function SquareRootOrFail(n, label, idx, tup)
    ok, r := IsSquare(n);
    if not ok then
        print "FAILED: expected a square.";
        print "tuple index:", idx;
        print "tuple:", tup;
        print "quantity:", label;
        print "value:", n;
        error "Stopping because one of the square conditions failed.";
    end if;
    return r;
end function;

procedure VerifyTuple(tup, idx)
    a := Qx!tup[1];
    b := Qx!tup[2];
    c := Qx!tup[3];
    d := Qx!tup[4];

    A0 := a^2;
    B0 := b^2;
    C0 := c^2;
    D0 := d^2;

    if A0 eq B0 or A0 eq C0 or A0 eq D0 or B0 eq C0 or B0 eq D0 or C0 eq D0 then
        print "FAILED: repeated branch point.";
        print "tuple index:", idx;
        print "tuple:", tup;
        error "Stopping because the curve is singular.";
    end if;

    u0 := SquareRootOrFail((A0-C0)*(A0-D0), "u^2", idx, tup);
    v0 := SquareRootOrFail((B0-C0)*(B0-D0), "v^2", idx, tup);
    w0 := SquareRootOrFail((A0-C0)*(B0-C0), "w^2", idx, tup);
    t0 := SquareRootOrFail((A0-D0)*(B0-D0), "t^2", idx, tup);

    // Adjust the sign of u0 if necessary so that uv = wt.
    // This sign is not used in the final Mumford formula below, but it
    // records the convention from the statement being checked.
    if u0*v0 ne w0*t0 then
        u0 := -u0;
    end if;

    if u0*v0 ne w0*t0 then
        print "FAILED: could not choose signs with uv = wt by changing u.";
        print "tuple index:", idx;
        print "tuple:", tup;
        print "u, v, w, t:", u0, v0, w0, t0;
        error "Stopping because the sign convention uv = wt failed.";
    end if;

    rho0 := a/b;
    sigma0 := (A0-C0)/w0;   // = w0/(B0-C0), up to the same sign choice
    tau0 := (A0-D0)/t0;     // = t0/(B0-D0), up to the same sign choice

    s10 := 1 + rho0 + sigma0 + tau0;
    s20 := rho0 + sigma0 + tau0 + rho0*sigma0 + rho0*tau0 + sigma0*tau0;
    s30 := rho0*sigma0 + rho0*tau0 + sigma0*tau0 + rho0*sigma0*tau0;
    s40 := rho0*sigma0*tau0;
    delta0 := 1 + s20 + s40;

    if delta0 eq 0 then
        print "FAILED: delta = 0, so the displayed monic U formula is undefined.";
        print "tuple index:", idx;
        print "tuple:", tup;
        error "Stopping because delta vanished.";
    end if;

    f0 := x0*(x0+A0)*(x0+B0)*(x0+C0)*(x0+D0);
    C0curve := HyperellipticCurve(f0);
    J0 := Jacobian(C0curve);

    U0 :=
    x0^2
    + ((2*A0 + s20*(A0+B0) + 2*B0*s40)/delta0)*x0
    + (A0^2 + A0*B0*s20 + B0^2*s40)/delta0;

    Lambda0 :=
        s10*s20 - s10*s40^2 + 3*s10*s40 + s20*s30*s40 + 3*s30*s40 - s30;

    M0 :=
        A0*s10*s20 + 2*A0*s10*s40 + A0*s30*s40 - A0*s30
        - B0*s10*s40^2 + B0*s10*s40 + B0*s20*s30*s40 + 2*B0*s30*s40;

    V0 := -(b*v0/delta0^2) * (Lambda0*x0 + M0);

    Dprime := J0![U0, V0];

    // target class = (P_a - infinity) + (P_b - infinity)
    Dtarget := J0![x0 + A0, 0] + J0![x0 + B0, 0];

    if 2*Dprime ne Dtarget then
        print "FAILED: Jacobian equality 2*Dprime = Dtarget did not hold.";
        print "tuple index:", idx;
        print "tuple:", tup;
        print "u, v, w, t:", u0, v0, w0, t0;
        print "rho, sigma, tau:", rho0, sigma0, tau0;
        print "U0:", U0;
        print "V0:", V0;
        error "Stopping because the halving check failed.";
    end if;
end procedure;

tuples := ReadTupleFile(tuple_file);
print "Read", #tuples, "tuples from", tuple_file;

for i in [1..#tuples] do
    VerifyTuple(tuples[i], i);

    if i mod progress_interval eq 0 then
        print "Verified", i, "of", #tuples, "tuples.";
    end if;
end for;

print "Section 3 passed: verified the halving formula on all", #tuples, "tuples from", tuple_file;
