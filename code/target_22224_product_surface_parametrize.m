//////////////////////////////////////////////////////////////////////
// Parametrize the genus-zero rational-3-torsion cover on the minus
// elliptic quotient of the product-surface fibers t=1/2 and t=2/3.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();

for family in [5,6] do
    AA := family eq 5 select 16 else 144;
    BB := family eq 5 select 9 else 25;
    Kt<r> := RationalFunctionField(Q);
    PU<u> := PolynomialRing(Kt);
    s := (BB-AA*r)/(AA-BB*r); k := r*s;
    A := -2*k; B := 1+k^2; C := r^2+s^2;
    a2 := A+B+C; a4 := A*B+A*C+B*C; a6 := A*B*C;
    psi := 3*u^4+4*a2*u^3+6*a4*u^2+12*a6*u+(4*a2*a6-a4^2);
    FF<ux> := FunctionField(psi);
    cubic := (ux+FF!A)*(ux+FF!B)*(ux+FF!C);
    PY<Y> := PolynomialRing(FF);
    EE<vy> := ext<FF|Y^2-cubic>;
    assert Genus(FF) eq 0 and Genus(EE) eq 0;

    candidates := [];
    try
        candidates cat:= [P:P in Zeros(EE!r)|Degree(P) eq 1];
    catch e
        dummy := 0;
    end try;
    try
        candidates cat:= [P:P in Poles(EE!r)|Degree(P) eq 1];
    catch e
        dummy := 0;
    end try;
    if #candidates eq 0 then
        try candidates := Places(EE,1); catch e dummy := 0; end try;
    end if;
    print "PRODUCT_PARAM_START","family",family,"genus",Genus(EE),
          "degree1_candidates",#candidates;
    if #candidates eq 0 then
        print "PRODUCT_PARAM_NO_DEGREE1_PLACE",family;
        try
            degree2 := Places(EE,2);
            print "PRODUCT_PARAM_DEGREE2_PLACES",family,#degree2;
            if #degree2 gt 0 then
                D2 := Divisor(degree2[1]);
                V1,h1 := RiemannRochSpace(D2);
                V2,h2 := RiemannRochSpace(2*D2);
                assert Dimension(V1) eq 3 and Dimension(V2) eq 5;
                bas := [h1(b):b in Basis(V1)];
                pairs := [];
                for i in [1..3] do for j in [i..3] do
                    Append(~pairs,<i,j>);
                end for; end for;
                rows := [Eltseq((bas[z[1]]*bas[z[2]]) @@ h2):z in pairs];
                ker := Kernel(Matrix(Q,rows));
                assert Dimension(ker) eq 1;
                coeff := Eltseq(Basis(ker)[1]);
                PP<X0,X1,X2> := ProjectiveSpace(Q,2);
                xx := [X0,X1,X2];
                conicpoly := &+[coeff[k]*xx[pairs[k][1]]*xx[pairs[k][2]]
                                :k in [1..6]];
                con := Curve(PP,conicpoly);
                print "PRODUCT_PARAM_ANTICANONICAL_CONIC",family,conicpoly;
                ok,pt := HasRationalPoint(con);
                print "PRODUCT_PARAM_ANTICANONICAL_RATIONAL",family,ok,pt;
            end if;
        catch e0
            print "PRODUCT_PARAM_ANTICANONICAL_FAILED",family,e0`Object;
        end try;
        try
            CE := Curve(EE);
            cmap := Conic(CE);
            con := Codomain(cmap);
            print "PRODUCT_PARAM_CONIC",family,con;
            try
                ok,pt := HasRationalPoint(con);
                print "PRODUCT_PARAM_CONIC_RATIONAL",family,ok,pt;
            catch e2
                print "PRODUCT_PARAM_CONIC_RATIONAL_TEST_FAILED",family,e2`Object;
            end try;
        catch e
            print "PRODUCT_PARAM_CONIC_MODEL_FAILED",family,e`Object;
        end try;
        continue;
    end if;
    D := Divisor(candidates[1]);
    z,images := Parametrization(EE,D);
    print "PRODUCT_PARAM","family",family,"parameter",z,
          "images_count",#images,"images",images;
end for;

quit;
