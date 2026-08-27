//////////////////////////////////////////////////////////////////////
// Local Jacobian data for the genus-6 trigonal quotient of the
// contact-30 C3-root cover.
//
// This reconstructs the quotient directly in characteristic p and asks
// Magma for the L-polynomial of its normalized function field.  An
// irreducible degree-12 Frobenius polynomial at one good prime proves that
// the rational Jacobian has no nontrivial Q-isogeny factor.
//
//   magma -b code/contact30_c3root_genus6_local_jacobian.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(8*10^9);

if not assigned Primes then
    // p=13 is a fast good-reduction certificate.  Several larger primes
    // are much slower in Magma's general function-field L-polynomial code.
    Primes := [13];
elif Type(Primes) eq MonStgElt then
    Primes := [ StringToInteger(s) : s in Split(Primes, ",") ];
end if;
if not assigned ComputeGalois then
    ComputeGalois := false;
elif Type(ComputeGalois) eq MonStgElt then
    ComputeGalois := StringToLower(ComputeGalois) in
        { "true", "t", "yes", "y", "1" };
end if;
if not assigned ComputeClasses then
    ComputeClasses := false;
elif Type(ComputeClasses) eq MonStgElt then
    ComputeClasses := StringToLower(ComputeClasses) in
        { "true", "t", "yes", "y", "1" };
end if;

function QuotientPolynomial(k)
    K<z> := FunctionField(k);
    PT<T> := PolynomialRing(K);
    M<r> := ext<K | T^2-z*T+(5*z-7)/3>;

    RR := r;
    t := (5*RR^2-20*RR+19)/(RR^2-5);
    Y := -2*(5*RR^2-22*RR+25)/(RR^2-5);
    u := t^3;
    s := t^5+t^4+(M!5/2)*t^3+(M!1/2)*t
       - t*(t-M!1/2)*(t+1)*Y;
    C := (u^2+1)/(2*u);
    den := u^6+6*u^4*s-2*u^4+15*u^3*s-u*s^3+u^2;
    num := 15*u^5+90*u^4+20*u^3*s-6*u^2*s^2+231*u^3
       +2*u^2*s-15*u*s^2+90*u^2-20*u*s+15*u-2*s;
    q := num/den;
    A := (s+q)/2;
    B := (15-s*q)/2;
    assert Eltseq(A)[2] eq 0 and Eltseq(B)[2] eq 0
        and Eltseq(C)[2] eq 0;
    aa := K!Eltseq(A)[1];
    bb := K!Eltseq(B)[1];
    cc := K!Eltseq(C)[1];
    PK<x> := PolynomialRing(K);
    return 2*x^3+(aa-3)*x^2+(bb+3)*x+(cc-1);
end function;

Z<X> := PolynomialRing(Integers());

print "CONTACT30 GENUS6 LOCAL JACOBIAN";
for p in Primes do
    if p in [2,3,5] then
        print "prime",p,"skipped_denominator";
        continue;
    end if;
    k := GF(p);
    try
        f := QuotientPolynomial(k);
        if not IsIrreducible(f) then
            print "prime",p,"bad_or_reducible_defining_cubic";
            continue;
        end if;
        F<a> := FunctionField(f);
        g := Genus(F);
        if g ne 6 then
            print "prime",p,"bad_reduction_genus",g;
            continue;
        end if;
        time L := LPolynomial(F);
        LZ := Z!L;
        // L(T)=det(1-Frob*T).  Reverse it to the monic characteristic
        // polynomial det(X-Frob).
        chi := &+[ Coefficient(LZ,12-i)*X^i : i in [0..12] ];
        fac := Factorization(chi);
        print "prime",p,"genus",g;
        print "L",LZ;
        print "frobenius",chi;
        print "frobenius_factor_degrees",
            [ <Degree(fe[1]),fe[2]> : fe in fac ];
        print "jacobian_order",Evaluate(LZ,1);
        if ComputeClasses then
            Kz := BaseRing(Parent(f));
            zz := Kz.1;
            known_z := [ k!2,k!5,k!(14/3),k!(32/7) ];
            known_places := [];
            print "known_fiber_places_degree_ramification";
            for z0 in known_z do
                dz := Divisor(F!(zz-z0));
                pls := [ pl : pl in Support(dz) |
                    Valuation(dz,pl) gt 0 and Degree(pl) eq 1 ];
                print z0,[ <Degree(pl),Valuation(dz,pl)> : pl in pls ];
                known_places cat:= pls;
            end for;
            assert #known_places eq 6;
            time Cl,fromCl,toCl := ClassGroup(F);
            print "full_divisor_class_group",Cl;
            P0 := known_places[1];
            classes := [ toCl(Divisor(known_places[i])-Divisor(P0))
                         : i in [2..#known_places] ];
            print "known_difference_class_coordinates",
                [ Eltseq(c) : c in classes ];
            print "known_difference_orders",[ Order(c) : c in classes ];
            // Exact trigonal-fiber relations reduce the first five points
            // to B-E and D-E.  Together with the extra point F-E these are
            // the three visible rational divisor directions.
            visible := [ classes[1]-classes[4],
                         classes[3]-classes[4],
                         classes[5]-classes[4] ];
            viscoords := [ Eltseq(c)[1] : c in visible ];
            print "visible_three_coordinates_BminusE_DminusE_FminusE",
                viscoords;
            N := Evaluate(LZ,1);
            if IsPrime(N) and &and[ GCD(Integers()!v,N) eq 1
                                    : v in viscoords ] then
                u1,u2,u3 := Explode([ Integers()!v : v in viscoords ]);
                iu1 := InverseMod(u1,N);
                relation_basis := Matrix(Integers(),3,3,[
                    N,0,0,
                    (-iu1*u2) mod N,1,0,
                    (-iu1*u3) mod N,0,1 ]);
                relation_lll := LLL(relation_basis);
                print "modular_relation_lattice_LLL_rows",relation_lll;
                relation_lattice := Lattice(relation_basis);
                print "shortest_modular_relation_squared_norm",
                    Minimum(relation_lattice);
                print "shortest_modular_relations",
                    [ Eltseq(v) : v in ShortestVectors(relation_lattice) ];
            end if;
        end if;
        if #fac eq 1 and fac[1][2] eq 1 then
            print "Q_simple_specialization_certificate",true;
            print "extension_factor_degrees",
                [ <d,[ <Degree(h[1]),h[2]> : h in
                    Factorization(Z!WeilPolynomialOverFieldExtension(chi,d)) ]>
                  : d in [2..12] ];
            if ComputeGalois then
                time Gal,Galmap := GaloisGroup(chi);
                print "frobenius_galois_group_order",#Gal;
                print "frobenius_galois_group",Gal;
            end if;
        end if;
    catch e
        print "prime",p,"failed",e`Object;
    end try;
end for;

quit;
