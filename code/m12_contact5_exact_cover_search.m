//////////////////////////////////////////////////////////////////////
//  Exact point-contact-5 cover of the full split M(12) surface.
//
//  Put
//      a = (1-z^2)/(4*(r+1)),  w = 2*(r+1)/(1-z),
//  in the M(12) model and move the rational Weierstrass point w to
//  infinity.  This gives an odd quintic f=f_{r,z} and an exact point
//  D12 of order 12.
//
//  A rational point (u,c) gives a contact-5 class [(u,c)-infinity]
//  exactly when
//      f(x) = h(x)^2 + lc(f)*(x-u)^5
//  for a quadratic h with h(u)=c.  Writing
//      Ai = f^(i)(u)/i!,  Q2 = 4*A0*A2-A1^2,
//  and eliminating the other two coefficients of h gives the exact
//  one-dimensional cover
//      E3 = 8*A0^2*A3-A1*Q2 = 0,
//      E4 = 64*A0^3*A4-Q2^2 = 0,
//      c^2 = A0 != 0.
//  Thus a rational cover point produces D5 of exact order 5, and
//  D12+D5 has exact order 60.
//
//  The finite masks below impose this exact point-contact cover, not
//  merely 5 | #J(F_p).  Bad/boundary parameter residues are passed.
//
//  Typical runs (from torsion_jac):
//      magma -b mode:="finite" prime_bound:=19 \
//          code/m12_contact5_exact_cover_search.m
//      magma -b mode:="search" height:=50 prime_bound:=19 \
//          code/m12_contact5_exact_cover_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then
    mode := "search";
end if;
if not assigned height then
    height := 50;
elif Type(height) eq MonStgElt then
    height := StringToInteger(height);
end if;
if not assigned prime_bound then
    prime_bound := 19;
elif Type(prime_bound) eq MonStgElt then
    prime_bound := StringToInteger(prime_bound);
end if;
if not assigned max_hits then
    max_hits := 20;
elif Type(max_hits) eq MonStgElt then
    max_hits := StringToInteger(max_hits);
end if;
if not assigned max_print then
    max_print := 30;
elif Type(max_print) eq MonStgElt then
    max_print := StringToInteger(max_print);
end if;
if not assigned progress_interval then
    progress_interval := 1000000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;
if not assigned simplicity_bound then
    simplicity_bound := 97;
elif Type(simplicity_bound) eq MonStgElt then
    simplicity_bound := StringToInteger(simplicity_bound);
end if;
if not assigned geometry_height then
    geometry_height := 50;
elif Type(geometry_height) eq MonStgElt then
    geometry_height := StringToInteger(geometry_height);
end if;
if not assigned max_genus_degree then
    max_genus_degree := 12;
elif Type(max_genus_degree) eq MonStgElt then
    max_genus_degree := StringToInteger(max_genus_degree);
end if;
if not assigned MemGB then
    MemGB := 8;
elif Type(MemGB) eq MonStgElt then
    MemGB := StringToInteger(MemGB);
end if;
SetMemoryLimit(MemGB*10^9);

Q := Rationals();
Z := Integers();
P<x> := PolynomialRing(Q);
PT<T> := PolynomialRing(Q);

function RationalParametersOfHeight(B)
    vals := [];
    seen := {};
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num,den) ne 1 then
                continue;
            end if;
            q := Q!num/Q!den;
            key := Sprint(q);
            if key notin seen then
                Include(~seen,key);
                Append(~vals,q);
            end if;
        end for;
    end for;
    return vals;
end function;

function OddM12Quintic(r,z)
    if r eq -1 or z in {Q!-1,Q!1} then
        return false,P!0,Q!0,Q!0,Q!0;
    end if;
    a := (1-z^2)/(4*(r+1));
    if a eq 0 then
        return false,P!0,a,Q!0,Q!0;
    end if;
    S := a*x^2-x+r;
    h12 := (x-r)*(S+1);
    W := h12^2+4*a*x^2*S*(S+1);
    w := 2*(r+1)/(1-z);
    if w eq 0 or Evaluate(S+1,w) ne 0 then
        return false,P!0,a,Q!0,Q!0;
    end if;
    f := P!0;
    for i in [0..Degree(W)] do
        for j in [0..i] do
            f +:= Coefficient(W,i)*Binomial(i,j)*w^(i-j)*x^(6-j);
        end for;
    end for;
    X12 := -1/w;
    Y12 := Evaluate(h12,0)*X12^3;
    return true,f,a,X12,Y12;
end function;

function OddM12QuinticFinite(rr,zz,F)
    PF<X> := PolynomialRing(F);
    if rr eq -1 or zz^2 eq 1 then
        return false,PF!0;
    end if;
    a := (1-zz^2)/(F!4*(rr+1));
    if a eq 0 then
        return false,PF!0;
    end if;
    S := a*X^2-X+rr;
    h12 := (X-rr)*(S+1);
    W := h12^2+4*a*X^2*S*(S+1);
    w := 2*(rr+1)/(1-zz);
    if w eq 0 or Evaluate(S+1,w) ne 0 then
        return false,PF!0;
    end if;
    f := PF!0;
    for i in [0..Degree(W)] do
        for j in [0..i] do
            f +:= Coefficient(W,i)*Binomial(i,j)*w^(i-j)*X^(6-j);
        end for;
    end for;
    return true,f;
end function;

function IsGoodQuintic(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function ContactCovariants(f)
    K := BaseRing(Parent(f));
    A0 := f;
    A1 := Derivative(f);
    A2 := Derivative(f,2)/(K!2);
    A3 := Derivative(f,3)/(K!6);
    A4 := Derivative(f,4)/(K!24);
    Q2 := 4*A0*A2-A1^2;
    E3 := 8*A0^2*A3-A1*Q2;
    E4 := 64*A0^3*A4-Q2^2;
    return E3,E4;
end function;

function ContactAbscissas(f)
    E3,E4 := ContactCovariants(f);
    g := GCD(E3,E4);
    us := [];
    if IsZero(g) then
        return us,g;
    end if;
    for rt in Roots(g) do
        // A repeated rational root is still a genuine cover point; its
        // multiplicity only says that the projection is ramified there.
        u := rt[1];
        A0 := Evaluate(f,u);
        if A0 ne 0 and IsSquare(A0) then
            Append(~us,u);
        end if;
    end for;
    return us,g;
end function;

function ContactQuadratic(f,u,c)
    K := BaseRing(Parent(f));
    A1 := Evaluate(Derivative(f),u);
    A2 := Evaluate(Derivative(f,2),u)/(K!2);
    d := A1/(2*c);
    e := (A2-d^2)/(2*c);
    h := c+d*(Parent(f).1-u)+e*(Parent(f).1-u)^2;
    return h,d,e;
end function;

procedure ContactSelfTest()
    // Planted positive control from the standard contact-5/order-20
    // family.  At t=2, u=0 and c=1 are a known smooth contact.
    t := Q!2;
    h := 1+t*x+((t^2-1)/2)*x^2;
    f := h^2-((t+1)^4/4)*x^5;
    assert IsGoodQuintic(f);
    us,g := ContactAbscissas(f);
    assert Q!0 in us;
    hrec,d,e := ContactQuadratic(f,Q!0,Q!1);
    assert hrec eq h;
    assert f-hrec^2 eq LeadingCoefficient(f)*x^5;
    print "CONTACT_SELF_TEST_PASS","t",t,"u",0,
          "gcd_degree",Degree(g),"root_multiplicity",
          [rt[2] : rt in Roots(g) | rt[1] eq 0][1];
end procedure;

function ResidueKey(rr,zz,p)
    return rr*p+zz;
end function;

function FiniteContactMask(p)
    F := GF(p);
    allowed := {};
    bad := {};
    good := 0;
    contacts := 0;
    samples := [];
    for ri in [0..p-1] do
        for zi in [0..p-1] do
            key := ResidueKey(ri,zi,p);
            ok,f := OddM12QuinticFinite(F!ri,F!zi,F);
            if not ok or not IsGoodQuintic(f) then
                Include(~bad,key);
                continue;
            end if;
            // Homogenize the contact coordinate u.  The generic degrees
            // are 12 and 16.  A degree drop is a point on the projective
            // cover boundary (not a valid finite contact), so pass it
            // conservatively instead of using it to exclude a lift.
            E3,E4 := ContactCovariants(f);
            if Degree(E3) lt 12 or Degree(E4) lt 16 then
                Include(~bad,key);
                continue;
            end if;
            good +:= 1;
            us,g := ContactAbscissas(f);
            if #us gt 0 then
                Include(~allowed,key);
                contacts +:= #us;
                if #samples lt 8 then
                    Append(~samples,<ri,zi,Z!us[1]>);
                end if;
            end if;
        end for;
    end for;
    return allowed,bad,good,contacts,samples;
end function;

function ParameterResidues(params,plist)
    table := [];
    for q in params do
        row := [];
        for p in plist do
            if Denominator(q) mod p eq 0 then
                Append(~row,-1);
            else
                F := GF(p);
                Append(~row,Z!(F!Numerator(q)/F!Denominator(q)));
            end if;
        end for;
        Append(~table,row);
    end for;
    return table;
end function;

function PassesContactMasks(rrs,zzs,plist,allowed,badres)
    good_checked := 0;
    boundary_primes := [];
    for i in [1..#plist] do
        p := plist[i];
        rr := rrs[i];
        zz := zzs[i];
        if rr eq -1 or zz eq -1 then
            Append(~boundary_primes,p);
            continue;
        end if;
        key := ResidueKey(rr,zz,p);
        if key in badres[p] then
            Append(~boundary_primes,p);
            continue;
        end if;
        good_checked +:= 1;
        if key notin allowed[p] then
            return false,p,good_checked,boundary_primes;
        end if;
    end for;
    return true,0,good_checked,boundary_primes;
end function;

function IntegralModelPolynomial(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L,Denominator(Coefficient(f,i)));
    end for;
    return P!(L^2*f),L;
end function;

function FrobeniusPolynomial(C,p)
    ef := EulerFactor(C,p);
    d := Degree(ef);
    return &+[Q!Coefficient(ef,i)*T^(d-i) : i in [0..d]];
end function;

function FullSimplicityCertificate(C,plist)
    for p in plist do
        try
            Phi := FrobeniusPolynomial(C,p);
            fac := Factorization(Phi);
            if #fac ne 1 or fac[1][2] ne 1 or Degree(fac[1][1]) ne 4 then
                continue;
            end if;
            Gal := GaloisGroup(Phi);
            desc := "";
            try
                desc := TransitiveGroupDescription(Gal);
            catch e2
                desc := "unknown";
            end try;
            if Order(Gal) eq 8 and desc eq "D(4)" then
                return true,"D4",p,Phi;
            end if;
        catch e
            continue;
        end try;
    end for;
    for p in plist do
        try
            Phi := FrobeniusPolynomial(C,p);
            fac := Factorization(Phi);
            if Degree(Phi) ne 4 or #fac ne 1 or fac[1][2] ne 1 or
               Degree(fac[1][1]) ne 4 then
                continue;
            end if;
            K<pi> := NumberField(Phi);
            power_ok := true;
            for n in [2..12] do
                if Degree(MinimalPolynomial(pi^n)) lt 4 then
                    power_ok := false;
                    break;
                end if;
            end for;
            if power_ok then
                return true,"root_power",p,Phi;
            end if;
        catch e
            continue;
        end try;
    end for;
    return false,"none",0,PT!0;
end function;

//////////////////////////////////////////////////////////////////////
//  Geometry in the normalized rational-root chart.
//
//  Put b=a*w and xi=w*X.  Then z=2*b-1, r=w*(1-b)-1,
//  D12 is based at xi=-1, and (after a square scaling of y)
//
//    F = L*(L*A^2+4*b*(1+xi)^2*(w*L-xi^2)),
//    L = b+(2*b-1)*xi,  A=xi+w*(1+b*xi).
//
//  The open smooth degree-five chart is
//      b*w*(b-1)*(2*b-1) != 0.
//////////////////////////////////////////////////////////////////////

function NormalizedRootChartQuintic(b,w)
    L := b+(2*b-1)*x;
    A := x+w*(1+b*x);
    return L*(L*A^2+4*b*(1+x)^2*(w*L-x^2));
end function;

function StripLinearBoundaries(f,qs)
    counts := [];
    if IsZero(f) then
        return f,[0 : q in qs];
    end if;
    for q in qs do
        count := 0;
        while TotalDegree(GCD(f,q)) gt 0 do
            f := ExactQuotient(f,q);
            count +:= 1;
        end while;
        Append(~counts,count);
    end for;
    return f,counts;
end function;

procedure RunNormalizedGeometry()
    Rbw<b0,w0> := PolynomialRing(Q,2,"grevlex");
    RX<xi> := PolynomialRing(Rbw);
    L := b0+(2*b0-1)*xi;
    A := xi+w0*(1+b0*xi);
    F := L*(L*A^2+4*b0*(1+xi)^2*(w0*L-xi^2));
    E3x,E4x := ContactCovariants(F);

    R3<B,W,U> := PolynomialRing(Q,3,"grevlex");
    to3 := hom<Rbw -> R3 | B,W>;
    E3 := &+[to3(Coefficient(E3x,i))*U^i : i in [0..Degree(E3x)]];
    E4 := &+[to3(Coefficient(E4x,i))*U^i : i in [0..Degree(E4x)]];

    print "NORMALIZED_CONTACT_GEOMETRY";
    print "F degree_xi",Degree(F),"terms",#Terms(F);
    print "E3 degree_u",Degree(E3,U),"total_degree",TotalDegree(E3),
          "terms",#Terms(E3);
    print "E4 degree_u",Degree(E4,U),"total_degree",TotalDegree(E4),
          "terms",#Terms(E4);

    E3open,c3 := StripLinearBoundaries(E3,[B,W,B-1,2*B-1]);
    E4open,c4 := StripLinearBoundaries(E4,[B,W,B-1,2*B-1]);
    print "OPEN_BOUNDARY_POWERS E3",c3,"E4",c4;
    common := GCD(E3open,E4open);
    print "COMMON_FACTOR total_degree",TotalDegree(common),
          "terms",#Terms(common);
    for fe in Factorization(common) do
        print "COMMON_FACTOR_PART","multiplicity",fe[2],
              "degree",TotalDegree(fe[1]),"terms",#Terms(fe[1]);
        if #Terms(fe[1]) le 80 then print fe[1]; end if;
    end for;
    if TotalDegree(common) gt 0 then
        E3open := ExactQuotient(E3open,common);
        E4open := ExactQuotient(E4open,common);
    end if;

    print "OPEN E3 factors";
    for fe in Factorization(E3open) do
        print "E3_FACTOR","multiplicity",fe[2],
              "degree",TotalDegree(fe[1]),"degree_u",Degree(fe[1],U),
              "terms",#Terms(fe[1]);
        if #Terms(fe[1]) le 80 then print fe[1]; end if;
    end for;
    print "OPEN E4 factors";
    for fe in Factorization(E4open) do
        print "E4_FACTOR","multiplicity",fe[2],
              "degree",TotalDegree(fe[1]),"degree_u",Degree(fe[1],U),
              "terms",#Terms(fe[1]);
        if #Terms(fe[1]) le 80 then print fe[1]; end if;
    end for;

    R2<BB,WW> := PolynomialRing(Q,2,"grevlex");
    to2 := hom<R3 -> R2 | BB,WW,R2!0>;
    bwto2 := hom<Rbw -> R2 | BB,WW>;
    Disc := R2!bwto2(Discriminant(F));
    DiscOpen,disc_counts :=
        StripLinearBoundaries(Disc,[BB,WW,BB-1,2*BB-1]);
    print "DISCRIMINANT_BOUNDARY_POWERS",disc_counts;
    print "DISCRIMINANT_OPEN_FACTORS";
    discfac := Factorization(DiscOpen);
    for fe in discfac do
        print "DISC_FACTOR","multiplicity",fe[2],
              "degree",TotalDegree(fe[1]),"terms",#Terms(fe[1]);
        if #Terms(fe[1]) le 80 then print fe[1]; end if;
    end for;

    print "SPECIAL_CONTACT_SLICES";
    for uv in [Q!0,Q!-1] do
        e3s := R2!to2(Evaluate(E3open,[B,W,R3!uv]));
        e4s := R2!to2(Evaluate(E4open,[B,W,R3!uv]));
        gs := GCD(e3s,e4s);
        print "SLICE","u",uv,"E3_degree",TotalDegree(e3s),
              "E4_degree",TotalDegree(e4s),
              "gcd_degree",TotalDegree(gs),"gcd_terms",#Terms(gs);
        for fe in Factorization(gs) do
            print "SLICE_FACTOR","u",uv,"multiplicity",fe[2],
                  "degree",TotalDegree(fe[1]),"terms",#Terms(fe[1]),
                  fe[1];
        end for;
    end for;
    print "SLICE_U0_SQUARE_COVER v^2 = w*(w+4*b)";

    // The equations have much smaller degree in w than in the contact
    // coordinate u.  Eliminate w and keep the useful (b,u) projection.
    Rbu<Bp,Up> := PolynomialRing(Q,2,"grevlex");
    toBU := hom<R3 -> Rbu | Bp,Rbu!0,Up>;
    print "RESULTANT_ELIMINATE_W_BEGIN",
          "degrees_w",<Degree(E3open,W),Degree(E4open,W)>;
    time Res3 := Resultant(E3open,E4open,W);
    if IsZero(Res3) then
        print "RESULTANT_ZERO_AFTER_COMMON_FACTOR_REMOVAL";
        return;
    end if;
    assert Degree(Res3,W) eq 0;
    Res := Rbu!toBU(Res3);
    ResOpen,res_counts := StripLinearBoundaries(Res,[Bp,Bp-1,2*Bp-1]);
    print "RESULTANT_BOUNDARY_POWERS",res_counts;
    print "RESULTANT_OPEN degree",TotalDegree(ResOpen),
          "terms",#Terms(ResOpen);
    resfac := Factorization(ResOpen);
    print "RESULTANT_FACTORS",#resfac;

    lc3 := Rbu!toBU(LeadingCoefficient(E3open,{2}));
    lc4 := Rbu!toBU(LeadingCoefficient(E4open,{2}));
    infinity_boundary := lc3*lc4;

    genuine := [];
    for fe in resfac do
        g := fe[1];
        divides_infinity := TotalDegree(GCD(g,infinity_boundary)) eq TotalDegree(g);
        print "RES_FACTOR","multiplicity",fe[2],
              "degree",TotalDegree(g),"degree_b",Degree(g,Bp),
              "degree_u",Degree(g,Up),"terms",#Terms(g),
              "w_infinity_boundary",divides_infinity;
        if #Terms(g) le 100 then print g; end if;
        if not divides_infinity and TotalDegree(g) gt 0 then
            Append(~genuine,g);
        end if;
    end for;

    print "PROJECTED_COMPONENT_GEOMETRY";
    for g in genuine do
        print "GENUINE_COMPONENT","degree",TotalDegree(g),
              "degree_b",Degree(g,Bp),"degree_u",Degree(g,Up),
              "terms",#Terms(g);
        if TotalDegree(g) gt max_genus_degree then
            print " geometry_skipped_degree_limit",max_genus_degree;
            continue;
        end if;
        try
            Caff := Curve(AffineSpace(Rbu),g);
            Cp := ProjectiveClosure(Caff);
            print " projective_degree",Degree(Cp),
                  "normalized_genus",Genus(Cp);
        catch eg
            print " geometry_failed",eg`Object;
        end try;
    end for;

    // Search the plane projection by specializing b.  First recover rational
    // u, then recover every common rational w-root of E3 and E4.  No height
    // bound is imposed on u or w; the in-box counter records H(b),H(w)<=B.
    bvals := RationalParametersOfHeight(geometry_height);
    QU<Y> := PolynomialRing(Q);
    QW<Zz> := PolynomialRing(Q);
    seen_bu := {};
    seen_bw := {};
    projected_bu := 0;
    recovered_bw := 0;
    recovered_bw_in_box := 0;
    exact_cover_points := 0;
    exact60 := 0;
    for bv in bvals do
        if bv in {Q!0,Q!1,Q!1/2} then continue; end if;
        for g in genuine do
            ev := hom<Rbu -> QU | bv,Y>;
            gy := ev(g);
            if IsZero(gy) then
                continue;
            end if;
            for rt in Roots(gy) do
                uv := Q!rt[1];
                keybu := Sprint(<bv,uv>);
                if keybu in seen_bu then continue; end if;
                Include(~seen_bu,keybu);
                projected_bu +:= 1;
                recover := hom<R3 -> QW | bv,Zz,uv>;
                e3w := recover(E3open);
                e4w := recover(E4open);
                gw := GCD(e3w,e4w);
                if IsZero(gw) then continue; end if;
                for wrt in Roots(gw) do
                    wv := Q!wrt[1];
                    if wv eq 0 then continue; end if;
                    keybw := Sprint(<bv,wv,uv>);
                    if keybw in seen_bw then continue; end if;
                    Include(~seen_bw,keybw);
                    recovered_bw +:= 1;
                    if Max(Abs(Numerator(wv)),Denominator(wv)) le
                       geometry_height then
                        recovered_bw_in_box +:= 1;
                    end if;
                    fn := NormalizedRootChartQuintic(bv,wv);
                    if not IsGoodQuintic(fn) then continue; end if;
                    us,gg := ContactAbscissas(fn);
                    if uv notin us then continue; end if;
                    exact_cover_points +:= 1;
                    rv := wv*(1-bv)-1;
                    zv := 2*bv-1;
                    print "NORMALIZED_EXACT_COVER_POINT","b",bv,
                          "w",wv,"u",uv,"r",rv,"z",zv;
                    ok,forig,av,X12,Y12 := OddM12Quintic(rv,zv);
                    if not ok or not IsGoodQuintic(forig) then continue; end if;
                    xcontact := uv/wv;
                    uorig,gorig := ContactAbscissas(forig);
                    if xcontact notin uorig then
                        print " NORMALIZATION_CONTACT_MISMATCH",xcontact,uorig;
                        continue;
                    end if;
                    corig := SquareRoot(Evaluate(forig,xcontact));
                    fI,scale := IntegralModelPolynomial(forig);
                    try
                        C := HyperellipticCurve(fI);
                        J := Jacobian(C);
                        D12 := J![x-X12,scale*Y12];
                        D5 := J![x-xcontact,scale*corig];
                        P60 := D12+D5;
                        orders := <Order(D12),Order(D5),Order(P60)>;
                        if orders eq <12,5,60> then
                            exact60 +:= 1;
                            G,phi := TorsionSubgroup(J);
                            invs := Invariants(G);
                            print "NORMALIZED_ORDER60_HIT","b",bv,"w",wv,
                                  "u",uv,"r",rv,"z",zv,
                                  "torsion",invs,"orders",orders;
                        else
                            print " NORMALIZED_ORDER_MISMATCH",orders;
                        end if;
                    catch ee
                        print " NORMALIZED_EXACT_ERROR",ee`Object;
                    end try;
                end for;
            end for;
        end for;
    end for;
    print "NORMALIZED_PROJECTION_SCAN_SUMMARY","height_b",geometry_height,
          "b_values",#bvals,"projected_bu",projected_bu,
          "recovered_bwu",recovered_bw,
          "recovered_bwu_in_Hbw_box",recovered_bw_in_box,
          "exact_cover_points",exact_cover_points,"exact60",exact60;
end procedure;

//////////////////////////////////////////////////////////////////////
//  The low-degree projected component P8 and its missing square cover.
//////////////////////////////////////////////////////////////////////

procedure RunP8SquareGeometry()
    Rbu<B,U> := PolynomialRing(Q,2,"grevlex");
    P8 :=
        B^4*U^4 - 4*B^4*U^3 - B^3*U^4 - 14*B^4*U^2
        + 12*B^3*U^3 + 1/4*B^2*U^4 - 12*B^4*U
        + 28*B^3*U^2 - 17/2*B^2*U^3 - 3*B^4
        + 18*B^3*U - 49/4*B^2*U^2 + 7/4*B*U^3
        + 3*B^3 - 3*B^2*U + 1/2*B*U^2 + 3/4*B^2
        - 3/4*B*U + 1/4*U^2;

    P2<Bh,Uh,Zh> := ProjectiveSpace(Q,2);
    P8h := Homogenization(Evaluate(P8,[Bh,Uh]),Zh);
    C8 := Curve(P2,P8h);
    print "P8_SQUARE_GEOMETRY";
    print "P8 degree",Degree(C8),"irreducible",IsIrreducible(C8),
          "normalized_genus",Genus(C8);
    sing := SingularPoints(C8);
    print "P8 rational_singular_points",#sing;
    for pp in sing do
        print " P8_SINGULAR",pp,"multiplicity",Multiplicity(C8,pp);
    end for;

    Ccon,mc := Conic(C8);
    print "P8_NORMALIZATION_CONIC",Ccon;
    locally := IsLocallySolvable(Ccon);
    haspt,pc := HasRationalPoint(Ccon);
    print "P8_CONIC locally_solvable",locally,
          "has_rational_point",haspt;
    if not haspt then
        print "P8_NO_Q_NORMALIZATION_POINT";
        print "The rational plane points at infinity do not lift to a degree-1 place.";
        print "P8_DIAGONAL_CONIC 13*X^2 + 169*Y^2 = 432*Z^2";
        print "P8_Q3_OBSTRUCTION left valuation even, right valuation odd";
        return;
    end if;
    print "P8_CONIC_POINT",pc;

    // Obtain a P1 parametrization of the original singular plane model.
    phi_found := false;
    phi := IdentityMap(C8);
    try
        phi := ParametrizeOrdinaryCurve(C8);
        phi_found := true;
    catch ep0
        pts8 := Points(C8 : Bound := 1000);
        print "P8_PLANE_POINTS_BOUND1000",#pts8,pts8;
        for pp in pts8 do
            try
                phi := Parametrization(C8,pp);
                phi_found := true;
                break;
            catch ep1
                continue;
            end try;
        end for;
    end try;
    if not phi_found then
        print "P8_PARAMETRIZATION_FAILED_DESPITE_CONIC_POINT";
        return;
    end if;
    eqs := DefiningEquations(phi);
    print "P8_PARAMETRIZATION coordinate_degrees",
          [Degree(e) : e in eqs];

    Kt<t> := FunctionField(Q);
    bt := Kt!Evaluate(eqs[1],[t,Kt!1])/
          (Kt!Evaluate(eqs[3],[t,Kt!1]));
    ut := Kt!Evaluate(eqs[2],[t,Kt!1])/
          (Kt!Evaluate(eqs[3],[t,Kt!1]));
    print "P8_b(t)",bt;
    print "P8_u(t)",ut;

    // Rebuild E3,E4 in (b,w,u) and recover the rational w branches over Q(t).
    R3<b,w,u> := PolynomialRing(Q,3,"grevlex");
    RX<xi> := PolynomialRing(R3);
    Ln := b+(2*b-1)*xi;
    An := xi+w*(1+b*xi);
    Fn := Ln*(Ln*An^2+4*b*(1+xi)^2*(w*Ln-xi^2));
    A0 := R3!Evaluate(Fn,u);
    A1 := R3!Evaluate(Derivative(Fn),u);
    A2 := R3!Evaluate(Derivative(Fn,2),u)/2;
    A3 := R3!Evaluate(Derivative(Fn,3),u)/6;
    A4 := R3!Evaluate(Derivative(Fn,4),u)/24;
    QQ2 := 4*A0*A2-A1^2;
    E3 := ExactQuotient(8*A0^2*A3-A1*QQ2,b);
    E4 := ExactQuotient(64*A0^3*A4-QQ2^2,b);

    KW<Wvar> := PolynomialRing(Kt);
    recover := hom<R3 -> KW | bt,Wvar,ut>;
    gw := GCD(recover(E3),recover(E4));
    print "P8_RECOVERY_W degree",Degree(gw),
          "factor_degrees",[Degree(fe[1]) : fe in Factorization(gw)];

    branches := [];
    for fe in Factorization(gw) do
        if Degree(fe[1]) eq 1 then
            Append(~branches,-Coefficient(fe[1],0)/Coefficient(fe[1],1));
        end if;
    end for;
    print "P8_RATIONAL_W_BRANCHES",#branches;
    for wv in branches do
        print "P8_w(t)",wv;
        Lv := bt+(2*bt-1)*ut;
        Av := ut+wv*(1+bt*ut);
        fv := Lv*(Lv*Av^2+4*bt*(1+ut)^2*(wv*Lv-ut^2));
        num := Numerator(fv);
        den := Denominator(fv);
        sqraw := num*den;
        fac := Factorization(sqraw);
        H := Parent(sqraw)!LeadingCoefficient(sqraw);
        for ff in fac do
            if IsOdd(ff[2]) then
                H *:= ff[1]/LeadingCoefficient(ff[1]);
            end if;
        end for;
        print "P8_SQUARE_RAW degrees",<Degree(num),Degree(den)>,
              "raw_factor_degrees",[<Degree(ff[1]),ff[2]> : ff in fac];
        print "P8_SQUARE_CLASS degree",Degree(H),"factorization",
              Factorization(H);
        if Degree(H) ge 3 then
            CH := HyperellipticCurve(H);
            print "P8_SQUARE_COVER genus",Genus(CH),
                  "points_bound_1000",Points(CH : Bound := 1000);
            for pp in [3,5,7,11,13,17,19,23,29,31] do
                try
                    CHp := ChangeRing(CH,GF(pp));
                    print " P8_SQUARE_LOCAL",pp,"points",#CHp;
                catch eloc
                    print " P8_SQUARE_LOCAL",pp,"bad";
                end try;
            end for;
        else
            print "P8_SQUARE_COVER_LOW_DEGREE",H;
        end if;
    end for;
end procedure;

ContactSelfTest();
if mode eq "geometry" then
    RunNormalizedGeometry();
    quit;
end if;
if mode eq "p8" then
    RunP8SquareGeometry();
    quit;
end if;

// Larger primes have much thinner contact masks, so apply them first.
// At p=19 the degree-16 covariant drops degree (95=0), hence that
// projective chart is boundary and contributes no rigorous filtering.
plist := Reverse([p : p in PrimesUpTo(prime_bound) | p gt 5 and p ne 19]);
allowed := AssociativeArray();
badres := AssociativeArray();

print "M(12) exact point-contact-5 cover search";
print "mode",mode,"height",height,"prime_bound",prime_bound,
      "primes",plist;
print "cover equations E3=E4=0, c^2=f(u) != 0";
for p in plist do
    good_allowed,bad,good,contacts,samples := FiniteContactMask(p);
    allowed[p] := good_allowed;
    badres[p] := bad;
    print "CONTACT_MASK",p,"good_pairs",good,
          "allowed_pairs",#good_allowed,"contact_points",contacts,
          "bad_boundary",#bad,"samples",samples;
end for;

if mode eq "finite" then
    print "DONE finite exact-contact masks";
    quit;
end if;

params := RationalParametersOfHeight(height);
residues := ParameterResidues(params,plist);
checked := 0;
mask_survivors := 0;
smooth_survivors := 0;
exact_cover_points := 0;
order60_hits := 0;
cyclic60_exact_hits := 0;
order60_containing_hits := 0;
simple_hits := 0;
printed := 0;
kill_counts := AssociativeArray();

for ir in [1..#params] do
    r := params[ir];
    for iz in [1..#params] do
        z := params[iz];
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "PROGRESS",checked,"mask_survivors",mask_survivors,
                  "smooth",smooth_survivors,
                  "cover_points",exact_cover_points,
                  "order60",order60_hits;
        end if;
        if r eq -1 or z in {Q!-1,Q!1} then
            continue;
        end if;
        pass,killp,good_checked,boundary_primes :=
            PassesContactMasks(residues[ir],residues[iz],plist,
                               allowed,badres);
        if not pass then
            if IsDefined(kill_counts,killp) then
                kill_counts[killp] +:= 1;
            else
                kill_counts[killp] := 1;
            end if;
            continue;
        end if;
        mask_survivors +:= 1;
        ok,f,a,X12,Y12 := OddM12Quintic(r,z);
        if not ok or not IsGoodQuintic(f) then
            continue;
        end if;
        smooth_survivors +:= 1;
        us,g := ContactAbscissas(f);
        if #us eq 0 then
            continue;
        end if;

        for u in us do
            c := SquareRoot(Evaluate(f,u));
            h5,d,e := ContactQuadratic(f,u,c);
            if f-h5^2 ne LeadingCoefficient(f)*(x-u)^5 then
                print "INTERNAL_CONTACT_MISMATCH",r,z,u;
                continue;
            end if;
            exact_cover_points +:= 1;
            if printed lt max_print then
                print "EXACT_COVER_POINT","r",r,"z",z,"a",a,
                      "u",u,"c",c,"good_checked",good_checked,
                      "boundary_primes",boundary_primes;
                printed +:= 1;
            end if;

            fI,L := IntegralModelPolynomial(f);
            try
                C := HyperellipticCurve(fI);
                J := Jacobian(C);
                D12 := J![x-X12,L*Y12];
                D5 := J![x-u,L*c];
                o12 := Order(D12);
                o5 := Order(D5);
                P60 := D12+D5;
                o60 := Order(P60);
                if o12 ne 12 or o5 ne 5 or o60 ne 60 then
                    print "ORDER_MISMATCH","r",r,"z",z,"u",u,
                          "orders",<o12,o5,o60>;
                    continue;
                end if;
                G,phi := TorsionSubgroup(J);
                invs := Invariants(G);
                order60_hits +:= 1;
                if invs eq [60] then
                    cyclic60_exact_hits +:= 1;
                    print "CYCLIC_Z60_HIT","r",r,"z",z,"u",u;
                else
                    order60_containing_hits +:= 1;
                    print "ORDER60_CONTAINING_HIT","r",r,"z",z,
                          "u",u,"torsion",invs;
                end if;
                splist := [p : p in PrimesUpTo(simplicity_bound) | p ge 3];
                simple,stype,sp,Phi := FullSimplicityCertificate(C,splist);
                if simple then
                    simple_hits +:= 1;
                end if;
                print "ORDER60_HIT","r",r,"z",z,"a",a,"u",u,
                      "torsion",invs,"orders",<o12,o5,o60>,
                      "simple",simple,"simple_type",stype,
                      "simple_prime",sp;
                print "  f_integral =",fI;
                print "  D12 =",D12;
                print "  D5 =",D5;
                print "  P60 =",P60;
                if simple then
                    print "  Frobenius =",Phi;
                end if;
                if order60_hits ge max_hits then
                    break ir;
                end if;
            catch err
                print "EXACT_ERROR","r",r,"z",z,"u",u,
                      "message",err`Object;
            end try;
        end for;
    end for;
end for;

print "DONE M12 exact point-contact-5 cover";
print "parameters",#params,"checked_pairs",checked,
      "mask_survivors",mask_survivors,
      "smooth_survivors",smooth_survivors,
      "exact_cover_points",exact_cover_points,
      "order60_hits",order60_hits,
      "cyclic60_exact_hits",cyclic60_exact_hits,
      "order60_containing_hits",order60_containing_hits,
      "simple_hits",simple_hits;
print "kill_counts";
for p in Sort([k : k in Keys(kill_counts)]) do
    print p,kill_counts[p];
end for;
quit;
