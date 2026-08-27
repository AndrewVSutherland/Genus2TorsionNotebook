//////////////////////////////////////////////////////////////////////
//  Exact point-contact-5 cover on the Q4-root rational-Weierstrass
//  surface inside M(12).
//
//  In the completed-square M(12) model,
//
//      W=(T+1)*Q4,
//      Q4=(x-r)^2*(T+1)+4*a*x^2*T,
//      T=a*x^2-x+r.
//
//  Force a rational root v of Q4.  Put d=v-r, alpha=a*v^2, and
//  k=alpha/d.  Then
//
//      Q4(v)=4*alpha^2+(d^2-4*d)*alpha+d^2*(1-d),
//
//  and, after alpha=k*d,
//
//      Q4(v)/d^2 = (2*k-1)^2-(1-k)*d.
//
//  Write D=1-k and q=2*k-1.  Away from the degenerate boundary,
//
//      d=q^2/D,  a=k*q^2/(D*v^2),  r=v-q^2/D.
//
//  Move v to infinity with X=v/(x-v).  Thus old x=v*(X+1)/X.
//  Put
//
//      P=4*k*D^2*X^2+(2*k*q^2-v*D)*X+k*q^2,
//      R=-q^2*D*X^2+(2*k*q^2-v*D)*X+k*q^2,
//      M=q^2*X+v*D,
//      Q=M^2*P+4*k*q^2*D*(X+1)^2*R.
//
//  Directly, P=D*X^2*(T+1), R=D*X^2*T, M=D*X*(old x-r),
//  and the square-equivalent odd model is
//
//      F=P*Q=D^4*X^6*W(v+v/X).
//
//  The X^4 coefficient of Q cancels identically, so Q is cubic.
//  Thus F has generic factor type [2,3], only one rational 2-direction,
//  and exact cyclic Z/60 remains possible.  The marked order-12 point is
//  X=-1 with y=-D^2*r*(r+1), up to sign.
//
//  Modes (from torsion_jac):
//
//      magma -b mode:="finite" prime_bound:=19 \
//          code/m12_q4root_contact5_exact_search.m
//
//      magma -b mode:="search" height:=30 prime_bound:=19 max_hits:=20 \
//          code/m12_q4root_contact5_exact_search.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned mode then
    mode := "search";
end if;
if mode notin {"finite","search"} then
    error "mode must be finite or search";
end if;
if not assigned height then
    height := 30;
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
    progress_interval := 500000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;
if not assigned simplicity_bound then
    simplicity_bound := 97;
elif Type(simplicity_bound) eq MonStgElt then
    simplicity_bound := StringToInteger(simplicity_bound);
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
PT<Tfrob> := PolynomialRing(Q);

function HeightRationals(B)
    vals := [];
    for den in [1..B] do
        for num in [-B..B] do
            if GCD(num,den) eq 1 then
                Append(~vals,Q!num/Q!den);
            end if;
        end for;
    end for;
    return vals;
end function;

function Q4RootQuintic(k,v)
    D := 1-k;
    q := 2*k-1;
    // k=0 or q=0 forces a=0; D=0 cannot solve the Q4-root
    // equation on this affine chart.  v=0 is a marked discriminant
    // boundary and cannot be used for the Mobius normalization.
    if k eq 0 or D eq 0 or q eq 0 or v eq 0 then
        return false,P!0,Q!0,Q!0,Q!0,Q!0,Q!0;
    end if;
    d := q^2/D;
    alpha := k*d;
    a := k*q^2/(D*v^2);
    r := v-d;

    P2 := 4*k*D^2*x^2+(2*k*q^2-v*D)*x+k*q^2;
    R2 := -q^2*D*x^2+(2*k*q^2-v*D)*x+k*q^2;
    M := q^2*x+v*D;
    G3 := M^2*P2+4*k*q^2*D*(1+x)^2*R2;
    f := P2*G3;
    y12 := -D^2*r*(r+1);
    return true,f,a,r,d,alpha,y12;
end function;

function Q4RootQuinticFinite(kk,vv,F)
    PF<X> := PolynomialRing(F);
    D := 1-kk;
    q := 2*kk-1;
    if kk eq 0 or D eq 0 or q eq 0 or vv eq 0 then
        return false,PF!0;
    end if;
    P2 := 4*kk*D^2*X^2+(2*kk*q^2-vv*D)*X+kk*q^2;
    R2 := -q^2*D*X^2+(2*kk*q^2-vv*D)*X+kk*q^2;
    M := q^2*X+vv*D;
    G3 := M^2*P2+4*kk*q^2*D*(1+X)^2*R2;
    return true,P2*G3;
end function;

function GoodQuintic(f)
    return Degree(f) eq 5 and Discriminant(f) ne 0;
end function;

function FactorDegrees(f)
    return Sort([Degree(fe[1]) : fe in Factorization(f)]);
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
    // Standard planted contact-5/order-20 positive control.
    t0 := Q!2;
    h0 := 1+t0*x+((t0^2-1)/2)*x^2;
    f0 := h0^2-((t0+1)^4/4)*x^5;
    assert GoodQuintic(f0);
    us,g := ContactAbscissas(f0);
    assert Q!0 in us;
    hrec,drec,erec := ContactQuadratic(f0,Q!0,Q!1);
    assert hrec eq h0;
    assert f0-hrec^2 eq LeadingCoefficient(f0)*x^5;
    print "CONTACT_SELF_TEST_PASS","u",0,"gcd_degree",Degree(g);
end procedure;

procedure ChartSelfTest()
    // Algebraic control for the Q4-root parametrization at (k,v)=(2,1).
    k0 := Q!2;
    v0 := Q!1;
    ok,f,a,r,d,alpha,y12 := Q4RootQuintic(k0,v0);
    assert ok and GoodQuintic(f);
    assert FactorDegrees(f) eq [2,3];
    assert Evaluate(f,-1) eq y12^2;
    assert <a,r,d,alpha> eq <Q!-18,Q!10,Q!-9,Q!-18>;

    S := a*x^2-x+r;
    Q4 := (x-r)^2*(S+1)+4*a*x^2*S;
    assert Evaluate(Q4,v0) eq 0;
    assert alpha eq a*v0^2 and d eq v0-r;
    print "Q4ROOT_SELF_TEST_PASS","k",k0,"v",v0,
          "a",a,"r",r,"factor_degrees",FactorDegrees(f),
          "D12_square",Evaluate(f,-1);
end procedure;

function ResidueKey(kk,vv,p)
    return kk*p+vv;
end function;

function FiniteContactMask(p)
    F := GF(p);
    allowed := {};
    bad := {};
    good := 0;
    contacts := 0;
    samples := [];
    for ki in [0..p-1] do
        for vi in [0..p-1] do
            key := ResidueKey(ki,vi,p);
            ok,f := Q4RootQuinticFinite(F!ki,F!vi,F);
            if not ok or not GoodQuintic(f) then
                Include(~bad,key);
                continue;
            end if;
            E3,E4 := ContactCovariants(f);
            // Homogenize E3 and E4 to their separate generic degrees 12
            // and 16.  A common point at u=infinity occurs only when both
            // leading coefficients vanish.  A degree drop in just one
            // equation is not a projective contact point.
            if Degree(E3) lt 12 and Degree(E4) lt 16 then
                Include(~bad,key);
                continue;
            end if;
            good +:= 1;
            us,g := ContactAbscissas(f);
            if #us gt 0 then
                Include(~allowed,key);
                contacts +:= #us;
                if #samples lt 8 then
                    Append(~samples,<ki,vi,Z!us[1]>);
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

function PassesMasks(kres,vres,plist,allowed,badres)
    good_checked := 0;
    boundary_primes := [];
    for i in [1..#plist] do
        p := plist[i];
        kv := kres[i];
        vv := vres[i];
        if kv eq -1 or vv eq -1 then
            Append(~boundary_primes,p);
            continue;
        end if;
        key := ResidueKey(kv,vv,p);
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

function IntegralModel(f)
    L := 1;
    for i in [0..Degree(f)] do
        L := LCM(L,Denominator(Coefficient(f,i)));
    end for;
    return P!(L^2*f),L;
end function;

function FrobeniusPolynomial(C,p)
    ef := EulerFactor(C,p);
    dd := Degree(ef);
    return &+[Q!Coefficient(ef,i)*Tfrob^(dd-i) : i in [0..dd]];
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

ContactSelfTest();
ChartSelfTest();

// The generic leading coefficients of (E3,E4) for a quintic are
// (5*c5^3,95*c5^4).  At p=19 only E4 drops degree; E3 stays nonzero at
// u=infinity, so p=19 still gives a valid and useful affine mask.
plist := Reverse([p : p in PrimesUpTo(prime_bound) | p gt 5]);
allowed := AssociativeArray();
badres := AssociativeArray();

print "M12 Q4-root exact point-contact5 search";
print "mode",mode,"height",height,"prime_bound",prime_bound,
      "primes",plist;
print "generic odd factor type [2,3]; exact cyclic [60] is permitted";
for p in plist do
    good_allowed,bad,good,contacts,samples := FiniteContactMask(p);
    allowed[p] := good_allowed;
    badres[p] := bad;
    print "CONTACT_MASK",p,"good_pairs",good,
          "allowed_pairs",#good_allowed,"contact_points",contacts,
          "bad_boundary",#bad,"samples",samples;
end for;

if mode eq "finite" then
    print "DONE finite Q4-root contact masks";
    quit;
end if;

params := HeightRationals(height);
residues := ParameterResidues(params,plist);
checked := 0;
mask_survivors := 0;
smooth_survivors := 0;
exact_cover_points := 0;
order60_hits := 0;
cyclic60_hits := 0;
containing_hits := 0;
simple_hits := 0;
printed := 0;
kill_counts := AssociativeArray();

for ik in [1..#params] do
    k := params[ik];
    for iv in [1..#params] do
        v := params[iv];
        checked +:= 1;
        if progress_interval gt 0 and checked mod progress_interval eq 0 then
            print "PROGRESS",checked,"mask",mask_survivors,
                  "smooth",smooth_survivors,"cover",exact_cover_points,
                  "order60",order60_hits;
        end if;
        if k eq 0 or k eq 1 or 2*k-1 eq 0 or v eq 0 then
            continue;
        end if;
        pass,killp,good_checked,boundary_primes :=
            PassesMasks(residues[ik],residues[iv],plist,allowed,badres);
        if not pass then
            if IsDefined(kill_counts,killp) then
                kill_counts[killp] +:= 1;
            else
                kill_counts[killp] := 1;
            end if;
            continue;
        end if;
        mask_survivors +:= 1;

        ok,f,a,r,d,alpha,y12 := Q4RootQuintic(k,v);
        if not ok or not GoodQuintic(f) then
            continue;
        end if;
        smooth_survivors +:= 1;
        us,g := ContactAbscissas(f);
        for u in us do
            c := SquareRoot(Evaluate(f,u));
            h5,dh,eh := ContactQuadratic(f,u,c);
            if f-h5^2 ne LeadingCoefficient(f)*(x-u)^5 then
                print "INTERNAL_CONTACT_MISMATCH",k,v,u;
                continue;
            end if;
            exact_cover_points +:= 1;
            if printed lt max_print then
                print "EXACT_COVER_POINT","k",k,"v",v,
                      "a",a,"r",r,"u",u,"c",c,
                      "factor_degrees",FactorDegrees(f),
                      "good_checked",good_checked,
                      "boundary_primes",boundary_primes;
                printed +:= 1;
            end if;

            fI,Lscale := IntegralModel(f);
            try
                C := HyperellipticCurve(fI);
                J := Jacobian(C);
                D12 := J![x+1,Lscale*y12];
                D5 := J![x-u,Lscale*c];
                o12 := Order(D12);
                o5 := Order(D5);
                P60 := D12+D5;
                o60 := Order(P60);
                if o12 ne 12 or o5 ne 5 or o60 ne 60 then
                    print "ORDER_MISMATCH","k",k,"v",v,"u",u,
                          "orders",<o12,o5,o60>;
                    continue;
                end if;
                Gtors,phi := TorsionSubgroup(J);
                invs := Invariants(Gtors);
                order60_hits +:= 1;
                if invs eq [60] then
                    cyclic60_hits +:= 1;
                    print "CYCLIC_Z60_HIT","k",k,"v",v,"u",u;
                else
                    containing_hits +:= 1;
                    print "ORDER60_CONTAINING_HIT","k",k,"v",v,
                          "u",u,"torsion",invs;
                end if;
                splist := [p : p in PrimesUpTo(simplicity_bound) | p ge 3];
                simple,stype,sp,Phi := FullSimplicityCertificate(C,splist);
                if simple then
                    simple_hits +:= 1;
                end if;
                print "ORDER60_HIT","k",k,"v",v,"a",a,"r",r,
                      "u",u,"torsion",invs,
                      "orders",<o12,o5,o60>,
                      "simple",simple,"simple_type",stype,
                      "simple_prime",sp;
                print "  f_integral =",fI;
                print "  D12 =",D12;
                print "  D5 =",D5;
                print "  P60 =",P60;
                if simple then print "  Frobenius =",Phi; end if;
                if order60_hits ge max_hits then
                    break ik;
                end if;
            catch err
                print "EXACT_ERROR","k",k,"v",v,"u",u,
                      "message",err`Object;
            end try;
        end for;
    end for;
end for;

print "DONE M12 Q4-root exact point-contact5 search";
print "parameters",#params,"checked_pairs",checked,
      "mask_survivors",mask_survivors,
      "smooth_survivors",smooth_survivors,
      "exact_cover_points",exact_cover_points,
      "order60_hits",order60_hits,
      "cyclic60_hits",cyclic60_hits,
      "order60_containing_hits",containing_hits,
      "simple_hits",simple_hits;
print "kill_counts";
for p in Sort([k : k in Keys(kill_counts)]) do
    print p,kill_counts[p];
end for;

quit;
