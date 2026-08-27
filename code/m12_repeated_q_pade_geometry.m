//////////////////////////////////////////////////////////////////////
// Low-memory Pade equations for a point of order dividing 10 on M(12).
//
// At a non-Weierstrass point P=(u,y0), write t=x-u and
//
//   y/y0 = c0+c1*t+... , c0=1.
//
// A function A(x)+B(x)y, deg(A)<=5 and deg(B)<=2, can vanish to order
// 10 at P iff a quadratic B(t)=1+r*t+s*t^2 makes the coefficients in
// degrees 6,7,8,9 vanish.  The first two equations recover r,s on the
// chart Delta=c5^2-c4*c6 != 0; the last two give just two equations in
// (b,w,u).  No square root occurs because y/y0=sqrt(F(u+t)/F(u)).
// Rational points P are the signed double cover y0^2=F(u).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned mode then mode := "summary"; end if;
require mode in {"summary","modular","search"}:
        "mode must be summary, modular, or search";
if not assigned p then p := 7;
elif Type(p) eq MonStgElt then p := StringToInteger(p); end if;
if not assigned MemGB then MemGB := 3;
elif Type(MemGB) eq MonStgElt then MemGB := StringToInteger(MemGB); end if;
if not assigned do_primary then do_primary := true;
elif Type(do_primary) eq MonStgElt then
    do_primary := do_primary in {"true","True","1"};
end if;
SetMemoryLimit(MemGB*10^9);

function PrimitivePolynomial(f)
    R := Parent(f);
    if IsFinite(BaseRing(R)) or f eq 0 then return f; end if;
    cs := Coefficients(f);
    den := LCM([Denominator(a):a in cs]);
    nums := [Integers()!(den*a):a in cs];
    con := GCD(nums); if con eq 0 then con:=1; end if;
    return R!((Rationals()!den/Rationals()!Abs(con))*f);
end function;

function PadeData(k)
    R<b,w,u> := PolynomialRing(k,3,"grevlex");
    K := FieldOfFractions(R);
    P<t> := PolynomialRing(K);
    x := u+t;
    L := b+(2*b-1)*x;
    H := x+w*(1+b*x);
    Fshift := P!(L*(L*H^2+4*b*(1+x)^2*(w*L-x^2)));
    f := [K!Coefficient(Fshift,i):i in [0..5]];
    F0 := f[1];
    c := [K!1];
    for n in [1..9] do
        fn := n le 5 select f[n+1]/F0 else K!0;
        conv := n eq 1 select K!0 else &+[c[i+1]*c[n-i+1]:i in [1..n-1]];
        Append(~c,(fn-conv)/2);
    end for;
    // n=6,7: c_n+r*c_(n-1)+s*c_(n-2)=0.
    Delta := c[6]^2-c[5]*c[7];       // c5^2-c4*c6
    rr := (-c[7]*c[6]+c[5]*c[8])/Delta;
    ss := (-c[6]*c[8]+c[7]^2)/Delta;
    assert c[7]+rr*c[6]+ss*c[5] eq 0;
    assert c[8]+rr*c[7]+ss*c[6] eq 0;
    E8 := c[9]+rr*c[8]+ss*c[7];
    E9 := c[10]+rr*c[9]+ss*c[8];
    N8 := PrimitivePolynomial(R!Numerator(E8));
    N9 := PrimitivePolynomial(R!Numerator(E9));
    ND := PrimitivePolynomial(R!Numerator(Delta));
    Fu := PrimitivePolynomial(R!Numerator(F0));
    return R,Fu,ND,N8,N9,c,rr,ss;
end function;

function IsUnitIdeal(I)
    return &or[g ne 0 and TotalDegree(g) eq 0:g in Basis(I)];
end function;

procedure Report(label,I)
    B := Basis(I);
    print "IDEAL",label,"unit",IsUnitIdeal(I),"basis_len",#B,
          "degrees",[TotalDegree(f):f in B];
    if IsUnitIdeal(I) then return; end if;
    try
        d,ds := Dimension(I);
        print "IDEAL",label,"dimension",d,"component_degrees",ds;
    catch err print "IDEAL",label,"dimension_failed",err`Object; end try;
    try print "IDEAL",label,"degree",Degree(I);
    catch err print "IDEAL",label,"degree_failed",err`Object; end try;
end procedure;

k := mode eq "summary" select Rationals() else GF(p);
if mode ne "summary" then
    require IsPrime(p) and p notin {2,3,5}: "use a prime away from 2,3,5";
end if;
R,Fu,ND,N8,N9,c,rr,ss := PadeData(k);
print "M12_REPEATED_Q_PADE_FORMULAS_PASS";
print "SHAPES","Fu",<TotalDegree(Fu),#Terms(Fu)>,
      "Delta",<TotalDegree(ND),#Terms(ND)>,
      "N8",<TotalDegree(N8),#Terms(N8)>,
      "N9",<TotalDegree(N9),#Terms(N9)>;

if mode eq "summary" then
    print "FACTORIZATION_N8",
          [<TotalDegree(a[1]),#Terms(a[1]),a[2]>:a in Factorization(N8)];
    print "FACTORIZATION_N9",
          [<TotalDegree(a[1]),#Terms(a[1]),a[2]>:a in Factorization(N9)];
    quit;
end if;

b:=R.1; w:=R.2; u:=R.3;
I := ideal<R|N8,N9>;
Report("raw",I);
factors := [b,w,b-1,2*b-1,Fu,ND];
names := ["b","w","b-1","2b-1","Fu","Delta"];
J := I;
for i in [1..#factors] do
    nf := NormalForm(factors[i],Basis(J));
    print "SAT_BEGIN",names[i],"degree",TotalDegree(nf),"terms",#Terms(nf);
    time J := Saturation(J,ideal<R|nf>);
    Report(Sprintf("sat_%o",names[i]),J);
    if IsUnitIdeal(J) then break; end if;
end for;
if do_primary and not IsUnitIdeal(J) then
    print "PRIMARY_BEGIN";
    try
        time comps := PrimaryDecomposition(J);
        print "PRIMARY_COUNT",#comps;
        for i in [1..#comps] do
            Report(Sprintf("component_%o",i),comps[i]);
            Bc := Basis(comps[i]);
            if #Bc le 12 and &and[#Terms(g) le 100:g in Bc] then
                print "COMPONENT_BASIS",i,Bc;
            end if;
        end for;
    catch err print "PRIMARY_FAILED",err`Object; end try;
end if;
print "M12_REPEATED_Q_PADE_MODULAR_DONE";
quit;
