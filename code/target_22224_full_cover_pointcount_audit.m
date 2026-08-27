// Exhaustive, contact-free check of the good-reduction local obstruction.
// It enumerates the full A(2,2,2,8) signed cover and tests #J(F_p) directly.
SetColumns(0);
if not assigned PrimeList then PrimeList:=[11,13,17]; end if;

function CoverRadicands(v)
    a,b,c,d:=Explode(v);
    return [a*b*c*d,
            a*(a+b)*(a+c)*(a+d),
            b*(b+a)*(b+c)*(b+d),
            c*(c+a)*(c+b)*(c+d)];
end function;

function IsProjectiveRepresentative(v)
    for z in v do
        if z ne 0 then return z eq 1; end if;
    end for;
    return false;
end function;

for p in PrimeList do
    F:=GF(p); P<x>:=PolynomialRing(F);
    proj:=0; cover:=0; smooth:=0; div3:=0;
    keys:={}; keys3:={}; hist:=AssociativeArray(); samples:=[];
    for a,b,c,d in F do
        v:=[a,b,c,d];
        if not IsProjectiveRepresentative(v) then continue; end if;
        proj+:=1;
        rad:=CoverRadicands(v);
        if not &and[z ne 0 and IsSquare(z):z in rad] then continue; end if;
        cover+:=1;
        sq:=[z^2:z in v];
        if 0 in sq or #Set(sq) ne 4 then continue; end if;
        smooth+:=1;
        e1:=&+sq;
        e2:=&+[sq[i]*sq[j]:i,j in [1..4]|i lt j];
        e3:=&+[sq[i]*sq[j]*sq[k]:i,j,k in [1..4]|i lt j and j lt k];
        e4:=&*sq;
        key:=<Integers()!e1,Integers()!e2,Integers()!e3,Integers()!e4>;
        Include(~keys,key);
        f:=x*&*[x+z:z in sq];
        n:=Integers()!Evaluate(LPolynomial(HyperellipticCurve(f)),1);
        r:=n mod 3;
        if not IsDefined(hist,r) then hist[r]:=0; end if;
        hist[r]+:=1;
        if r eq 0 then
            div3+:=1; Include(~keys3,key);
            if #samples lt 12 then Append(~samples,<v,n,key>); end if;
        end if;
    end for;
    print "FULL_COVER_POINTCOUNT_AUDIT","p",p,"projective",proj,
          "cover",cover,"smooth",smooth,"curvekeys",#keys,
          "presentations_3div",div3,"curvekeys_3div",#keys3,
          "mod3_hist",hist,"samples",samples;
end for;
