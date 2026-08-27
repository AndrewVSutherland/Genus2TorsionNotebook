//////////////////////////////////////////////////////////////////////
// Fast rank-2 lattice walk on the P8 elliptic quotient.  Uses Magma's
// explicit Generators(E), avoiding the expensive generic MW-group call.
//////////////////////////////////////////////////////////////////////
SetColumns(0);
if not assigned N then N:=30; elif Type(N) eq MonStgElt then N:=StringToInteger(N); end if;
if not assigned log_file then log_file:=Sprintf("results/target_22224_p8_fast_lattice_N%o.log",N); end if;
if not assigned output_file then output_file:=Sprintf("results/target_22224_p8_fast_lattice_N%o.tsv",N); end if;
SetLogFile(log_file : Overwrite:=true);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);
fixed:=[Q!528,-726,-891]; d0:=Q!50;

function Primitive(vals)
    L:=LCM([Denominator(z):z in vals]); v:=[Z!(L*z):z in vals];
    g:=GCD([Abs(z):z in v]); if g gt 1 then v:=[z div g:z in v]; end if;
    for z in v do if z ne 0 then if z lt 0 then v:=[-w:w in v]; end if; break; end if; end for;
    return v;
end function;

function ThreeBound(v)
    for p in [17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127] do
        k:=GF(p); P<X>:=PolynomialRing(k);
        fp:=X*&*[X+(k!(z mod p))^2:z in v];
        if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
        nj:=Z!Evaluate(LPolynomial(HyperellipticCurve(fp)),1);
        if nj mod 3 ne 0 then return false,p,nj; end if;
    end for;
    return true,0,0;
end function;

Ra:=(fixed[1]+d0*T^2)/(fixed[1]+d0);
Rb:=(fixed[2]+d0*T^2)/(fixed[2]+d0);
C:=HyperellipticCurve(Ra*Rb); P0:=C![Q!1,Q!1,Q!1];
E,phi:=EllipticCurve(C,P0); Einv:=Inverse(phi);
Emin,minmap:=MinimalModel(E);
mingens:=Generators(Emin); allgens:=[Inverse(minmap)(g):g in mingens];
free:=[g:g in allgens|Order(g) eq 0];
TG,tmap:=TorsionSubgroup(Emin); tors:=[Inverse(minmap)(tmap(g)):g in TG];
assert #free eq 2;
print "P8_FAST_START","N",N,"E",E,"generators",allgens;
mult1:=[i*free[1]:i in [-N..N]]; mult2:=[j*free[2]:j in [-N..N]];
seen:={}; full:=0; survivors:=0; mapped:=0;
out:=Open(output_file,"w");
fprintf out,"i\tj\ttorsion\tt_num\tt_den\ta\tb\tc\td\tthree_survivor\tkill_prime\tkill_order\n";
for ii in [1..#mult1] do for jj in [1..#mult2] do for kk in [1..#tors] do
    ep:=mult1[ii]+mult2[jj]+tors[kk];
    try cp:=Einv(ep); catch err continue; end try;
    if cp[3] eq 0 then continue; end if;
    tt:=Q!(cp[1]/cp[3]); mapped+:=1; dd:=d0*tt^2;
    good:=true;
    for z in fixed do ok,s:=IsSquare((z+dd)/(z+d0)); if not ok then good:=false; break; end if; end for;
    if not good then continue; end if;
    v:=Primitive(fixed cat [dd]); key:=<v[1],v[2],v[3],v[4]>;
    if key in seen then continue; end if; Include(~seen,key); full+:=1;
    survives,p,nj:=ThreeBound(v); if survives then survivors+:=1; end if;
    fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
        ii-N-1,jj-N-1,kk,Numerator(tt),Denominator(tt),v[1],v[2],v[3],v[4],
        survives select 1 else 0,p,nj;
    print "P8_FULL","ij",ii-N-1,jj-N-1,"tor",kk,"t",tt,"tuple",v,
          "three",survives,"kill",p,nj;
    if survives then
        PP<X>:=PolynomialRing(Q); f:=X*&*[X+(Q!z)^2:z in v];
        inv:=Invariants(TorsionSubgroup(Jacobian(HyperellipticCurve(f))));
        print "P8_EXACT_SURVIVOR","torsion",inv,"f",f;
    end if;
end for; end for; end for;
delete out;
print "P8_FAST_DONE","mapped",mapped,"full",full,"three_survivors",survivors;
UnsetLogFile(); quit;
