//////////////////////////////////////////////////////////////////////
// Compute the image of F3's rank-three MW generators in E(Q_13)
// modulo 13^2, and impose v_13(den(t)) >= 2 on the exact rational
// function
//
//   t = N/(13 D),
//   N=-2023 X+130 Y+1077103 Z,
//   D=  429 X+ 10 Y+  22399 Z.
//////////////////////////////////////////////////////////////////////
SetColumns(0);SetSeed(3); Q:=Rationals();Z:=Integers();p:=13;modulus:=p^2;
R<T>:=PolynomialRing(Q);fixed:=[Q!-612,Q!34,Q!289];d0:=Q!-338;
Rx:=(fixed[1]+d0*T^2)/(fixed[1]+d0);Ry:=(fixed[2]+d0*T^2)/(fixed[2]+d0);
C:=HyperellipticCurve(Rx*Ry);P0:=C![Q!1,Q!1,Q!1];Eraw,phi:=EllipticCurve(C,P0);
E,minmap:=MinimalModel(Eraw);gens:=Generators(E);free:=[g:g in gens|Order(g) eq 0];
TG,tmap:=TorsionSubgroup(E);tors:=[tmap(g):g in TG];assert #free eq 3 and #tors eq 2;

function PrimitiveCoords(P)
  if P eq E!0 then return [Z!0,Z!1,Z!0];end if;
  q:=[Q!P[i]:i in [1..3]];den:=LCM([Denominator(x):x in q]);
  v:=[Z!(x*den):x in q];g:=GCD([Abs(x):x in v]);if g gt 1 then v:=[x div g:x in v];end if;
  return v;
end function;
function Key(P,m)
  v:=PrimitiveCoords(P);j:=0;
  for i in [1..3] do if GCD(v[i],p) eq 1 then j:=i;break;end if;end for;
  assert j ne 0;u:=InverseMod(v[j] mod m,m);
  return <(v[1]*u) mod m,(v[2]*u) mod m,(v[3]*u) mod m>;
end function;
function Val(n)
  if n eq 0 then return 99;end if;v:=0;n:=Abs(n);while n mod p eq 0 do v+:=1;n div:=p;end while;return v;
end function;

function KeyPeriod(H,m)
  A:=AssociativeArray();per:=0;
  for j in [0..1000] do k:=Key(j*H,m);if IsDefined(A,k) then per:=j-A[k];break;end if;A[k]:=j;end for;
  return per,A;
end function;
candidates:=free cat [free[i]+free[j]:i,j in [1..3]|i lt j] cat [free[i]+tors[2]:i in [1..3]];
period:=0;G:=candidates[1];keys:=AssociativeArray();
for H in candidates do per,A:=KeyPeriod(H,modulus);if per gt period then period:=per;G:=H;keys:=A;end if;end for;
print "F3_P13_MW_MOD169_START","E",E,"gens",gens,"base",G,"period",period;
assert period gt 0;
logs:=[];
for H in free cat tors do
  k:=Key(H,modulus);assert IsDefined(keys,k);Append(~logs,keys[k]);
end for;
print "LOCAL_LOGS","free",logs[1..3],"torsion_cosets",logs[4..5];
allowed:=[];ambiguous:=[];
for j in [0..period-1] do
  P:=j*G;v:=PrimitiveCoords(P);N:=-2023*v[1]+130*v[2]+1077103*v[3];
  D:=429*v[1]+10*v[2]+22399*v[3];vn:=Val(N);vd:=Val(D);
  denval:=Maximum(0,1+vd-vn);
  if denval ge 2 then Append(~allowed,j);end if;
  if vn ge 2 and vd ge 2 then Append(~ambiguous,<j,vn,vd>);end if;
end for;
print "ALLOWED_LOCAL_CLASSES",allowed,"count",#allowed;
print "AMBIGUOUS_HIGH_V",ambiguous;
// Resolve the sole 13^2-indeterminate class recursively in the formal
// direction.  At depth e, points in one scalar class differ by the local
// period 16*13^(e-1); if min(v(N),v(D))<e, both valuations are stable on
// the whole class.  Otherwise split it into its thirteen children.
maxdepth:=8; frontier:=[<x[1],2,period>:x in ambiguous]; resolved_allowed:=[];
for depth in [3..maxdepth] do
  nextfrontier:=[]; hist:=AssociativeArray();
  for c in frontier do
    j0:=c[1];oldperiod:=c[3];
    for q in [0..p-1] do
      j:=j0+oldperiod*q;P:=j*G;v:=PrimitiveCoords(P);
      N:=-2023*v[1]+130*v[2]+1077103*v[3];D:=429*v[1]+10*v[2]+22399*v[3];
      vn:=Val(N);vd:=Val(D);tag:="stable_fail";
      if Minimum(vn,vd) ge depth then
        tag:="ambiguous";Append(~nextfrontier,<j,depth,oldperiod*p>);
      elif Maximum(0,1+vd-vn) ge 2 then
        tag:="allowed";Append(~resolved_allowed,<j,oldperiod*p,depth,vn,vd>);
      end if;
      if IsDefined(hist,tag) then hist[tag]+:=1;else hist[tag]:=1;end if;
      print "FORMAL_CHILD","depth",depth,"j",j,"vn",vn,"vd",vd,"tag",tag;
    end for;
  end for;
  print "FORMAL_DEPTH_DONE",depth,"hist",hist,"next",#nextfrontier;
  frontier:=nextfrontier;if #frontier eq 0 then break;end if;
end for;
print "FORMAL_RESOLVED_ALLOWED",resolved_allowed;
print "FORMAL_UNRESOLVED",frontier;
if not assigned output_file then output_file:="results/target_22224_f3_p13_mw_mod169.tsv";end if;
out:=Open(output_file,"w");fprintf out,"period\tlog_g1\tlog_g2\tlog_g3\tlog_t0\tlog_t1\tallowed\n";
fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\t%o\n",period,logs[1],logs[2],logs[3],logs[4],logs[5],Join([Sprint(x):x in allowed],",");delete out;
print "F3_P13_MW_MOD169_DONE",output_file;quit;
