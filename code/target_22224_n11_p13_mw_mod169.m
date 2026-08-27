//////////////////////////////////////////////////////////////////////
// N11 rank-two quotient: local MW image modulo 13^2 for the exact
// parameter t=N/(13D).  This refines the forced denominator-depth chart.
//////////////////////////////////////////////////////////////////////
SetColumns(0);SetSeed(8);Q:=Rationals();Z:=Integers();p:=13;modulus:=p^2;
R<T>:=PolynomialRing(Q);fixed:=[Q!-960,Q!1800,Q!2535];d0:=Q!-1352;t0:=Q!1/13;
Rx:=(fixed[2]+d0*T^2)/(fixed[2]+d0*t0^2);Ry:=(fixed[3]+d0*T^2)/(fixed[3]+d0*t0^2);
C:=HyperellipticCurve(Rx*Ry);P0:=C![t0,Q!1,Q!1];Eraw,phi:=EllipticCurve(C,P0);
E,minmap:=MinimalModel(Eraw);gens:=Generators(E);free:=[g:g in gens|Order(g) eq 0];
TG,tmap:=TorsionSubgroup(E);tors:=[tmap(g):g in TG];assert #free eq 2 and #tors eq 4;

function PrimitiveCoords(P)
  if P eq E!0 then return [Z!0,Z!1,Z!0];end if;
  q:=[Q!P[i]:i in [1..3]];den:=LCM([Denominator(x):x in q]);v:=[Z!(x*den):x in q];
  g:=GCD([Abs(x):x in v]);if g gt 1 then v:=[x div g:x in v];end if;return v;
end function;
function Key(P,m)
  v:=PrimitiveCoords(P);j:=0;for i in [1..3] do if GCD(v[i],p) eq 1 then j:=i;break;end if;end for;
  assert j ne 0;u:=InverseMod(v[j] mod m,m);return <(v[1]*u) mod m,(v[2]*u) mod m,(v[3]*u) mod m>;
end function;
function Val(n)
  if n eq 0 then return 99;end if;v:=0;n:=Abs(n);while n mod p eq 0 do v+:=1;n div:=p;end while;return v;
end function;
function KeyPeriod(H,m)
  A:=AssociativeArray();per:=0;P:=E!0;
  for j in [0..1000] do k:=Key(P,m);if IsDefined(A,k) then per:=j-A[k];break;end if;A[k]:=j;P+:=H;end for;
  return per,A;
end function;

candidates:=free cat [free[1]+free[2]] cat [free[i]+tors[j]:i in [1..2],j in [2..4]];
period:=0;G:=candidates[1];keys:=AssociativeArray();
for H in candidates do per,A:=KeyPeriod(H,modulus);print "BASE_TRY",H,per;
  if per gt period then period:=per;G:=H;keys:=A;end if;end for;
print "N11_P13_MW_MOD169_START","E",E,"gens",gens,"base",G,"period",period;
logs:=[];cyclic:=true;
for H in free cat tors do k:=Key(H,modulus);if not IsDefined(keys,k) then cyclic:=false;Append(~logs,-1);else Append(~logs,keys[k]);end if;end for;
print "LOCAL_LOGS","free",logs[1..2],"torsion",logs[3..6],"cyclic",cyclic;
if not cyclic then
  k2:=Key(2*free[1],modulus);print "TWO_G1_IN_BASE",IsDefined(keys,k2),IsDefined(keys,k2) select keys[k2] else -1;
  if IsDefined(keys,k2) and logs[2] ne -1 and &and[logs[i] ne -1:i in [3..6]] then
    log2g1:=keys[k2];out2:=Open("results/target_22224_n11_p13_mw_mod169.tsv","w");
    fprintf out2,"period\tlog_2g1\tlog_g2\tlog_t0\tlog_t1\tlog_t2\tlog_t3\tepsilon\tallowed\n";
    for eps in [0..1] do allowed2:=[];amb2:=[];P:=(eps eq 0) select E!0 else free[1];
      for j in [0..period-1] do
        v:=PrimitiveCoords(P);N:=-162345*v[1]+76*v[2]+938508450*v[3];D:=-617*v[1]+76*v[2]+3411910*v[3];
        vn:=Val(N);vd:=Val(D);if Maximum(0,1+vd-vn) ge 2 then Append(~allowed2,j);end if;
        if vn ge 2 and vd ge 2 then Append(~amb2,<j,vn,vd>);end if;P+:=G;
      end for;
      print "COSET_ALLOWED","eps",eps,allowed2,"count",#allowed2,"ambiguous",amb2;
      fprintf out2,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",period,log2g1,logs[2],logs[3],logs[4],logs[5],logs[6],eps,Join([Sprint(x):x in allowed2],",");
    end for;delete out2;
  end if;
end if;
if cyclic then
  allowed:=[];ambiguous:=[];
  P:=E!0;
  for j in [0..period-1] do
    v:=PrimitiveCoords(P);N:=-162345*v[1]+76*v[2]+938508450*v[3];
    D:=-617*v[1]+76*v[2]+3411910*v[3];vn:=Val(N);vd:=Val(D);denval:=Maximum(0,1+vd-vn);
    if denval ge 2 then Append(~allowed,j);end if;
    if vn ge 2 and vd ge 2 then Append(~ambiguous,<j,vn,vd>);end if;P+:=G;
  end for;
  print "ALLOWED_LOCAL_CLASSES",allowed,"count",#allowed;print "AMBIGUOUS",ambiguous;
  out:=Open("results/target_22224_n11_p13_mw_mod169.tsv","w");
  fprintf out,"period\tlog_g1\tlog_g2\tlog_t0\tlog_t1\tlog_t2\tlog_t3\tallowed\n";
  fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",period,logs[1],logs[2],logs[3],logs[4],logs[5],logs[6],Join([Sprint(x):x in allowed],",");delete out;
end if;
print "N11_P13_MW_MOD169_DONE";quit;
