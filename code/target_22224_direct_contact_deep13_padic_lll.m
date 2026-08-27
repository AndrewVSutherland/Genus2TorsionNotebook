//////////////////////////////////////////////////////////////////////
// Simultaneous/common-denominator reconstruction on the certified
// 11^5 x 13^5 direct-contact incidence bank.
//////////////////////////////////////////////////////////////////////
SetColumns(0);SetSeed(1);Z:=Integers();Q:=Rationals();
root:="results/";
if not assigned IncludedPrimes then IncludedPrimes:=[];
elif Type(IncludedPrimes) eq MonStgElt then IncludedPrimes:=[StringToInteger(s):s in Split(IncludedPrimes,",")];end if;
if not assigned output_file then output_file:=root cat "target_22224_direct_contact_deep13_padic_lll.tsv";end if;
if not assigned log_file then log_file:=root cat "target_22224_direct_contact_deep13_padic_lll.log";end if;
SetLogFile(log_file:Overwrite:=true);

function ReadTriples(path,offset)
 L:=Split(Read(path),"\n");S:={};
 for i in [2..#L] do if #L[i] eq 0 then continue;end if;z:=Split(L[i],"\t");
  if #z ge offset+3 then Include(~S,<StringToInteger(z[offset+1]),StringToInteger(z[offset+2]),StringToInteger(z[offset+3])>);end if;
 end for;return Setseq(S);
end function;
function CRT2(a,m,b,n)
 return (a+m*(((b-a)*InverseMod(m,n)) mod n)) mod (m*n);
end function;
function Allowed(p)
 B:=ReadTriples(Sprintf("%otarget_22224_direct_contact_deep13_boundary_p%o.tsv",root,p),1);
 T:=ReadTriples(Sprintf("%otarget_22224_direct_contact_deep13_p%o.tsv",root,p),1);
 return Seqset(B cat T);
end function;
function PassAt(q,p,A)
 if &or[Denominator(z) mod p eq 0:z in q] then return true,true;end if;
 r:=<Z!((Numerator(z)*InverseMod(Denominator(z),p)) mod p):z in q>;
 return r in A,false;
end function;

m11:=11^5;m13:=13^5;M:=m11*m13;
A11:=ReadTriples(root cat "target_22224_direct_contact_deep13_padic_tangent_p11.tsv",0);
A13:=ReadTriples(root cat "target_22224_direct_contact_deep13_padic_tangent_p13.tsv",0);
ps:=[17,19,23,29,31,37,41,43];AA:=AssociativeArray();for p in ps do AA[p]:=Allowed(p);end for;
good:=[<1,<0,0,0>>];
for p in IncludedPrimes do
 T:=ReadTriples(Sprintf("%otarget_22224_direct_contact_deep13_p%o.tsv",root,p),1);next:=[];
 for g in good do for t in T do
  Append(~next,<g[1]*p,<CRT2(g[2][j],g[1],t[j],p):j in [1..3]>>);
 end for;end for;good:=next;
end for;
Gmod:=good[1][1];TotalM:=M*Gmod;
coeffs:=[];
for c1,c2,c3,c4 in [-1..1] do c:=[c1,c2,c3,c4];
 if c eq [0,0,0,0] then continue;end if;
 j:=Min([i:i in [1..4]|c[i] ne 0]);if c[j] eq 1 then Append(~coeffs,c);end if;
end for;

out:=Open(output_file,"w");fprintf out,"u_num\tu_den\tt_num\tt_den\tw_num\tw_den\n";
pairs:=0;vectors:=0;nonzeroD:=0;unique:=0;survive:=0;seen:={};passes:=[0:i in ps];nonunits:=[0:i in ps];
print "PADIC_LLL_START","masks",#A11,#A13,"included",IncludedPrimes,"good_states",#good,
      "lattices",#A11*#A13*#good,"M",TotalM,"coeffs",#coeffs;
for a in A11 do for b in A13 do
 rdeep:=[CRT2(a[j],m11,b[j],m13):j in [1..3]];
 for gs in good do pairs+:=1;r:=[CRT2(rdeep[j],M,gs[2][j],Gmod):j in [1..3]];
 B:=Matrix(Z,4,4,[TotalM,0,0,0, 0,TotalM,0,0, 0,0,TotalM,0, r[1],r[2],r[3],1]);R:=LLL(B);
 for c in coeffs do vectors+:=1;v:=Vector(Z,c)*R;if v[4] eq 0 then continue;end if;nonzeroD+:=1;
  if v[4] lt 0 then v:=-v;end if;g:=GCD([Abs(v[i]):i in [1..4]]);v:=Vector(Z,[v[i] div g:i in [1..4]]);
  q:=[Q!v[i]/v[4]:i in [1..3]];if q[2] eq 0 then continue;end if;k:=<q[1],q[2],q[3]>;
  if k in seen then continue;end if;Include(~seen,k);unique+:=1;ok:=true;
  for i in [1..#ps] do yes,nu:=PassAt(q,ps[i],AA[ps[i]]);if not yes then ok:=false;break;end if;
   passes[i]+:=1;if nu then nonunits[i]+:=1;end if;
  end for;
  if not ok then continue;end if;survive+:=1;
  fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\n",Numerator(q[1]),Denominator(q[1]),Numerator(q[2]),Denominator(q[2]),Numerator(q[3]),Denominator(q[3]);
 end for;
 end for;
end for;end for;delete out;
print "PADIC_LLL_DONE","pairs",pairs,"vectors",vectors,"nonzeroD",nonzeroD,"unique",unique,
      "passes",[<ps[i],passes[i],nonunits[i]>:i in [1..#ps]],"survivors",survive,"output",output_file;
UnsetLogFile();quit;
