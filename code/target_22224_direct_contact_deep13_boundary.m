//////////////////////////////////////////////////////////////////////
// Fresh p=11,13 boundary incidence for the direct contact quotient.
// No legacy Aaux equation is used.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(1);
if not assigned PrimeList then PrimeList:=[11,13];
elif Type(PrimeList) eq MonStgElt then
  PrimeList:=[StringToInteger(s):s in Split(PrimeList,",")];
end if;
if not assigned output_prefix then
  output_prefix:="results/target_22224_direct_contact_deep13_boundary";
end if;
if not assigned log_file then log_file:=output_prefix cat ".log"; end if;
SetLogFile(log_file:Overwrite:=true); Z:=Integers();

function DirectData(K,u,t,w)
  P<X>:=PolynomialRing(K); q:=X^2+u*X+t^2;
  h:=X^3+(1+3*u)/2*X^2+w*X+t^3;
  F:=ExactQuotient(h^2-q^3,X); return P,F,q,h;
end function;
function CoverRadicands(a,b,c,d)
  return [a*b*c*d,a*(a+b)*(a+c)*(a+d),
          b*(b+a)*(b+c)*(b+d),c*(c+a)*(c+b)*(c+d)];
end function;
function RootMagnitudesWithMultiplicity(F)
  fac:=Factorization(F); vals:=[];
  if not &and[Degree(z[1]) eq 1:z in fac] then return false,[]; end if;
  for z in fac do
    rr:=-Coefficient(z[1],0)/Coefficient(z[1],1);
    if not IsSquare(-rr) then return false,[]; end if;
    mag:=SquareRoot(-rr);
    for j in [1..z[2]] do Append(~vals,mag); end for;
  end for;
  return #vals eq 4,vals;
end function;
function SignChoices(mags)
  ans:=[[]];
  for z in mags do
    choices:=z eq 0 select [z] else [z,-z]; new:=[];
    for a in ans do for b in choices do Append(~new,Append(a,b)); end for; end for;
    ans:=new;
  end for;
  return ans;
end function;
procedure Bump(~A,k)
  if IsDefined(A,k) then A[k]+:=1; else A[k]:=1; end if;
end procedure;

print "TARGET_22224_DIRECT_CONTACT_BOUNDARY_START",PrimeList;
for p in PrimeList do
  K:=GF(p); outname:=Sprintf("%o_p%o.tsv",output_prefix,p);out:=Open(outname,"w");
  fprintf out,"p\tu\tt\tw\ta\tb\tc\td\tflags\tr0\tr1\tr2\tr3\n";
  triples:=0; split_boundary_triples:={}; cover_boundary_triples:={};
  incidences:=0; sigcounts:=AssociativeArray();
  for u in K do for t in K do for w in K do
    triples+:=1;P,F,q,h:=DirectData(K,u,t,w);X:=P.1;
    ok,mags:=RootMagnitudesWithMultiplicity(F); if not ok then continue; end if;
    split_boundary:=Discriminant(X*F) eq 0 or t eq 0 or
                    Discriminant(q) eq 0 or Degree(GCD(q,X*F)) gt 0;
    if not split_boundary then continue; end if;
    Include(~split_boundary_triples,<Z!u,Z!t,Z!w>);
    for vals in SignChoices(mags) do
      a,b,c,d:=Explode(vals);rads:=CoverRadicands(a,b,c,d);
      if not &and[IsSquare(r):r in rads] then continue; end if;
      flags:=[];
      if t eq 0 then Append(~flags,"t0"); end if;
      if Discriminant(q) eq 0 then Append(~flags,"qrep"); end if;
      if Degree(GCD(q,X*F)) gt 0 then Append(~flags,"qgcd"); end if;
      if &or[z eq 0:z in mags] then Append(~flags,"root0"); end if;
      if #Set(mags) lt 4 then Append(~flags,"rootrep"); end if;
      if &or[r eq 0:r in rads] then Append(~flags,"rad0"); end if;
      sig:=#flags eq 0 select "other" else Join(flags,"+"); Bump(~sigcounts,sig);
      Include(~cover_boundary_triples,<Z!u,Z!t,Z!w>);incidences+:=1;
      fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
              p,Z!u,Z!t,Z!w,Z!a,Z!b,Z!c,Z!d,sig,
              Z!rads[1],Z!rads[2],Z!rads[3],Z!rads[4];
    end for;
  end for;end for;end for;delete out;
  print "DIRECT_CONTACT_BOUNDARY","p",p,"triples",triples,
        "split_boundary_triples",#split_boundary_triples,
        "cover_boundary_triples",#cover_boundary_triples,
        "signed_incidences",incidences,
        "signatures",Sort([<k,sigcounts[k]>:k in Keys(sigcounts)]),
        "file",outname;
end for;
print "TARGET_22224_DIRECT_CONTACT_BOUNDARY_DONE";
UnsetLogFile();quit;
