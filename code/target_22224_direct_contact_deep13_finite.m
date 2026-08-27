//////////////////////////////////////////////////////////////////////
// Fresh finite masks for the direct cubic-contact [2,2,2,24] search.
//
// IMPORTANT: this derives F directly from h^2-q^3 and does not use the
// superseded Aaux equations or their p=13 boundary masks.
//
// q = x^2+u*x+t^2,
// h = x^3+(1+3u)/2*x^2+w*x+t^3,
// h^2-q^3 = x*F.
//
// A target row must pass all of:
//   * F has four distinct roots -a_i^2;
//   * the exact four-radicand A(2,2,2,8) cover for some signs;
//   * D3=[q,h mod q] has exact order 3 in J(F_p);
//   * the corresponding marked D4 has order 4 and belongs to 2J(F_p).
//////////////////////////////////////////////////////////////////////

SetColumns(0); SetSeed(1);
if not assigned PrimeList then PrimeList:=[11,13,17,19,23,29,31];
elif Type(PrimeList) eq MonStgElt then
  PrimeList:=[StringToInteger(s):s in Split(PrimeList,",")];
end if;
if not assigned output_prefix then
  output_prefix:="results/target_22224_direct_contact_deep13";
end if;
if not assigned log_file then log_file:=output_prefix cat "_finite.log"; end if;
SetLogFile(log_file:Overwrite:=true);
Z:=Integers();

function DirectData(K,u,t,w)
  P<X>:=PolynomialRing(K);
  q:=X^2+u*X+t^2;
  h:=X^3+(1+3*u)/2*X^2+w*X+t^3;
  F:=ExactQuotient(h^2-q^3,X);
  return P,F,q,h;
end function;

function SquareRootBranches(F)
  fac:=Factorization(F);
  if #fac ne 4 or not &and[Degree(z[1]) eq 1 and z[2] eq 1:z in fac] then
    return false,[];
  end if;
  vals:=[]; roots:=[];
  for z in fac do
    rr:=-Coefficient(z[1],0)/Coefficient(z[1],1);
    if rr eq 0 or rr in roots or not IsSquare(-rr) then return false,[]; end if;
    Append(~roots,rr); Append(~vals,SquareRoot(-rr));
  end for;
  return true,vals;
end function;

function CoverRadicands(a,b,c,d)
  return [a*b*c*d,
          a*(a+b)*(a+c)*(a+d),
          b*(b+a)*(b+c)*(b+d),
          c*(c+a)*(c+b)*(c+d)];
end function;

function FiniteDivisibleBy2(J,D)
  G,mp:=AbelianGroup(J); z:=D@@mp;
  cs:=Eltseq(z); invs:=Invariants(G);
  return &and[GCD(2,invs[i]) eq 1 or cs[i] mod 2 eq 0:i in [1..#cs]];
end function;

print "TARGET_22224_DIRECT_CONTACT_FINITE_START",PrimeList;
for p in PrimeList do
  K:=GF(p); outname:=Sprintf("%o_p%o.tsv",output_prefix,p);
  out:=Open(outname,"w");
  fprintf out,"p\tu\tt\tw\tpairing\ta\tb\tc\td\tordD3\tordD4\tD4div2\n";
  triples:=0; tzero:=0; repeatedq:=0; singularF:=0; split_square:=0;
  cover_signs:=0; jac_tests:=0; d3order3:=0; d4order4:=0;
  group_divisible:=0; target_triples:={}; target_curves:={};
  bad_identity:=0;
  for u in K do for t in K do for w in K do
    triples+:=1;
    if t eq 0 then tzero+:=1; continue; end if;
    P,F,q,h:=DirectData(K,u,t,w); X:=P.1;
    if h^2-q^3 ne X*F then bad_identity+:=1; continue; end if;
    if Discriminant(q) eq 0 then repeatedq+:=1; continue; end if;
    if Discriminant(X*F) eq 0 or Degree(GCD(q,X*F)) gt 0 then
      singularF+:=1; continue;
    end if;
    ok,mags:=SquareRootBranches(F); if not ok then continue; end if;
    split_square+:=1;
    found_target:=false;
    for mask in [0..15] do
      vals:=[((mask div 2^(i-1)) mod 2 eq 0) select mags[i] else -mags[i]
            :i in [1..4]];
      a0,b0,c0,d0:=Explode(vals);
      arrangements:=[<a0,b0,c0,d0>,<a0,c0,b0,d0>,<a0,d0,b0,c0>];
      for pairing in [1..3] do
        ar:=arrangements[pairing]; a:=ar[1];b:=ar[2];c:=ar[3];d:=ar[4];
        rads:=CoverRadicands(a,b,c,d);
        if not &and[r ne 0 and IsSquare(r):r in rads] then continue; end if;
        cover_signs+:=1;
        f:=X*F; J:=Jacobian(HyperellipticCurve(f));
        D3:=J![q,h mod q]; ord3:=Order(D3); jac_tests+:=1;
        if ord3 ne 3 then continue; end if; d3order3+:=1;
        g4:=(X-a*b)*(X-c*d)-X*(a+b)*(c+d);
        L4:=(X-a*b)*(c+d)+(a+b)*(X-c*d);
        if Degree(g4) ne 2 then continue; end if;
        D4:=J![g4/LeadingCoefficient(g4),(-X*L4) mod g4];
        ord4:=Order(D4); if ord4 ne 4 then continue; end if; d4order4+:=1;
        div2:=FiniteDivisibleBy2(J,D4); if not div2 then continue; end if;
        group_divisible+:=1; found_target:=true;
        Include(~target_curves,Sprint(F));
        fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t1\n",
                p,Z!u,Z!t,Z!w,pairing,Z!a,Z!b,Z!c,Z!d,ord3,ord4;
      end for;
    end for;
    if found_target then Include(~target_triples,<Z!u,Z!t,Z!w>); end if;
  end for; end for; end for;
  delete out;
  print "DIRECT_CONTACT_FINITE","p",p,"triples",triples,
        "tzero",tzero,"repeatedq",repeatedq,"singular_or_gcd",singularF,
        "split_negative_squares",split_square,"cover_signs",cover_signs,
        "jac_tests",jac_tests,"D3_order3",d3order3,"D4_order4",d4order4,
        "D4_divisible",group_divisible,"target_triples",#target_triples,
        "target_curves",#target_curves,"bad_identity",bad_identity,
        "file",outname;
end for;
print "TARGET_22224_DIRECT_CONTACT_FINITE_DONE";
UnsetLogFile(); quit;
