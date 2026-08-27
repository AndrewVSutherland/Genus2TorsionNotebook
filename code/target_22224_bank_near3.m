// Rank the existing full A(2,2,2,8) bank by how long its Jacobian orders
// retain a factor 3 at successive good primes.  These are deformation seeds,
// not candidate solutions: the first good prime with 3 not dividing #J is a
// rigorous obstruction.
SetColumns(0);
if not assigned bank_file then bank_file:="data/tor2228_bank.txt"; end if;
if not assigned top_n then top_n:=40; end if;
if Type(top_n) eq MonStgElt then top_n:=StringToInteger(top_n); end if;
if not assigned exclude_rectangles then exclude_rectangles:=true; end if;
if Type(exclude_rectangles) eq MonStgElt then
    exclude_rectangles:=exclude_rectangles in {"true","True","1","yes"};
end if;

Z:=Integers(); Q:=Rationals(); P<x>:=PolynomialRing(Q);
primes:=[11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,
         73,79,83,89,97,101,103,107,109,113,127,131,137,139,
         149,151,157,163,167,173,179,181,191,193,197,199];

function ReadRows(filename)
    out:=[]; seen:={};
    for raw in Split(Read(filename),"\n") do
        s:=raw;
        if #s gt 0 and s[#s] eq "\r" then s:=s[1..#s-1]; end if;
        if #s lt 2 or s[1] ne "[" or s[#s] ne "]" then continue; end if;
        parts:=Split(s[2..#s-1],",");
        if #parts ne 4 then continue; end if;
        v:=[StringToInteger(z):z in parts];
        g:=GCD(v); v:=[z div g:z in v];
        key:=Sprint(v);
        if key in seen then continue; end if;
        Include(~seen,key); Append(~out,v);
    end for;
    return out;
end function;

function Score(v)
    f:=x*&*[x+(Q!z)^2:z in v];
    used:=[]; good3:=0; kill:=0; gcdJ:=0;
    for p in primes do
        fp:=ChangeRing(f,GF(p));
        if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue; end if;
        n:=Z!Evaluate(LPolynomial(HyperellipticCurve(fp)),1);
        Append(~used,<p,n>);
        gcdJ:=gcdJ eq 0 select n else GCD(gcdJ,n);
        if n mod 3 ne 0 then kill:=p; break; end if;
        good3+:=1;
    end for;
    return good3,kill,gcdJ,used;
end function;

function IsMultiplicativeRectangle(v)
    return v[1]*v[2] eq v[3]*v[4]
        or v[1]*v[3] eq v[2]*v[4]
        or v[1]*v[4] eq v[2]*v[3];
end function;

rows:=ReadRows(bank_file); ranked:=[];
for i in [1..#rows] do
    if exclude_rectangles and IsMultiplicativeRectangle(rows[i]) then
        continue;
    end if;
    depth,kill,g,used:=Score(rows[i]);
    Append(~ranked,<depth,kill,rows[i],g,used>);
end for;
Sort(~ranked,func<a,b|
    a[1] gt b[1] select -1 else
    a[1] lt b[1] select 1 else
    Sprint(a[3]) lt Sprint(b[3]) select -1 else 1>);
print "TARGET_22224_BANK_NEAR3","unique_bank",#rows,
      "ranked_after_rectangle_filter",#ranked,
      "exclude_rectangles",exclude_rectangles,"top",Min(top_n,#ranked);
for i in [1..Min(top_n,#ranked)] do
    r:=ranked[i];
    print "NEAR3","rank",i,"depth",r[1],"kill",r[2],"tuple",r[3],
          "gcd",r[4],"used",r[5];
end for;
