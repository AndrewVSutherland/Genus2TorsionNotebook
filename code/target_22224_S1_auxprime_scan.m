//////////////////////////////////////////////////////////////////////
// Search for a prime giving an intrinsic full-cover/3-primary obstruction
// on every residue disk of S1.
//////////////////////////////////////////////////////////////////////
SetColumns(0); Z:=Integers();
fixed:=[Z!-980,Z!2205,Z!7350]; d0:=Z!-2166;
if not assigned Bound then Bound:=499;
elif Type(Bound) eq MonStgElt then Bound:=StringToInteger(Bound); end if;
log_file:="results/target_22224_S1_auxprime_scan.log";
output_file:="results/target_22224_S1_auxprime_scan.tsv";
SetLogFile(log_file:Overwrite:=true);

function Radicands(vals)
    a,b,c,d:=Explode(vals);
    return [a*b*c*d,
            a*(a+b)*(a+c)*(a+d),
            b*(b+a)*(b+c)*(b+d),
            c*(c+a)*(c+b)*(c+d)];
end function;

out:=Open(output_file,"w");
fprintf out,"p\tcover_residues\tsingular_cover\tJ3_cover\tinfinity_full\tcomplete_kill\n";
hits:=[];
for p in PrimesInInterval(11,Bound) do
    if p in {13,19} or 245 mod p eq 0 then continue; end if;
    F:=GF(p); P<X>:=PolynomialRing(F); A,B,C:=Explode([F!z:z in fixed]); D:=F!d0;
    cover:=[]; singular:=[]; j3:=[];
    for tt in F do
        vals:=[A,B,C,D*tt^2]; rr:=Radicands(vals);
        if not &and[IsSquare(r):r in rr] then continue; end if;
        Append(~cover,Z!tt);
        f:=X*&*[X+z^2:z in vals];
        if Degree(f) ne 5 or Discriminant(f) eq 0 then
            Append(~singular,Z!tt); continue;
        end if;
        nj:=Z!Evaluate(LPolynomial(HyperellipticCurve(f)),1);
        if nj mod 3 eq 0 then Append(~j3,<Z!tt,nj>); end if;
    end for;
    inf:=Radicands([A*0,B*0,C*0,D]);
    // Remove the common z^6 before z=0.
    inf:=[A*B*C*D,
          A*(A+B)*(A+C)*D,
          B*(B+A)*(B+C)*D,
          C*(C+A)*(C+B)*D];
    infFull:=&and[IsSquare(r):r in inf];
    kill:=#singular eq 0 and #j3 eq 0 and not infFull;
    fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\n",p,cover,singular,
            [x[1]:x in j3],infFull select 1 else 0,kill select 1 else 0;
    print "S1_AUXPRIME",p,"cover",cover,"singular",singular,
          "j3",j3,"infinity_units",[Z!r:r in inf],
          "infinity_full",infFull,"complete_kill",kill;
    if kill then Append(~hits,p); end if;
end for;
delete out;
print "S1_AUXPRIME_SCAN_DONE","bound",Bound,"complete_kill_primes",hits,
      "output",output_file;
UnsetLogFile(); quit;
