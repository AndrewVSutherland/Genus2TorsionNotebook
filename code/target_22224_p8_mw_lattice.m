//////////////////////////////////////////////////////////////////////
// Rank-2 Mordell--Weil lattice walk on all three P8 elliptic quotients.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(8);
if not assigned Bound then Bound:=20; end if;
if Type(Bound) eq MonStgElt then Bound:=StringToInteger(Bound); end if;
if not assigned log_file then log_file:="results/target_22224_p8_mw_lattice_B20.log"; end if;
if not assigned output_file then output_file:="results/target_22224_p8_mw_lattice_B20.tsv"; end if;
SetLogFile(log_file : Overwrite:=true);
Q:=Rationals(); Z:=Integers(); R<T>:=PolynomialRing(Q);
fixed:=[Q!528,-726,-891]; d0:=Q!50;
out:=Open(output_file,"w");
fprintf out,"pair\tm\tn\ttorsion_translate\tt\ta\tb\tc\td\tprimitive_a\tprimitive_b\tprimitive_c\tprimitive_d\tthree_survivor\tgcd_bound\n";

function ExactSquare(x)
    x:=Q!x;if x lt 0 then return false,Q!0; end if;
    a,ra:=IsSquare(Z!Numerator(x));b,rb:=IsSquare(Z!Denominator(x));
    return a and b,(a and b) select Q!ra/rb else Q!0;
end function;
function Primitive(vals)
    den:=LCM([Denominator(z):z in vals]); v:=[Z!(z*den):z in vals];
    g:=GCD(v);if g ne 0 then v:=[z div g:z in v];end if;
    for z in v do if z ne 0 then if z lt 0 then v:=[-x:x in v];end if;break;end if;end for;
    return v;
end function;
function FullAt(t)
    vals:=fixed cat [d0*t^2];
    if #Set([z^2:z in vals]) ne 4 then return false,vals;end if;
    for z in fixed do ok,r:=ExactSquare((z+d0*t^2)/(z+d0));if not ok then return false,vals;end if;end for;
    return true,vals;
end function;
prime_list:=[11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,149,151,157];
function ThreeBound(vals)
    f:=T*&*[T+z^2:z in vals];g:=0;used:=[];
    for p in prime_list do
        try
            fp:=ChangeRing(f,GF(p));if Degree(fp) ne 5 or Discriminant(fp) eq 0 then continue;end if;
            nn:=Z!Evaluate(LPolynomial(HyperellipticCurve(fp)),1);
            g:=(g eq 0) select nn else GCD(g,nn);Append(~used,<p,nn,g>);
            if g mod 3 ne 0 then return false,g,used;end if;
        catch e
            continue;
        end try;
    end for;
    return g ne 0 and g mod 3 eq 0,g,used;
end function;

print "P8_MW_LATTICE_START","Bound",Bound,"fixed",fixed,"d0",d0;
seenT:={Q|};full:=0;newfull:=0;three:=0;tested:=0;
for ij in [<1,2>,<1,3>,<2,3>] do
    i,j:=Explode(ij);x:=fixed[i];y:=fixed[j];
    Rx:=(x+d0*T^2)/(x+d0);Ry:=(y+d0*T^2)/(y+d0);
    C:=HyperellipticCurve(Rx*Ry);P0:=C![Q!1,Q!1,Q!1];
    E,phi:=EllipticCurve(C,P0);Einv:=Inverse(phi);
    Emin,minmap:=MinimalModel(E);mininv:=Inverse(minmap);
    gg:=Generators(Emin);assert #gg eq 3 and Order(gg[1]) eq 2;
    tor:=gg[1];P:=gg[2];Qp:=gg[3];
    print "P8_MW_QUOTIENT",ij,"minimal",Emin,"generators",gg,
          "height_pairing",HeightPairingMatrix([P,Qp]);
    for m in [-Bound..Bound] do
      for n in [-Bound..Bound] do
        for e in [0..1] do
          tested+:=1;epmin:=m*P+n*Qp+e*tor;
          try cp:=Einv(mininv(epmin)); catch err continue;end try;
          if cp[3] eq 0 then continue;end if;
          tt:=Q!(cp[1]/cp[3]);if tt in seenT then continue;end if;Include(~seenT,tt);
          good,vals:=FullAt(tt);if not good then continue;end if;
          full+:=1;known:=tt in {Q!1,Q!-1,Q!19/5,Q!-19/5};if not known then newfull+:=1;end if;
          survives,g,used:=ThreeBound(vals);if survives then three+:=1;end if;
          prim:=Primitive(vals);
          fprintf out,"%o-%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
                  i,j,m,n,e,tt,vals[1],vals[2],vals[3],vals[4],
                  prim[1],prim[2],prim[3],prim[4],survives select 1 else 0,g;
          print "P8_MW_FULL","pair",ij,"coeffs",m,n,e,"t",tt,"known",known,
                "primitive",prim,"three_survivor",survives,"gcd",g,"used",used;
          if survives then
              f:=T*&*[T+z^2:z in vals];G,mp:=TorsionSubgroup(Jacobian(HyperellipticCurve(f)));
              print "P8_MW_EXACT_TORSION",Invariants(G);
          end if;
        end for;
      end for;
    end for;
end for;
delete out;
print "P8_MW_LATTICE_DONE","tested",tested,"distinct_t",#seenT,"full",full,
      "new_full",newfull,"three_survivors",three,"output",output_file;
UnsetLogFile();quit;
