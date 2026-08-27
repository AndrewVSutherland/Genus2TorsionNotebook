//////////////////////////////////////////////////////////////////////
// Elliptic-fibration search on the p=31 external boundary.
//
// Fix A.  The first double-square equation becomes the genus-one quartic
//
//   W^2 = B^4 + (A^2-3) B^2 + A^-2,   W=B*rho,
//
// with the rational boundary point B=1.  We convert it to an elliptic
// curve, search for rational points and take short multiples, then impose
// the second square exactly.  A is restricted by A == +/-1 (mod 31), so
// every fibre searched lies in the cubic-thick external boundary chamber.
//////////////////////////////////////////////////////////////////////

SetColumns(0); SetSeed(1);
if not assigned AHeight then AHeight:=100; end if;
if not assigned PointBound then PointBound:=300; end if;
if not assigned MultipleBound then MultipleBound:=30; end if;
if not assigned BoundaryMode then BoundaryMode:="external"; end if;
if not assigned log_file then
  log_file:="results/target_22224_p31_boundary_fibration.log";
end if;
if not assigned output_file then
  output_file:="results/target_22224_p31_boundary_fibration_candidates.tsv";
end if;
SetLogFile(log_file:Overwrite:=true);

Q:=Rationals(); Z:=Integers(); PR<B>:=PolynomialRing(Q);
if Type(AHeight) eq MonStgElt then AHeight:=StringToInteger(AHeight); end if;
if Type(PointBound) eq MonStgElt then PointBound:=StringToInteger(PointBound); end if;
if Type(MultipleBound) eq MonStgElt then MultipleBound:=StringToInteger(MultipleBound); end if;
out:=Open(output_file,"w");
fprintf out,"A_num\tA_den\tB_num\tB_den\trho_num\trho_den\tsigma_num\tsigma_den\n";

function IsExactSquareQ(x)
  x:=Q!x;
  return x ge 0 and IsSquare(Numerator(x)) and IsSquare(Denominator(x));
end function;

function ExactCollision(A,B)
  C:=1/(A*B);
  vals:=[A^2,B^2,C^2];
  if #Seqset(vals) lt 3 then return true; end if;
  R:=A^2+B^2+C^2-3; S:=1/A^2+1/B^2+1/C^2-3;
  return R eq 0 or S eq 0 or A^2 eq 1 or B^2 eq 1 or C^2 eq 1;
end function;

function Mod31Value(x)
  if Numerator(x) mod 31 eq 0 or Denominator(x) mod 31 eq 0 then
    return false,GF(31)!0;
  end if;
  return true,(GF(31)!Numerator(x))/(GF(31)!Denominator(x));
end function;

function InternalBoundary31(A,B)
  okA,a:=Mod31Value(A); okB,b:=Mod31Value(B);
  if not okA or not okB then return false; end if;
  c:=1/(a*b); x:=a^2; y:=b^2; z:=c^2;
  return x eq y or x eq z or y eq z;
end function;

fibres:=0; singular_fibres:=0; conversion_fail:=0; pointsearch_fail:=0;
searched_E_points:=0; inverse_fail:=0; first_square_points:=0;
second_square_points:=0; mode_boundary_points:=0; nontrivial_points:=0; written:=0;
seenA:={}; seenAB:={};

print "TARGET_22224_P31_BOUNDARY_FIBRATION_START",
      "mode",BoundaryMode,"AHeight",AHeight,
      "PointBound",PointBound,"MultipleBound",MultipleBound;

for den in [1..AHeight] do for num in [1..AHeight] do
  if GCD(num,den) ne 1 or num mod 31 eq 0 or den mod 31 eq 0 then continue; end if;
  A:=Q!num/den;
  ar:=GF(31)!num/(GF(31)!den);
  if BoundaryMode eq "external" and ar ne 1 and ar ne -1 then continue; end if;
  if A in seenA then continue; end if; Include(~seenA,A);
  quartic:=B^4+(A^2-3)*B^2+1/A^2;
  if Discriminant(quartic) eq 0 then singular_fibres+:=1; continue; end if;
  fibres+:=1;
  Cq:=HyperellipticCurve(quartic);
  base:=Cq![1,A-1/A,1];
  try
    E,mp:=EllipticCurve(Cq,base);
  catch e
    conversion_fail+:=1; continue;
  end try;
  eps:={E!0};
  try
    found:=PointSearch(E,PointBound);
    for ep in found do
      Include(~eps,ep);
      if Order(ep) eq 0 then
        for n in [-MultipleBound..MultipleBound] do Include(~eps,n*ep); end for;
      end if;
    end for;
  catch e
    pointsearch_fail+:=1;
  end try;
  searched_E_points +:= #eps;
  for ep in eps do
    try cp:=ep@@mp; catch e inverse_fail+:=1; continue; end try;
    if cp[3] eq 0 then continue; end if;
    b:=Q!(cp[1]/cp[3]); w:=Q!(cp[2]/cp[3]^2);
    if b eq 0 then continue; end if;
    rho:=Q!(w/b);
    if rho^2 ne A^2+b^2+1/(A*b)^2-3 then continue; end if;
    key:=Sprint(<A,b>); if key in seenAB then continue; end if; Include(~seenAB,key);
    first_square_points+:=1;
    S:=Q!(1/A^2+1/b^2+(A*b)^2-3);
    if not IsExactSquareQ(S) then continue; end if;
    oksq,sigma:=IsSquare(S); assert oksq;
    second_square_points+:=1; sigma:=Q!sigma;
    if BoundaryMode eq "internal" and not InternalBoundary31(A,b) then continue; end if;
    mode_boundary_points+:=1;
    if ExactCollision(A,b) then continue; end if;
    nontrivial_points+:=1;
    fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\t%o\t%o\n",
            Numerator(A),Denominator(A),Numerator(b),Denominator(b),
            Numerator(rho),Denominator(rho),Numerator(sigma),Denominator(sigma);
    written+:=1;
    print "NONTRIVIAL_DOUBLE_SQUARE",
          "A",A,"B",b,"C",1/(A*b),"rho",rho,"sigma",sigma;
  end for;
end for; end for;

delete out;
print "P31_BOUNDARY_FIBRATION_SUMMARY",
      "fibres",fibres,"singular",singular_fibres,
      "conversion_fail",conversion_fail,"pointsearch_fail",pointsearch_fail,
      "E_points",searched_E_points,"inverse_fail",inverse_fail,
      "first_square_points",first_square_points,
      "second_square_points",second_square_points,
      "mode_boundary_points",mode_boundary_points,
      "nontrivial",nontrivial_points,"written",written;
print "OUTPUT",output_file;
print "TARGET_22224_P31_BOUNDARY_FIBRATION_DONE";
UnsetLogFile(); quit;
