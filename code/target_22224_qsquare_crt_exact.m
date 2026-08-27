//////////////////////////////////////////////////////////////////////
// Exact verifier for rational survivors of the boundary-aware q-square
// CRT sieve for [2,2,2,24].
//
// Run from the repository root:
//   magma -b code/target_22224_qsquare_crt_exact.m
//
// Options:
//   candidate_file:=results/target_22224_qsquare_crt_sieve_h187_candidates.tsv
//   log_file:=results/target_22224_qsquare_crt_exact_h187.log
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if not assigned candidate_file then
    candidate_file := "results/target_22224_qsquare_crt_sieve_h187_candidates.tsv";
end if;
if not assigned log_file then
    log_file := "results/target_22224_qsquare_crt_exact_h187.log";
end if;
SetLogFile(log_file : Overwrite:=true);

Q:=Rationals(); Z:=Integers();
P<x>:=PolynomialRing(Q);

function MakeMonic(g)
    return g/LeadingCoefficient(g);
end function;

function DirectData(aa,bb,cc,dd,uu,tt,vv)
    f:=x*(x+aa^2)*(x+bb^2)*(x+cc^2)*(x+dd^2);
    q:=x^2+uu*x+tt^2;
    h:=x^3+(1+3*uu)/2*x^2+vv*x+tt^3;
    return f,q,h;
end function;

function MarkedD12Data(aa,bb,cc,dd,uu,tt,vv)
    f,q,h:=DirectData(aa,bb,cc,dd,uu,tt,vv);
    if h^2-f-q^3 ne 0 then return false,f,P!0,P!0; end if;
    g4:=(x-aa*bb)*(x-cc*dd)-x*(aa+bb)*(cc+dd);
    L4:=(x-aa*bb)*(cc+dd)+(aa+bb)*(x-cc*dd);
    if Degree(g4) ne 2 then return false,f,P!0,P!0; end if;
    v3:=h mod q; v4:=(-x*L4) mod g4;
    Rm:=(q*x) mod g4; Rn:=q mod g4; Rt:=(v4-v3) mod g4;
    a0:=Coefficient(Rm,0); a1:=Coefficient(Rm,1);
    b0:=Coefficient(Rn,0); b1:=Coefficient(Rn,1);
    c0:=Coefficient(Rt,0); c1:=Coefficient(Rt,1);
    det:=a0*b1-a1*b0;
    if det eq 0 then return false,f,P!0,P!0; end if;
    mc:=(c0*b1-c1*b0)/det; nc:=(a0*c1-a1*c0)/det;
    ell:=v3+q*(mc*x+nc);
    if mc eq 0 or (ell^2-f) mod (q*g4) ne 0 then
        return false,f,P!0,P!0;
    end if;
    Uraw:=ExactQuotient(ell^2-f,q*g4);
    if Degree(Uraw) ne 2 then return false,f,P!0,P!0; end if;
    U:=MakeMonic(Uraw); V:=(-ell) mod U;
    return (V^2-f) mod U eq 0,f,U,V;
end function;

function IntegralSquareModel(f)
    den:=LCM([Denominator(Coefficient(f,i)):i in [0..Degree(f)]]);
    return P!(den^2*f),den;
end function;

function RootPowerWitness(f)
    C:=HyperellipticCurve(f);
    for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97] do
        try
            fp:=ChangeRing(f,GF(p));
            if Degree(fp) notin {5,6} or Discriminant(fp) eq 0 then continue; end if;
            chi:=LPolynomial(ChangeRing(C,GF(p)));
            if not IsIrreducible(chi) then continue; end if;
            R<T>:=PolynomialRing(Q); chQ:=R!chi;
            K<pi>:=NumberField(chQ);
            degs:=[Degree(MinimalPolynomial(pi^n)):n in [2..12]];
            if &and[d eq 4:d in degs] then return true,p,chQ,degs; end if;
        catch e
            continue;
        end try;
    end for;
    R<T>:=PolynomialRing(Q);
    return false,0,R!0,[];
end function;

function ParseRational(num,den)
    return Q!StringToInteger(num)/StringToInteger(den);
end function;

lines:=Split(Read(candidate_file),"\n");
candidate_rows:=0; smooth_presentations:=0; exact_tests:=0; hits:=0;
seen:={};

print "TARGET_22224_QSQUARE_CRT_EXACT_START";
print "candidate_file",candidate_file,"data_lines",Maximum(0,#lines-1);

for li in [2..#lines] do
    if #lines[li] eq 0 then continue; end if;
    fields:=Split(lines[li],"\t");
    if #fields lt 10 then continue; end if;
    candidate_rows+:=1;
    A0:=ParseRational(fields[1],fields[2]);
    B0:=ParseRational(fields[3],fields[4]);
    C0:=ParseRational(fields[5],fields[6]);
    rho:=ParseRational(fields[7],fields[8]);
    sigma:=ParseRational(fields[9],fields[10]);
    assert A0*B0*C0 eq 1;
    assert A0^2+B0^2+C0^2-3 eq rho^2;
    assert 1/A0^2+1/B0^2+1/C0^2-3 eq sigma^2;

    ss:=1/(2*rho); tt:=ss^2; uu:=2*tt;
    mags:=[ss*A0,ss*B0,ss*C0]; ddmag:=2*ss^2*sigma;
    vv:=3*tt^2+ddmag^2/2;
    choices:=[<mags[1],mags[2],mags[3]>,
             <mags[1],mags[3],mags[2]>,
             <mags[2],mags[3],mags[1]>];
    for pairing in [1..3] do
      ch:=choices[pairing];
      for sa in [-1,1] do for sb in [-1,1] do for sc in [-1,1] do
        for sd in [-1,1] do
          aa:=sa*ch[1]; bb:=sb*ch[2]; cc:=sc*ch[3]; dd:=sd*ddmag;
          ok,f,U,V:=MarkedD12Data(aa,bb,cc,dd,uu,tt,vv);
          if not ok or Discriminant(f) eq 0 then continue; end if;
          key:=Sprint(<f,U,V>);
          if key in seen then continue; end if;
          Include(~seen,key); smooth_presentations+:=1;
          fI,scale:=IntegralSquareModel(f);
          J:=Jacobian(HyperellipticCurve(fI));
          D:=J![U,scale*V];
          ord:=Order(D); exact_tests+:=1;
          if ord ne 12 then
              print "EXACT_BAD_MARKED_ORDER",ord,A0,B0,C0,pairing,sa,sb,sc,sd;
              continue;
          end if;
          divisible,half:=IsDivisibleBy(D,2);
          if not divisible then continue; end if;
          T,mp:=TorsionSubgroup(J);
          inv:=Invariants(T);
          simple,pcert,chi,degs:=RootPowerWitness(fI);
          hits+:=1;
          print "TARGET_22224_EXACT_HIT",
                "A",A0,"B",B0,"C",C0,"rho",rho,"sigma",sigma,
                "pairing",pairing,"signs",[sa,sb,sc,sd],
                "a",aa,"b",bb,"c",cc,"d",dd,"u",uu,"t",tt,"v",vv,
                "marked_order",ord,"torsion",inv,
                "simple",simple,"pcert",pcert,"chi",chi,"power_degrees",degs,
                "curve",fI;
        end for;
      end for; end for; end for;
    end for;
end for;

print "TARGET_22224_QSQUARE_CRT_EXACT_DONE",
      "candidate_rows",candidate_rows,
      "unique_smooth_presentations",smooth_presentations,
      "exact_tests",exact_tests,"hits",hits;
UnsetLogFile();
quit;

