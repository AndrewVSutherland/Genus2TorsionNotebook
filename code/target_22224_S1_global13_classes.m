//////////////////////////////////////////////////////////////////////
// Global 13-adic coefficient-class pass for the rank-one S1 fibre
//
//   (-4,9,30,(-2166/245)*T^2),
//   E: y^2=x^3+x^2+3392*x+62288,
//   E(Q)=<G=(152,-2028)> + <W=(-17,0)>.
//
// The reduction at 13 is nonsplit I_4.  Direct exact computation gives
// z(14G)=-x(14G)/y(14G) in 13^2 Z_13^x; hence the order of G modulo
// the depth-d formal subgroup E_d is 14*13^(d-2), d>=2.
//
// We enumerate both torsion cosets modulo that exact coefficient period
// and compare the projective S1 coefficient tuple with the corrected
// sampled direct-contact bank modulo 13^5.  The T-map is evaluated in
// homogeneous form, avoiding division at its pole/base point.
//////////////////////////////////////////////////////////////////////
SetColumns(0); SetSeed(181305);
Q:=Rationals(); Z:=Integers(); p:=13; mod5:=p^5;
if not assigned FormalDepth then FormalDepth:=5;
elif Type(FormalDepth) eq MonStgElt then FormalDepth:=StringToInteger(FormalDepth); end if;
if not assigned Precision then Precision:=80;
elif Type(Precision) eq MonStgElt then Precision:=StringToInteger(Precision); end if;
if not assigned BlockSize then BlockSize:=1;
elif Type(BlockSize) eq MonStgElt then BlockSize:=StringToInteger(BlockSize); end if;
if not assigned log_file then
    log_file:=Sprintf("results/target_22224_S1_global13_classes_depth%o.log",FormalDepth);
end if;
if not assigned output_file then
    output_file:=Sprintf("results/target_22224_S1_global13_classes_depth%o.tsv",FormalDepth);
end if;
SetLogFile(log_file:Overwrite:=true);

E:=EllipticCurve([Q!0,Q!1,Q!0,Q!3392,Q!62288]);
G:=E![Q!152,Q!-2028,Q!1]; W:=E![Q!-17,Q!0,Q!1];
tors:=[E!0,W];
fixed:=[Q!-4,Q!9,Q!30]; d0:=Q!-2166/245;

// After clearing the common denominator from the Magma quartic map,
// T=N/D with these two integral linear forms on the plane cubic.
nvec:=[Z!-89425,Z!1463,Z!4228700];
dvec:=[Z!-86279,Z!1463,Z!-7215668];

function AbsRes(a,m)
    a:=a mod m; return Min(a,(m-a) mod m);
end function;

function DeepKey(vals,p,m)
    best:=<m,m,m,m>;
    for i in [1..4] do
        x:=vals[i] mod m; if x mod p eq 0 then continue; end if;
        q:=InverseMod(x,m);
        zz:=Sort([AbsRes((z mod m)*q,m):z in vals]);
        key:=<zz[1],zz[2],zz[3],zz[4]>;
        if key lt best then best:=key; end if;
    end for;
    return best;
end function;

function ReadBank(path,p,m)
    lines:=Split(Read(path),"\n"); dummy:=<Z!m,Z!m,Z!m,Z!m>;
    out:={dummy}; Exclude(~out,dummy);
    for i in [2..#lines] do
        if #lines[i] eq 0 then continue; end if;
        z:=Split(lines[i],"\t"); if #z lt 8 then continue; end if;
        Include(~out,DeepKey([StringToInteger(z[j]):j in [5..8]],p,m));
    end for;
    return out;
end function;

bank_path:="results/target_22224_direct_contact_deep13_padic_tangent_p13.tsv";
bank:=ReadBank(bank_path,p,mod5);

K:=pAdicField(p,Precision); EK:=BaseChange(E,K);
function KP(P)
    if P eq Parent(P)!0 then return EK!0; end if;
    return EK![K!P[1],K!P[2],K!P[3]];
end function;
GK:=KP(G); torsK:=[KP(T):T in tors];

function RatMod(q,m)
    q:=Q!q; return (Numerator(q) mod m)*InverseMod(Denominator(q) mod m,m) mod m;
end function;

// The displayed N/D formula has one removable base point B=3G+W.
// Both lines meet E simply there.  Cancelling the common local parameter
// gives the exact value from the tangent derivatives.
BQ:=3*G+W; xb:=Q!BQ[1]; yb:=Q!BQ[2];
slope:=(3*xb^2+2*xb+Q!3392)/(2*yb);
tbase:=(Q!nvec[1]+Q!nvec[2]*slope)/(Q!dvec[1]+Q!dvec[2]*slope);
baseVals:=[RatMod(fixed[i],mod5):i in [1..3]] cat [RatMod(d0*tbase^2,mod5)];
baseKey:=DeepKey(baseVals,p,mod5);

function LinEval(v,P)
    // The identity is (X:Y:Z)=(0:1:0), not an affine Z=1 point.
    if P eq EK!0 then return K!v[2]; end if;
    return K!v[1]*P[1]+K!v[2]*P[2]+K!v[3]*P[3];
end function;

// Return the canonical coefficient key.  At the unique simultaneous zero
// N=D=0 (the rational point 3G+W), D has the extra local zero and T has a
// pole, so the limiting projective coefficient tuple is [0,0,0,1].
function PointKey(P)
    nn:=LinEval(nvec,P); dd:=LinEval(dvec,P);
    vn:=nn eq 0 select Precision else Valuation(nn);
    vd:=dd eq 0 select Precision else Valuation(dd);
    if vn ge Precision div 2 and vd ge Precision div 2 then
        return baseKey,vn,vd,true;
    end if;
    valsK:=[K!fixed[i]*dd^2:i in [1..3]] cat [K!d0*nn^2];
    nzvals:=[Valuation(x):x in valsK|x ne 0];
    minv:=Minimum(nzvals); scale:=K!(p^(-minv));
    valsZ:=[x eq 0 select Z!0 else (Z!(scale*x)) mod mod5:x in valsK];
    return DeepKey(valsZ,p,mod5),
           vn,vd,false;
end function;

assert FormalDepth ge 2;
period:=14*p^(FormalDepth-2);
H:=period*GK;
hz:=H eq EK!0 select K!0 else -H[1]/H[2];
assert H ne EK!0 and Valuation(hz) eq FormalDepth;
Hprev:=(period div p)*GK;
assert FormalDepth eq 2 or Valuation(-Hprev[1]/Hprev[2]) eq FormalDepth-1;
H0:=14*GK; assert Valuation(-H0[1]/H0[2]) eq 2;
for n in [1..13] do
    PP:=n*GK;
    assert PP ne EK!0 and not (Valuation(PP[1]) le -2 and Valuation(PP[2]) le -3);
end for;

out:=Open(output_file,"w");
fprintf out,"m_residue\ttorsion_coset\tkey\tv13_N\tv13_D\tbasepoint\n";
tested:=0; hits:=0; basepoints:=0; unique:={};
vnpairs:=AssociativeArray();
print "S1_GLOBAL13_CLASSES_START","depth",FormalDepth,"period",period,
      "classes",2*period,"formal_step_v",Valuation(hz),
      "local_information",LocalInformation(E,p),"bank_keys",#bank,
      "precision",Precision,"block_size",BlockSize,
      "base_point",BQ,"base_t",tbase,"base_key",baseKey,
      "base_key_in_bank",baseKey in bank;
parentPeriod:=FormalDepth ge 3 select period div p else period;
parentKeys:=[* [* *] : ti in [1..#torsK] *];
parentMismatch:=0;
for ti in [1..#torsK] do
  for blockStart in [0..period-1 by BlockSize] do
    blockEnd:=Min(period-1,blockStart+BlockSize-1);
    // Re-anchor each short block by binary scalar multiplication.  A long
    // repeated-addition walk loses 13-adic precision on this bad-reduction
    // model and is not acceptable for a residue certificate.
    P:=blockStart*GK+torsK[ti];
    for m in [blockStart..blockEnd] do
        key,vn,vd,isbase:=PointKey(P); tested+:=1; Include(~unique,key);
        if m lt parentPeriod then
            Append(~parentKeys[ti],key);
        elif FormalDepth ge 6 and key ne parentKeys[ti][(m mod parentPeriod)+1] then
            parentMismatch+:=1;
            if parentMismatch le 20 then
                print "S1_GLOBAL13_PARENT_MISMATCH","m",m,"ti",ti,
                      "key",key,"parent_key",parentKeys[ti][(m mod parentPeriod)+1];
            end if;
        end if;
        tag:=Sprintf("%o,%o",vn,vd);
        if not IsDefined(vnpairs,tag) then vnpairs[tag]:=0; end if;
        vnpairs[tag]+:=1;
        if isbase then basepoints+:=1; end if;
        if key in bank then
            hits+:=1;
            fprintf out,"%o\t%o\t%o\t%o\t%o\t%o\n",m,ti,key,vn,vd,isbase select 1 else 0;
            print "S1_GLOBAL13_BANK_HIT","m",m,"ti",ti,"key",key,
                  "vN",vn,"vD",vd,"base",isbase;
        end if;
        P+:=GK;
    end for;
  end for;
end for;
delete out;
print "S1_GLOBAL13_VALUATION_DISTRIBUTION";
for tag in Sort(Setseq(Keys(vnpairs))) do print tag,vnpairs[tag]; end for;
print "S1_GLOBAL13_CLASSES_DONE","tested",tested,"unique_keys",#unique,
      "basepoints",basepoints,"bank_hits",hits,
      "parent_period",parentPeriod,"parent_mismatches",parentMismatch,
      "output",output_file;
UnsetLogFile(); quit;
