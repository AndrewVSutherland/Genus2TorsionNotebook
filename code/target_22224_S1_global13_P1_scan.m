//////////////////////////////////////////////////////////////////////
// Exhaustive P^1(Z/13^5) scan for the S1 coefficient line.
//
// This is stronger than a Mordell--Weil coefficient scan relative to the
// current 144-row sampled tangent bank: it tests every possible 13-adic
// value of T, including the nonintegral/infinity chart.
//////////////////////////////////////////////////////////////////////
SetColumns(0);
Z:=Integers(); Q:=Rationals(); p:=13; M:=p^5;
log_file:="results/target_22224_S1_global13_P1_scan.log";
output_file:="results/target_22224_S1_global13_P1_scan.tsv";
SetLogFile(log_file:Overwrite:=true);

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
function RatMod(q,m)
    q:=Q!q; return (Numerator(q) mod m)*InverseMod(Denominator(q) mod m,m) mod m;
end function;

bank_path:="results/target_22224_direct_contact_deep13_padic_tangent_p13.tsv";
bank:=ReadBank(bank_path,p,M);
a:=Z!-4; b:=Z!9; c:=Z!30; d0:=RatMod(Q!-2166/245,M);
out:=Open(output_file,"w");
fprintf out,"chart\tparameter\tkey\n";
finite:=0; infinity:=0; hits:=0; uniq:={};
print "S1_GLOBAL13_P1_SCAN_START","modulus",M,"bank_keys",#bank,
      "finite_classes",M,"infinity_chart_classes",p^4;

// Diagnose the first depth at which the S1 line separates from all bank
// rows.  This is already depth two; the depth-five scan below is retained
// as a direct same-modulus comparison.
for e in [1..4] do
    Me:=p^e; banke:=ReadBank(bank_path,p,Me); d0e:=RatMod(Q!-2166/245,Me);
    he:=0; uke:={};
    for t in [0..Me-1] do
        key:=DeepKey([a,b,c,d0e*t^2],p,Me);
        if key in banke then he+:=1; Include(~uke,key); end if;
    end for;
    for u in [0..p^(e-1)-1] do
        s:=p*u; key:=DeepKey([a*s^2,b*s^2,c*s^2,d0e],p,Me);
        if key in banke then he+:=1; Include(~uke,key); end if;
    end for;
    print "S1_GLOBAL13_P1_LEVEL","exponent",e,"modulus",Me,
          "bank_keys",#banke,"projective_classes",Me+p^(e-1),
          "hit_classes",he,"hit_keys",#uke;
end for;

// D is a unit: [N:D]=[t:1], t modulo p^5.
for t in [0..M-1] do
    key:=DeepKey([a,b,c,d0*t^2],p,M); finite+:=1; Include(~uniq,key);
    if key in bank then
        hits+:=1; fprintf out,"finite\t%o\t%o\n",t,key;
        if hits le 20 then print "S1_GLOBAL13_P1_HIT","finite",t,key; end if;
    end if;
end for;

// N is a unit and D is not: [N:D]=[1:13*u], u modulo p^4.
// Multiplying the coefficient tuple by D^2 gives the integral chart
// [a(13u)^2,b(13u)^2,c(13u)^2,d0].
for u in [0..p^4-1] do
    s:=p*u; key:=DeepKey([a*s^2,b*s^2,c*s^2,d0],p,M);
    infinity+:=1; Include(~uniq,key);
    if key in bank then
        hits+:=1; fprintf out,"infinity\t%o\t%o\n",u,key;
        if hits le 20 then print "S1_GLOBAL13_P1_HIT","infinity",u,key; end if;
    end if;
end for;
delete out;
print "S1_GLOBAL13_P1_SCAN_DONE","tested",finite+infinity,
      "finite",finite,"infinity",infinity,"unique_keys",#uniq,
      "bank_hits",hits,"output",output_file;
UnsetLogFile(); quit;
