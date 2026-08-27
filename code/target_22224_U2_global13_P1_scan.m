//////////////////////////////////////////////////////////////////////
// Exhaustive P^1(Z/13^e), 1<=e<=5, comparison for the reciprocal
// U2/U1 repeated-fibre coefficient lines against the current corrected
// sampled p13 tangent bank.
//
// Primary line:
//   U2 = [-1071,-1054,1116,1134*T^2].
// Reciprocal cross-check (integrally scaled presentation):
//   U1 = [-1071,1116,1134,-1054*T^2].
//////////////////////////////////////////////////////////////////////
SetColumns(0);
Z:=Integers(); p:=13;
bank_path:="results/target_22224_direct_contact_deep13_padic_tangent_p13.tsv";
log_file:="results/target_22224_U2_global13_P1_scan.log";
output_file:="results/target_22224_U2_global13_P1_scan.tsv";
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

families:=[*
    [*"U2",[Z!-1071,Z!-1054,Z!1116],Z!1134*],
    [*"U1",[Z!-1071,Z!1116,Z!1134],Z!-1054*]
*];
out:=Open(output_file,"w");
fprintf out,"family\texponent\tchart\tparameter\tkey\n";
print "U2_GLOBAL13_P1_SCAN_START","bank",bank_path;

for fam in families do
    name:=fam[1]; abc:=fam[2]; d0:=fam[3];
    print "U2_GLOBAL13_FAMILY",name,"fixed",abc,"moving_coefficient",d0;
    for e in [1..5] do
        M:=p^e; bank:=ReadBank(bank_path,p,M);
        hits:=0; hitkeys:={}; allkeys:={}; finite:=0; infinity:=0;

        // D is a unit: [N:D]=[t:1].
        for t in [0..M-1] do
            key:=DeepKey(abc cat [d0*t^2],p,M);
            finite+:=1; Include(~allkeys,key);
            if key in bank then
                hits+:=1; Include(~hitkeys,key);
                fprintf out,"%o\t%o\tfinite\t%o\t%o\n",name,e,t,key;
            end if;
        end for;

        // N is a unit, D is not: [N:D]=[1:13*u].
        for u in [0..p^(e-1)-1] do
            s:=p*u;
            key:=DeepKey([abc[i]*s^2:i in [1..3]] cat [d0],p,M);
            infinity+:=1; Include(~allkeys,key);
            if key in bank then
                hits+:=1; Include(~hitkeys,key);
                fprintf out,"%o\t%o\tinfinity\t%o\t%o\n",name,e,u,key;
            end if;
        end for;

        print "U2_GLOBAL13_LEVEL","family",name,"exponent",e,
              "modulus",M,"bank_keys",#bank,
              "finite_classes",finite,"infinity_classes",infinity,
              "projective_classes",finite+infinity,
              "line_keys",#allkeys,"hit_classes",hits,"hit_keys",#hitkeys;
        if e eq 1 then
            print "U2_GLOBAL13_MOD13_BANK_KEYS",Sort(Setseq(bank));
            print "U2_GLOBAL13_MOD13_LINE_KEYS",name,Sort(Setseq(allkeys));
        end if;
    end for;
end for;
delete out;
print "U2_GLOBAL13_P1_SCAN_DONE","output",output_file;
UnsetLogFile(); quit;
