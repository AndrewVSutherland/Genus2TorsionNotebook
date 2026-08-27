//////////////////////////////////////////////////////////////////////
// Mod-125 continuation of every smooth compensated omega-pole residue.
//
// The rank-5 Jacobian has a 2-dimensional correction kernel, so each
// residue has exactly 25 lifts at every digit.  We enumerate those affine
// kernels (rather than 5^7 ambient corrections) through modulus 125 and
// test whether U+2*nu finally leaves the degenerate contact component.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Z:=Integers(); k:=GF(5);
R<S,O,m,D,L,U,nu>:=PolynomialRing(Z,7,"grevlex");
a:=S*O^2-3; b:=25*S-3; M:=L^2;
c5:=50*S; c4:=b^2+2*a-15; c3:=2*a*b+22;
c2:=a^2+2*b-15; c1:=2*(a+3);
G:=m^2-c4-5*D;
H:=(8*S*(c3-20*S*O)-D^2)^2
   -256*S^2*m^2*(c2+O*D);
B3:=c5*M+3*U;
Delta3:=4*c4*M+12*(U^2+nu^2)-B3^2;
F3:=B3*Delta3+16*nu^3-8*c3*M-8*U^3-48*U*nu^2;
F2:=Delta3^2+64*B3*nu^3-64*c2*M
    -192*(U^2*nu^2+nu^4);
F1:=Delta3*nu^3-4*c1*M-12*U*nu^4;
eqs:=[G,H,F1,F2,F3];
derivs:=[[Derivative(f,j):j in [1..7]]:f in eqs];

function LiftVectors(A,rhs)
    B:=Transpose(A);
    part:=Solution(B,rhs);
    K:=Kernel(B); assert Dimension(K) eq 2;
    bas:=Basis(K);
    return [part+a0*bas[1]+b0*bas[2]:a0 in k,b0 in k];
end function;

residues:=0; lifts25:=0; lifts125:=0; escape125:=0;
by_SO:=AssociativeArray(); samples25:=[]; samples125:=[];
for s0 in [0..4] do for o0 in [0..4] do for m0 in [0..4] do
for d0 in [0..4] do for l0 in [0..4] do for u0 in [0..4] do
for n0 in [0..4] do
    x0:=[s0,o0,m0,d0,l0,u0,n0];
    if s0*o0*m0*l0*n0*(u0-2*n0) eq 0 then continue; end if;
    f0:=[Z!Evaluate(f,x0):f in eqs];
    if not &and[z mod 5 eq 0:z in f0] then continue; end if;
    A:=Matrix(k,5,7,
        &cat[[k!((Z!Evaluate(derivs[i][j],x0)) mod 5):j in [1..7]]
              :i in [1..5]]);
    if Rank(A) ne 5 then continue; end if;
    residues+:=1;
    rhs1:=Vector(k,[k!((-(f0[i] div 5)) mod 5):i in [1..5]]);
    for v1 in LiftVectors(A,rhs1) do
        d1:=[Z!z:z in Eltseq(v1)];
        x25:=[x0[i]+5*d1[i]:i in [1..7]];
        assert &and[(Z!Evaluate(f,x25)) mod 25 eq 0:f in eqs];
        assert (x25[6]+2*x25[7]) mod 25 eq 0;
        lifts25+:=1;
        if #samples25 lt 8 then Append(~samples25,<x0,d1,x25>); end if;
        f25:=[Z!Evaluate(f,x25):f in eqs];
        rhs2:=Vector(k,[k!((-(f25[i] div 25)) mod 5):i in [1..5]]);
        for v2 in LiftVectors(A,rhs2) do
            d2:=[Z!z:z in Eltseq(v2)];
            x125:=[x25[i]+25*d2[i]:i in [1..7]];
            assert &and[(Z!Evaluate(f,x125)) mod 125 eq 0:f in eqs];
            lifts125+:=1;
            um:=(x125[6]+2*x125[7]) mod 125;
            if um eq 0 then continue; end if;
            escape125+:=1;
            key:=Sprint(<s0,o0>);
            if IsDefined(by_SO,key) then by_SO[key]+:=1; else by_SO[key]:=1; end if;
            if #samples125 lt 32 then
                Append(~samples125,<x0,d1,d2,x125,um>);
            end if;
        end for;
    end for;
end for; end for; end for; end for; end for; end for; end for;

print "CONTACT6_M612_T0_LOCAL5_INTERNAL_MOD125";
print "RANK5_RESIDUES",residues,"LIFTS_MOD25",lifts25,
      "LIFTS_MOD125",lifts125,"ESCAPE_U_PLUS_2NU_MOD125",escape125;
print "ESCAPE_COUNTS_BY_S_O",Sort([<x,by_SO[x]>:x in Keys(by_SO)]);
print "MOD25_SAMPLES_x0_d1_x25",samples25;
print "ESCAPING_MOD125_SAMPLES_x0_d1_d2_x125_Um",samples125;
print "CONTACT6_M612_T0_LOCAL5_INTERNAL_MOD125_DONE";
quit;
