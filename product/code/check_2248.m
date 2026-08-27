// check_2248.m — mod-p profile meets for the ten HLP [2,2,4,8] anchor models
// y^2 = x(x+a^2)(x+b^2)(x+c^2)(x+d^2) in data/ten2248models_abcd.txt.
// Goal: upper bound J(Q)_tors via J(F_p) meets; the quintic is fully split so
// (Z/2)^4 <= J(Q)[2] is automatic.  OUTCOME (logs/check2248{,b}.log): the meets
// do NOT pin the exact group -- [2,2,4,16] stays compatible for all ten rows
// (the gcd #J(F_p) = 256 is automatic here: both elliptic factors carry
// rational [2,8]-torsion, so 16*16 divides #J(F_p) at every good p), and the
// HLP containment [2,2,4,8] <= J(Q)_tors does not exclude that supergroup.
// The table's exact-[2,2,4,8] claim therefore rests on the small lane_hlp37
// witness (logs/verify_witnesses2.log), NOT on these anchors; the anchors are
// recorded as containment, with |J(Q)_tors| in {128, 256}.
// Usage: magma -b check_2248.m > ../logs/check2248.log
SetColumns(0);
SetMemoryLimit(4*10^9);
load "split_lab.m";  // run from product/code/

lines := Split(Read("../../data/ten2248models_abcd.txt"), "\n");
row := 0;
for line in lines do
    if #line lt 5 then continue; end if;
    row +:= 1;
    vals := eval line;   // [a,b,c,d]
    a := vals[1]; b := vals[2]; c := vals[3]; d := vals[4];
    g := x*(x+a^2)*(x+b^2)*(x+c^2)*(x+d^2) where x := RQx.1;
    invsList, ps := ProfileCurve(g, 30);
    printf "ROW %o primes %o\n", row, ps;
    for i in [1..#invsList] do printf "  p=%o invs=%o\n", ps[i], invsList[i]; end for;
    compat := CompatibleGroups(invsList);
    // impose the known lower bound (Z/2)^4 (five rational Weierstrass roots)
    big := [ v : v in compat | #v eq 4 and IsDivisibleBy(v[1],2) ];
    printf "ROW %o compat-with-2rank4: %o\n", row, big;
    ords := [ IsEmpty(iv) select 1 else &*iv : iv in invsList ];
    printf "ROW %o gcd #J(F_p) = %o\n", row, GCD(ords);
end for;
printf "SEARCH_DONE check2248\n";
quit;
