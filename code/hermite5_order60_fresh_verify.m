// Exact fresh-prime verification of the strong refined Hermite-sieve rows.

Q := Rationals();
P<x> := PolynomialRing(Q);

function HermiteData(s,t,m)
    e := t^2-s^2-1;
    q := x^2+e*x+s^2;
    A := x^5+m*x^4
        + (-2*m-s^5/2+5*s^3*t^2/2-5*s^3/2-5*s^2*t^3/2
           +t^5/2+5*t^3/2-3)*x^3
        + (m+2*s^5-5*s^3*t^2+5*s^3+5*s^2*t^3/2
           +t^5/2-5*t^3/2+2)*x^2
        +(5*s^3*e/2)*x+s^5;
    F,r := Quotrem(A^2-q^5,x^2*(x-1)^2);
    assert r eq 0;
    return F,q,A;
end function;

function CheckRows(label,rows,primes)
    kills := [0 : p in primes];
    bad := [0 : p in primes];
    survivors := [];
    for row in rows do
        s,t,m := Explode(row);
        F,q,A := HermiteData(s,t,m);
        assert A^2-x^2*(x-1)^2*F eq q^5;
        assert Degree(F) eq 5 and IsSquarefree(F);
        assert IsSquarefree(q) and Degree(GCD(q,F)) eq 0;
        pass := true;
        for j in [1..#primes] do
            p := primes[j];
            if Denominator(s)*Denominator(t)*Denominator(m) mod p eq 0 then
                bad[j] +:= 1;
                continue;
            end if;
            fp := PolynomialRing(GF(p))!F;
            if Degree(fp) ne 5 or not IsSquarefree(fp) then
                bad[j] +:= 1;
                continue;
            end if;
            n := #Jacobian(HyperellipticCurve(fp));
            r2 := #Factorization(fp)-1;
            if n mod 60 ne 0 or Valuation(n,2) le r2 then
                kills[j] +:= 1;
                pass := false;
                break;
            end if;
        end for;
        if pass then Append(~survivors,row); end if;
    end for;
    print "VERIFY",label,"rows",#rows;
    print "PRIMES",primes;
    print "KILL_COUNTS",kills;
    print "BAD_BEFORE_TEST_COUNTS",bad;
    print "FRESH_SURVIVORS",survivors;
    return kills,bad,survivors;
end function;

integer_rows := [
    [Q|-972,711,-828],[-942,-456,851],[-914,848,419],[-848,914,-424],
    [-796,200,375],[-795,-715,-143],[-753,-59,-830],[-736,252,924],
    [-711,972,823],[-710,-762,807],[-681,637,-588],[-654,-480,-661],
    [-645,246,416],[-637,681,583],[-603,378,-745],[-585,-479,-265],
    [-519,-985,587],[-487,-438,-825],[-473,-504,-227],[-473,-384,923],
    [-470,-885,-669],[-435,-692,68],[-378,603,740],[-332,-924,-680],
    [-330,190,-968],[-268,-757,265],[-252,736,-929],[-246,645,-421],
    [-223,128,-39],[-202,-948,215],[-200,796,-380],[-190,330,963],
    [-128,223,34],[59,753,825],[125,-206,-956],[178,-694,179],
    [206,-125,951],[208,-927,69],[226,-255,35],[243,-441,-297],
    [251,-479,15],[255,-226,-40],[384,473,-928],[405,-996,876],
    [438,487,820],[441,-243,292],[456,942,-856],[479,-251,-20],
    [479,585,260],[480,654,656],[504,473,222],[545,-649,122],
    [567,-922,142],[649,-545,-127],[692,435,-73],[694,-178,-184],
    [715,795,138],[757,268,-270],[762,710,-812],[804,-930,-961],
    [811,-978,-129],[832,-906,-949],[885,470,664],[906,-832,944],
    [922,-567,-147],[924,332,675],[927,-208,-74],[930,-804,956],
    [948,202,-220],[978,-811,124],[985,519,-592],[996,-405,-881]
];

rational_rows := [
    [Q|10,-4,-27/13], [Q|21,-12/7,-4/3], [Q|-2/7,-16/3,29/15],
    [Q|12/7,-21,-11/3], [Q|-27/16,13/3,-7/27],
    [Q|15/22,-4/27,-12/13], [Q|-9/23,4/23,-15/28],
    [Q|4/25,16/27,-14/15], [Q|4/27,11/25,11/10]
];

fresh_primes := [43,47,53,59,61,67,71,73,79,83,89,97,101];
ki,bi,si := CheckRows("INTEGER_H1000",integer_rows,fresh_primes);
kr,br,sr := CheckRows("RATIONAL_H30",rational_rows,fresh_primes);
assert ki[1..3] eq [58,12,2] and &+ki eq 72 and #si eq 0;
assert kr[1..2] eq [6,3] and &+kr eq 9 and #sr eq 0;
print "HERMITE5_ORDER60_FRESH_VERIFY_PASS";
