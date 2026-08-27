//////////////////////////////////////////////////////////////////////
// Complete-chart direct M(2,2,4,8) cover sweep on tor2244.txt.
//
// Unlike the standard-chart driver, this tests all 30 ordered choices of
// branch points sent to 0 and infinity.  For every resulting square model
// it applies all finite-root permutations and all cover sign choices from
// m2248_sieve.m.  This is the normalization used by the known HPL positive
// controls and closes the chart gap of the raw tuple test.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetSeed(1);

if assigned input_file then
    bank_file := input_file;
    delete input_file;
elif not assigned bank_file then
    bank_file := "data/tor2244_bank.txt";
end if;
if not assigned output_file then
    output_file := "results/target_2248_tor2244_all_charts_candidates.txt";
end if;
if not assigned log_file then
    log_file := "results/target_2248_tor2244_all_charts.log";
end if;
if not assigned first_row then first_row:=1;
elif Type(first_row) eq MonStgElt then first_row:=StringToInteger(first_row); end if;
if not assigned max_rows then max_rows:=0;
elif Type(max_rows) eq MonStgElt then max_rows:=StringToInteger(max_rows); end if;
if not assigned progress_interval then progress_interval:=250;
elif Type(progress_interval) eq MonStgElt then
    progress_interval:=StringToInteger(progress_interval);
end if;

load "code/m2248_sieve.m";

Q:=Rationals();

function FiniteBranchRoots(vals)
    return [Q!0] cat [-Q!v^2:v in vals];
end function;

function MobiusValue(idx,zero_idx,inf_idx,roots)
    if zero_idx eq 6 then
        w:=roots[inf_idx];
        return Q!1/(roots[idx]-w);
    elif inf_idx eq 6 then
        z:=roots[zero_idx];
        return roots[idx]-z;
    else
        z:=roots[zero_idx]; w:=roots[inf_idx];
        if idx eq 6 then return Q!1; end if;
        return (roots[idx]-z)/(roots[idx]-w);
    end if;
end function;

function NormalizedSquareTuple(vals,zero_idx,inf_idx)
    roots:=FiniteBranchRoots(vals);
    remaining:=[i:i in [1..6]|i ne zero_idx and i ne inf_idx];
    lambdas:=[MobiusValue(i,zero_idx,inf_idx,roots):i in remaining];
    if &or[z eq 0:z in lambdas] then return false,[]; end if;
    base:=lambdas[1]; tuple:=[Q!1];
    for j in [2..4] do
        ok,rt:=IsSquare(lambdas[j]/base);
        if not ok then return false,[]; end if;
        Append(~tuple,rt);
    end for;
    return true,tuple;
end function;

tuples:=ReadM2248TupleFile(bank_file);
first:=Max(1,first_row);
last:=max_rows gt 0 select Min(#tuples,first+max_rows-1) else #tuples;

SetLogFile(log_file : Overwrite:=true);
out:=Open(output_file,"w");
print "TARGET_2248_TOR2244_ALL_CHARTS_START";
print "CONFIG","input_file",bank_file,"rows_available",#tuples,
      "first",first,"last",last,"ordered_zero_infinity_charts",30;

rows_seen:=0; chart_trials:=0; square_charts:=0; full_charts:=0;
candidate_rows:=0; full_witnesses:=0;
for ri in [first..last] do
    rows_seen+:=1; vals:=tuples[ri]; row_has_full:=false;
    for zi in [1..6] do for ii in [1..6] do
        if zi eq ii then continue; end if;
        chart_trials+:=1;
        ok,nt:=NormalizedSquareTuple(vals,zi,ii);
        if not ok then continue; end if;
        square_charts+:=1;
        fw:=M2248WitnessesForTupleAllPermutations(nt,true);
        if #fw eq 0 then continue; end if;
        full_charts+:=1; full_witnesses+:=#fw; row_has_full:=true;
        print "FULL_CHART","row",ri,"source_tuple",vals,
              "zero_idx",zi,"infinity_idx",ii,"normalized_tuple",nt,
              "witnesses",#fw;
        for W in fw do
            fprintf out,
              "row=%o source=%o zero=%o infinity=%o normalized=%o ordered=%o perm=%o eps=[%o,%o,%o] rho=%o sigma=%o tau=%o F=[%o,%o,%o,%o]\n",
              ri,vals,zi,ii,nt,W`tuple,W`permutation,W`eps_rho,
              W`eps_sigma,W`eps_tau,W`rho,W`sigma,W`tau,
              W`F1,W`F2,W`F3,W`F4;
        end for;
    end for; end for;
    if row_has_full then candidate_rows+:=1; end if;
    if progress_interval gt 0 and rows_seen mod progress_interval eq 0 then
        print "PROGRESS","rows_seen",rows_seen,"current_row",ri,
              "chart_trials",chart_trials,"square_charts",square_charts,
              "candidate_rows",candidate_rows,"full_charts",full_charts,
              "full_witnesses",full_witnesses;
    end if;
end for;

delete out;
print "SUMMARY","rows_seen",rows_seen,"chart_trials",chart_trials,
      "square_charts",square_charts,"candidate_rows",candidate_rows,
      "full_charts",full_charts,"full_witnesses",full_witnesses;
print "TARGET_2248_TOR2244_ALL_CHARTS_DONE";
UnsetLogFile();
quit;
