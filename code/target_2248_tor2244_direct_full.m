//////////////////////////////////////////////////////////////////////
// Direct all-permutation full-cover sweep for
//     M(2,2,4,8) -> A(2,2,4,4)
// on the existing tor2244 tuple bank.
//
// This is a logged driver around code/m2248_sieve.m.  It complements the
// family-level Richelot source search by testing the target cover directly,
// without assuming an intermediate one-split [2,4,8] source.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned input_file then
    input_file := "data/tor2244_bank.txt";
end if;
if not assigned output_file then
    output_file := "results/target_2248_tor2244_direct_full.txt";
end if;
if not assigned log_file then
    log_file := "results/target_2248_tor2244_direct_full.log";
end if;
if not assigned max_rows then
    max_rows := 0;
elif Type(max_rows) eq MonStgElt then
    max_rows := StringToInteger(max_rows);
end if;
if not assigned progress_interval then
    progress_interval := 1000;
elif Type(progress_interval) eq MonStgElt then
    progress_interval := StringToInteger(progress_interval);
end if;

full_check := true;
all_permutations := true;

SetLogFile(log_file : Overwrite := true);
print "TARGET_2248_TOR2244_DIRECT_FULL_START";
print "CONFIG", "input_file", input_file, "output_file", output_file,
      "max_rows", max_rows, "all_permutations", all_permutations,
      "full_check", full_check;
load "code/m2248_sieve.m";
print "TARGET_2248_TOR2244_DIRECT_FULL_DONE";
UnsetLogFile();
quit;
