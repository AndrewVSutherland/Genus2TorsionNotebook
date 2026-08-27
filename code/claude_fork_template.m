// claude_fork_template.m -- VERIFIED working idiom for Magma's Fork/WaitForAllChildren
// (Magma 2.29-8 on claude-box and on the aws-spot boxes).
//
// Fork() is a RAW POSIX fork: it returns the child's PID in the parent and 0 in the child.
// There is no return-value plumbing -- children must write their results to per-child FILES
// and the parent must read them back after WaitForAllChildren().  Three rules:
//   1. the child MUST end with quit;  otherwise it falls through into the parent's loop
//      and you get 2^n processes;
//   2. do all expensive shared setup in the PARENT before forking -- children inherit it
//      copy-on-write, so N children of a parent holding X GB do not each copy X GB unless
//      they write to it;
//   3. Flush (or delete) the child's file handle before quit;, or you lose buffered output.
// Memory: budget SetMemoryLimit per child, not per box.  On a 370 GB / 192-core spot box,
// ~2*10^9 per child when fanning out 192 ways; fewer/fatter children (24 x 12 GB, 8 x 40 GB)
// for memory-hungry symbolic work.  On claude-box (128 GB / 32 vCPU) use ~16-24 children
// at 4*10^9.
//
// Demo below: shard 74 primes over 8 children, each counting roots of a shared polynomial.
// Runs in well under a second and prints FORK_DEMO_OK.

SetColumns(0);
// Magma 2.29-8: Fork() is a POSIX fork -- returns child PID in the parent, 0 in the child.
// Expensive shared setup happens ONCE in the parent and is inherited copy-on-write.
P<x> := PolynomialRing(Rationals());
shared := &*[x - i : i in [1..40]];            // pretend this was expensive
primes := PrimesInInterval(11, 400);
NCH := 8;
dir := "/tmp/claude-1000/-home-claude-torsion-jac/1c22f408-33bb-486e-870f-c5fca84bb522/scratchpad/forkout";
System("mkdir -p " cat dir);
for c in [0..NCH-1] do
  pid := Fork();
  if pid eq 0 then                              // ---- child ----
    F := Open(dir cat "/part" cat IntegerToString(c) cat ".txt", "w");
    for i in [c+1..#primes by NCH] do
      p := primes[i];
      Puts(F, IntegerToString(p) cat " " cat IntegerToString(#Roots(PolynomialRing(GF(p))!shared)));
    end for;
    Flush(F); delete F;
    quit;                                       // child MUST exit
  end if;
end for;
WaitForAllChildren();                           // ---- parent ----
tot := 0; n := 0;
for c in [0..NCH-1] do
  for l in Split(Read(dir cat "/part" cat IntegerToString(c) cat ".txt"), "\n") do
    if #l gt 0 then s := Split(l, " "); tot +:= StringToInteger(s[2]); n +:= 1; end if;
  end for;
end for;
printf "FORK_DEMO_OK children=%o primes_done=%o total_roots=%o\n", NCH, n, tot;
quit;
