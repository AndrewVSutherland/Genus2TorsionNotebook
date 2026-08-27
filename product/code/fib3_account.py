#!/usr/bin/env python3
"""fib3_account.py - stitch per-fiber dispositions of the saturated
completion pass (sigfib3_si<k>.log) across crash-resumed segments.
A fiber's class is decided by its marker lines, deduped by u:
  GENS2 present & not GENSHORT/RBOPEN/SATWARN -> complete
  GENSHORT -> genshort; RBOPEN -> rbopen; SATWARN -> satwarn
  NOBASE -> nobase; MWSKIP -> mwskip; CRASHSKIP (driver) -> crashskip
Also cross-checks coverage against the USET sizes in fib2_usets.m."""
import re, sys

USIZES = {}
cur = 2
for line in open("../data/fib2_usets.m"):
    m = re.match(r"\[ Rationals\(\) \| (.*) \]", line.strip().rstrip(","))
    if m:
        USIZES[cur] = len([x for x in m.group(1).split(",") if x.strip()])
        cur += 1

tot = {k: 0 for k in ("fibers", "complete", "genshort", "rbopen", "satwarn",
                      "nobase", "mwskip", "crashskip", "surv", "deck", "instance")}
print(f"{'sig':>4} {'uset':>5} {'seen':>5} {'compl':>6} {'short':>6} {'rbopen':>6} "
      f"{'satw':>5} {'nobase':>6} {'mwskip':>6} {'crash':>5} {'miss':>5}")
for k in range(2, 7):
    gens2, short, rbopen, satwarn, nobase, mwskip = set(), set(), set(), set(), set(), set()
    crashskip = 0
    surv, deck, inst = set(), set(), set()
    for line in open(f"../logs/sigfib3_si{k}.log"):
        m = re.match(r"GENS2 u=([^\s]+) ", line)
        if m: gens2.add(m.group(1)); continue
        m = re.match(r"GENSHORT u=([^\s]+) ", line)
        if m: short.add(m.group(1)); continue
        m = re.match(r"RBOPEN u=([^\s]+) ", line)
        if m: rbopen.add(m.group(1)); continue
        m = re.match(r"SATWARN u=([^\s]+) ", line)
        if m: satwarn.add(m.group(1)); continue
        m = re.match(r"NOBASE u=([^\s]+)", line)
        if m: nobase.add(m.group(1)); continue
        m = re.match(r"MWSKIP u=([^\s]+) ", line)
        if m: mwskip.add(m.group(1)); continue
        if line.startswith("CRASHSKIP"): crashskip += 1; continue
        # survivors/deck/instances dedupe by (t,u) - resumed segments can
        # overlap a few fibers and repeat marker lines (codex round 6)
        m = re.match(r"FIBSURV Sig=\d+ t=([^\s]+) u=([^\s]+)", line)
        if m: surv.add((m.group(1), m.group(2))); continue
        m = re.match(r"SKIPISO t=([^\s]+) u=([^\s]+)", line)
        if m: deck.add((m.group(1), m.group(2))); continue
        m = re.match(r"INSTANCE \S+ t=([^\s]+) u=([^\s]+)", line)
        if m: inst.add((m.group(1), m.group(2))); continue
    complete = gens2 - short - rbopen - satwarn
    seen = len(gens2 | nobase | mwskip) + crashskip
    miss = USIZES.get(k, 0) - seen
    print(f"si{k:<2} {USIZES.get(k,0):>5} {seen:>5} {len(complete):>6} {len(short):>6} "
          f"{len(rbopen):>6} {len(satwarn):>5} {len(nobase):>6} {len(mwskip):>6} "
          f"{crashskip:>5} {miss:>5}")
    if rbopen: print(f"      RBOPEN u: {sorted(rbopen)}")
    if satwarn: print(f"      SATWARN u: {sorted(satwarn)}")
    if nobase: print(f"      NOBASE u: {sorted(nobase)}")
    if mwskip: print(f"      MWSKIP u: {sorted(mwskip)}")
    if inst: print(f"      !!! INSTANCE pairs: {sorted(inst)}")
    tot["fibers"] += USIZES.get(k, 0); tot["complete"] += len(complete)
    tot["genshort"] += len(short); tot["rbopen"] += len(rbopen)
    tot["satwarn"] += len(satwarn); tot["nobase"] += len(nobase)
    tot["mwskip"] += len(mwskip); tot["crashskip"] += crashskip
    tot["surv"] += len(surv); tot["deck"] += len(deck); tot["instance"] += len(inst)
print("TOTALS:", {k: v for k, v in tot.items()})
acc = tot["complete"] + tot["genshort"] + tot["rbopen"] + tot["satwarn"] + \
      tot["nobase"] + tot["mwskip"] + tot["crashskip"]
print(f"accounted {acc} of {tot['fibers']} fibers; survivors {tot['surv']} "
      f"(deck-skipped {tot['deck']}), genuine instances {tot['instance']}")
