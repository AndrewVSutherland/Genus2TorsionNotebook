#!/usr/bin/env python3
"""sat_account.py - stitch the full-saturation certificate results
(sigsat_si<k>.log) across resumed segments.  Per fiber (deduped by u),
class priority: SATCERT/SATCERTM (certified idx=1) beats SATOPENB/SATOPENM
(uncertified) beats SATSKIP - a fiber retried in a later segment keeps its
best outcome.  SATVIOLATION is never absorbed: always reported."""
import re, sys

ULISTS = {}
cur = 2
for line in open("../data/sat_usets.m"):
    m = re.match(r"\[ Rationals\(\) \| (.*) \]", line.strip().rstrip(","))
    if m:
        ULISTS[cur] = [x.strip() for x in m.group(1).split(",") if x.strip()]
        cur += 1
USIZES = {k: len(v) for k, v in ULISTS.items()}

PRI = {"cert": 3, "certm": 3, "openb": 2, "openm": 2, "skip": 1, "crash": 0.5}
tot = {k: 0 for k in ("fibers", "cert", "certm", "openb", "openm", "skip", "viol", "crash")}
print(f"{'sig':>4} {'uset':>5} {'seen':>5} {'cert':>5} {'certm':>6} {'openB':>6} {'openM':>6} {'skip':>5} {'viol':>5} {'miss':>5}")
openb_lists = {}
for k in range(2, 7):
    cls = {}
    viol = set()
    # read the kill-proof Rec file first (same marker grammar), then the log
    # (for CRASHSKIP/VIOLATION and any pre-Rec-era segments)
    srcs = []
    for cand in (f"../logs/sigsatR_si{k}.txt", f"../logs/sigsat_si{k}.log"):
        try: srcs.append(open(cand))
        except FileNotFoundError: pass
    import itertools
    if srcs:
        fh = itertools.chain(*srcs)
    try:
        if not srcs: raise FileNotFoundError
    except FileNotFoundError:
        # absent log = entirely unrun surface: keep it in the denominator
        # and carry ALL its fibers into the follow-up list (codex)
        print(f"si{k:<2} {USIZES.get(k,0):>5} {0:>5} {0:>5} {0:>6} {0:>6} {0:>6} {0:>5} {0:>5} {USIZES.get(k,0):>5}  (no log - unrun)")
        tot["fibers"] += USIZES.get(k, 0)
        openb_lists[k] = list(ULISTS.get(k, []))
        continue
    for line in fh:
        for tag, c in (("SATCERTM", "certm"), ("SATCERT2", "certm"), ("SATCERT", "cert"),
                       ("SATOPENB", "openb"), ("SATOPENM", "openm"), ("SATSKIP", "skip")):
            m = re.match(tag + r" u=([^\s]+)", line)
            if m:
                u = m.group(1)
                if PRI[c] > PRI.get(cls.get(u, ""), 0):
                    cls[u] = c
                break
        else:
            m = re.match(r"SATVIOLATION u=([^\s]+)", line)
            if m: viol.add(m.group(1))
            else:
                # resolve CRASHSKIP to its fiber's u; lowest priority so any
                # later successful classification supersedes it (codex)
                m = re.match(r"CRASHSKIP fiber (\d+) ", line)
                if m:
                    idx = int(m.group(1))
                    if 1 <= idx <= len(ULISTS.get(k, [])):
                        u = ULISTS[k][idx-1]
                        if u not in cls:
                            cls[u] = "crash"
    n = {c: sum(1 for v in cls.values() if v == c) for c in ("cert", "certm", "openb", "openm", "skip", "crash")}
    crash = n["crash"]
    seen = len(cls)
    miss = USIZES.get(k, 0) - seen
    print(f"si{k:<2} {USIZES.get(k,0):>5} {seen:>5} {n['cert']:>5} {n['certm']:>6} "
          f"{n['openb']:>6} {n['openm']:>6} {n['skip']:>5} {len(viol):>5} {miss:>5}")
    if viol: print(f"      VIOLATIONS: {sorted(viol)}")
    # follow-up list = openb-classified PLUS every USET fiber not yet seen
    # (a partial/interrupted run must not drop its unprocessed tail - codex)
    unseen = [u for u in ULISTS.get(k, []) if u not in cls]
    openb_lists[k] = sorted(set(u for u, c in cls.items() if c == "openb") | set(unseen))
    for key in ("cert", "certm", "openb", "openm", "skip"):
        tot[key] += n[key]
    tot["fibers"] += USIZES.get(k, 0); tot["viol"] += len(viol); tot["crash"] += crash
print("TOTALS:", tot)
certified = tot["cert"] + tot["certm"]
print(f"certified idx=1: {certified}/{tot['fibers']} "
      f"({100.0*certified/max(tot['fibers'],1):.1f}%); "
      f"uncertified: openB {tot['openb']} openM {tot['openm']} skip {tot['skip']}")
if len(sys.argv) > 1 and sys.argv[1] == "--write-openb":
    out = ["[*"]
    for k in range(2, 7):
        us = openb_lists.get(k, [])
        out.append("[ Rationals() | " + ", ".join(us) + " ]" + ("," if k < 6 else ""))
    out.append("*]")
    open("../data/sat_usets_openb.m", "w").write("\n".join(out) + "\n")
    print("wrote ../data/sat_usets_openb.m")
