#!/usr/bin/env python3
"""
Compare projected PC coordinates from a persistent-model projection against a
stock smartpca baseline .evec.

Usage:
  compare_evec.py BASELINE.evec PROJECTED.evec [--tol 1e-4] [--npc 20]

Matches individuals by ID (only those present in the projected file, i.e. the
test individuals), aligns the first --npc principal components, and reports:
  - per-PC maximum absolute difference
  - overall maximum absolute difference
  - overall mean absolute difference

Exit status is non-zero if the maximum absolute difference exceeds --tol.

No numpy dependency (pure Python) so it runs anywhere.
"""
import argparse
import sys


def read_evec(path):
    """Return dict id -> list[float] of PC coordinates (drops eigval header)."""
    out = {}
    with open(path) as f:
        for line in f:
            s = line.split()
            if not s:
                continue
            if s[0].startswith("#"):     # #eigvals header
                continue
            iid = s[0]
            # last column is the population/qtmode label; coords are between
            vals = []
            for tok in s[1:]:
                try:
                    vals.append(float(tok))
                except ValueError:
                    break                # hit the trailing label
            out[iid] = vals
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("baseline")
    ap.add_argument("projected")
    ap.add_argument("--tol", type=float, default=1e-4)
    ap.add_argument("--npc", type=int, default=20)
    args = ap.parse_args()

    base = read_evec(args.baseline)
    proj = read_evec(args.projected)

    ids = [i for i in proj if i in base]
    if not ids:
        print("ERROR: no shared individual IDs between files")
        return 2

    npc = args.npc
    for i in ids:
        npc = min(npc, len(base[i]), len(proj[i]))
    if npc <= 0:
        print("ERROR: no comparable PCs")
        return 2

    per_pc_max = [0.0] * npc
    total_abs = 0.0
    total_n = 0
    global_max = 0.0
    worst = None
    for iid in ids:
        for j in range(npc):
            d = abs(base[iid][j] - proj[iid][j])
            per_pc_max[j] = max(per_pc_max[j], d)
            total_abs += d
            total_n += 1
            if d > global_max:
                global_max = d
                worst = (iid, j + 1, base[iid][j], proj[iid][j])

    mean_abs = total_abs / total_n if total_n else 0.0

    print("compared %d individuals over PCs 1-%d" % (len(ids), npc))
    print("per-PC maximum absolute difference:")
    for j in range(npc):
        print("  PC%-2d : %.3e" % (j + 1, per_pc_max[j]))
    print("-" * 40)
    print("maximum absolute difference : %.3e" % global_max)
    print("mean absolute difference    : %.3e" % mean_abs)
    if worst:
        print("worst: %s PC%d baseline=%.6f projected=%.6f"
              % (worst[0], worst[1], worst[2], worst[3]))
    print("tolerance                   : %.3e" % args.tol)

    if global_max > args.tol:
        print("RESULT: FAIL (max diff exceeds tolerance)")
        return 1
    print("RESULT: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
