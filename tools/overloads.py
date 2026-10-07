"""Add `b` -- what each violated contingency overloads -- to an existing corpus.

prep.py computes this on a full rebuild; this applies the same `over_of` to the
runs already in docs/data, which is cheaper and leaves their verified `j` alone.
Idempotent: a run that already has `b` on every violated frame is skipped.

    python3 tools/overloads.py docs/data [run_id ...]
"""
import os, sys, json, glob, time
from concurrent.futures import ProcessPoolExecutor

OUT = sys.argv[1] if len(sys.argv) > 1 else 'docs/data'
ONLY = set(sys.argv[2:])                 # read before argv is rewritten, not after
sys.argv = [sys.argv[0], OUT]            # prep.py reads argv[1] as its output dir
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import prep


def main():
    only = ONLY
    files = sorted(glob.glob(os.path.join(OUT, 'runs', '*.json')))
    recs, todo = [], {}
    for fn in files:
        r = json.load(open(fn))
        if only and r['id'] not in only:
            continue
        need = [f for p in r['phases'] for f in p['frames'] if f['v'] and 'b' not in f]
        if not need:
            continue
        recs.append((fn, r))
        for f in need:
            todo.setdefault(r['system'], set()).add((tuple(f['o']), tuple(f['v'])))
    if not recs:
        print('nothing to do')
        return

    ovl, t0 = {}, time.time()
    for s, pairs in todo.items():
        pairs = sorted(pairs)
        print('  %s: %d distinct (topology, v_ctg) pairs' % (s, len(pairs)), flush=True)
        with ProcessPoolExecutor() as ex:
            vals = list(ex.map(prep.over_of, [(s, t, v) for t, v in pairs], chunksize=100))
        ovl[s] = dict(zip(pairs, vals))
    print('  %d solves in %.0f s' % (sum(len(v) for p in ovl.values() for v in p), 
                                     time.time() - t0), flush=True)

    empty = nb = 0
    for fn, r in recs:
        ov = ovl[r['system']]
        for p in r['phases']:
            for f in p['frames']:
                if f['v'] and 'b' not in f:
                    f['b'] = ov[(tuple(f['o']), tuple(f['v']))]
                    for x in f['b']:
                        nb += len(x)
                        empty += not x
        json.dump(r, open(fn, 'w'), separators=(',', ':'))
    print('runs %d  overloaded branches %d  contingencies with no flows %d'
          % (len(recs), nb, empty))


if __name__ == '__main__':
    main()
