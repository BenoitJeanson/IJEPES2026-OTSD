"""Compile tmp/IJEPES/ablation into the static corpus the replay viewer reads.

Frames carry no objective in the callback logs, so every distinct topology is
re-evaluated with coordedit's validated `pf.risk` (eval_risk) -- the campaign's
objective is lost load (ov_lostload_coef=1, ov_n1violation_coef=null).
"""
import os, re, sys, json, glob, csv
import datetime as dt
from concurrent.futures import ProcessPoolExecutor

ROOT = os.environ.get('TNR_ROOT', '/Users/benoitjeanson/vsCode/TUD/tnr')
ABL = os.path.join(ROOT, 'tmp/IJEPES/ablation')
OUT = sys.argv[1]
sys.path.insert(0, os.path.join(ROOT, 'tools/coordedit'))
import pf

SYS = {                       # system -> (coordedit case name, dump, coord csv, TLF)
    'ieee118': ('case118_ieee', 'case118', 'case118.csv', 1.5),
    'ieee57':  ('case57_ieee',  'case57_ieee', 'case57_IEEE.csv', 1.0),
}
E = re.compile(r'\("(\d+)", "(\d+)"\)')
ITER = re.compile(r'^Iteration (\d+)\s*\t(.*?)\topen branches: Set\(\[(.*?)\]\)\s*$', re.M)
DUR = re.compile(r'(\d+)\s+(hour|minute|second|millisecond)s?')
UNIT = {'hour': 3600000, 'minute': 60000, 'second': 1000, 'millisecond': 1}


def ms(s):
    return sum(int(n) * UNIT[u] for n, u in DUR.findall(s))


def load_net(sysname):
    _, dump, coordf, tlf = SYS[sysname]
    d = json.load(open(os.path.join(ROOT, 'tools/coordedit/cases', dump + '.json')))
    edges, idx = [], {}
    for e in d['edges']:
        k = (min(int(e['f']), int(e['t'])), max(int(e['f']), int(e['t'])))
        if k in idx:                       # parallel pair already collapsed by the dump
            continue
        idx[k] = len(edges)
        edges.append({'f': k[0], 't': k[1], 'pmax': e['p_max'] * tlf})
    xy = {}
    with open(os.path.join(ROOT, 'data/exp_raw/coord', coordf)) as fh:
        for r in csv.DictReader(fh):
            xy[int(r['bus_i'])] = [float(r['x']), float(r['y'])]
    buses = [{'i': int(b), 'p': p, 'xy': xy.get(int(b), [0, 0])} for b, p in d['p'].items()]
    buses.sort(key=lambda b: b['i'])
    return {'case': d['case'], 'slack': int(d['slack']), 'tlf': tlf,
            'buses': buses, 'edges': edges}, idx


def parse_cb(path, idx):
    """-> [{t, o:[edge idx], v:[edge idx], c:cuts}]; guards against a concatenated re-solve."""
    txt = open(path).read()
    frames, prev = [], 0
    blocks = txt.split('Iteration ')
    for blk in blocks[1:]:
        m = re.match(r'(\d+)\s*\t(.*?)\topen branches: Set\(\[(.*?)\]\)', blk, re.S)
        if not m:
            continue
        n = int(m.group(1))
        if n <= prev:                      # numbering reset: a second solve in one file
            break
        prev = n
        o = sorted({idx[(min(int(a), int(b)), max(int(a), int(b)))] for a, b in E.findall(m.group(3))
                    if (min(int(a), int(b)), max(int(a), int(b))) in idx})
        tail = blk[m.end():]
        vm = re.search(r'v_ctg Set(?:\{[^}]*\})?\((?:\[(.*?)\])?\)', tail, re.S)
        v = sorted({idx[(min(int(a), int(b)), max(int(a), int(b)))]
                    for a, b in E.findall(vm.group(1) or '')
                    if (min(int(a), int(b)), max(int(a), int(b))) in idx}) if vm else []
        cm = re.search(r'applied cuts: (\d+)', tail)
        frames.append({'t': ms(m.group(2)), 'o': o, 'v': v, 'c': int(cm.group(1)) if cm else 0})
    return frames


NETS = {}


def _net(sysname):
    if sysname not in NETS:
        argv, sys.argv = sys.argv, sys.argv[:1]   # serve.py reads argv[1] as a port
        try:
            import serve
        finally:
            sys.argv = argv
        cname, dump, _, _ = SYS[sysname]
        NETS[sysname] = serve.net_for(cname, dump)
    return NETS[sysname]


REV = {}


def _rev(sysname):
    """edge index -> the (bus, bus) key pf.Net uses."""
    if sysname not in REV:
        net, idx = load_net(sysname)
        REV[sysname] = {v: k for k, v in idx.items()}
    return REV[sysname]


def risk_of(args):
    sysname, topo = args
    rev = _rev(sysname)
    return round(pf.risk(_net(sysname), frozenset(rev[i] for i in topo)), 6)


def main():
    # tools/sbs.jl writes this: the SBS membership, which no text log records.
    # `used[n]` is the set phase n searched in; `after[n]` is what it grew to.
    sbs_path = os.environ.get('SBS_JSON', os.path.join(os.path.dirname(OUT), 'sbs.json'))
    SBS = json.load(open(sbs_path)) if os.path.isfile(sbs_path) else {}
    if not SBS:
        print('  no sbs.json - run tools/sbs.jl first, or the SBS layer will be absent')

    os.makedirs(os.path.join(OUT, 'runs'), exist_ok=True)
    os.makedirs(os.path.join(OUT, 'net'), exist_ok=True)
    nets, idxs = {}, {}
    for s in SYS:
        nets[s], idxs[s] = load_net(s)
        json.dump(nets[s], open(os.path.join(OUT, 'net', s + '.json'), 'w'), separators=(',', ':'))

    runs, todo = [], {s: set() for s in SYS}
    for d in sorted(os.listdir(ABL)):
        rd = os.path.join(ABL, d)
        if d.startswith('_') or not os.path.isfile(os.path.join(rd, 'result.json')):
            continue
        res = json.load(open(os.path.join(rd, 'result.json')))
        man = json.load(open(os.path.join(rd, 'manifest.json')))
        s = res['system']
        if s not in SYS:
            continue
        idx = idxs[s]
        # some dirs hold two campaigns' logs (the run was executed twice); the
        # manifest's start time picks the one result.json actually describes
        # Files are stamped yyyy-mm-dd_HHMMSS from a clock read just before the
        # manifest is written, so the two can differ by a second; and a directory
        # may hold two campaigns. Take the stamp nearest the manifest's start time.
        t0 = dt.datetime.fromisoformat(man['started_at'].split('.')[0])
        stamps = sorted({m for f in os.listdir(rd)
                         for m in re.findall(r'\d{4}-\d{2}-\d{2}_\d{6}', f)})
        stamp = min(stamps, key=lambda x: abs(
            dt.datetime.strptime(x, '%Y-%m-%d_%H%M%S') - t0)) if stamps else ''
        key = lambda lst: sorted({idx[(min(int(a), int(b)), max(int(a), int(b)))]
                                  for a, b in (o.split('-') for o in (lst or []))
                                  if (min(int(a), int(b)), max(int(a), int(b))) in idx})
        phases = []
        for it in res.get('iterations', []):
            g = sorted(glob.glob(os.path.join(rd, '*_i%03d_cb.log' % it['iteration'])))
            g = [x for x in g if os.path.basename(x).startswith(stamp)] or g
            frames = parse_cb(g[-1], idx) if g else []
            sol = key(it['openings'])
            sb = SBS.get(d, {})
            used = key(sb.get('used', {}).get(str(it['iteration'])) or [])
            prev = key(sb.get('used', {}).get(str(it['iteration'] - 1)) or []) \
                if it['iteration'] > 1 else []
            phases.append({'sbs': used,
                           'sbs_new': sorted(set(used) - set(prev)) if it['iteration'] > 1 else [],
                           'sbs_after': key(sb.get('after', {}).get(str(it['iteration'])) or []),
                           'i': it['iteration'], 'obj': it['objective'], 'secure': it['secure'],
                           'lp': it['lp_solves'], 'sec': it['seconds'],
                           'cum': it['cumulative_seconds'], 'capped': it['capped'],
                           'sol': sol, 'frames': frames})
            for f in frames:
                todo[s].add(tuple(f['o']))
            todo[s].add(tuple(sol))
        p = man.get('parameters', {})
        rec = {'id': d, 'system': s, 'config': res['config'], 'block': res.get('block'),
               'H': res['H'], 'd_viol': res['d_viol'], 'd_sol': res['d_sol'],
               'seed': res.get('seed'), 'status': res['status'], 'obj': res['objective'],
               'wall': res['wall_seconds'], 'lp': res['lp_solves'],
               'bi': res['benders_iterations'], 'init_sbs': res.get('init_sbs'),
               'tlf': man.get('tlf'), 'disabled': man.get('disabled_component'),
               'git_sha': man.get('git_sha'), 'cpu': man.get('cpu'),
               'seed_start': key(p.get('heuristic_openings', [])),
               'ctg_master': key(p.get('contingencies_in_master', [])),
               'params': p, 'phases': phases}
        rec['nframes'] = sum(len(ph['frames']) for ph in phases)
        runs.append(rec)

    # --- the column the logs do not have: eval_risk for every distinct topology ---
    obj = {}
    for s, topos in todo.items():
        topos = sorted(topos)
        print('  %s: %d distinct topologies' % (s, len(topos)), flush=True)
        with ProcessPoolExecutor() as ex:
            vals = list(ex.map(risk_of, [(s, t) for t in topos], chunksize=200))
        obj[s] = dict(zip(topos, vals))

    for rec in runs:
        o = obj[rec['system']]
        for ph in rec['phases']:
            ph['solobj'] = o[tuple(ph['sol'])]
            for f in ph['frames']:
                f['j'] = o[tuple(f['o'])]
        json.dump(rec, open(os.path.join(OUT, 'runs', rec['id'] + '.json'), 'w'),
                  separators=(',', ':'))

    index = [{k: r[k] for k in ('id', 'system', 'config', 'block', 'H', 'd_viol', 'd_sol',
                                'seed', 'status', 'obj', 'wall', 'lp', 'bi', 'init_sbs',
                                'disabled', 'nframes')} | {'nphases': len(r['phases'])}
             for r in runs]
    meta = {'runs': index, 'generated': '2026-09-21',
            'git_sha': runs[0]['git_sha'] if runs else None,
            'nframes': sum(r['nframes'] for r in runs),
            'systems': {s: {'tlf': SYS[s][3], 'case': SYS[s][0]} for s in SYS}}
    json.dump(meta, open(os.path.join(OUT, 'index.json'), 'w'), separators=(',', ':'))
    print('runs %d  phases %d  frames %d' % (len(runs), sum(len(r['phases']) for r in runs),
                                             meta['nframes']))


if __name__ == '__main__':
    main()
