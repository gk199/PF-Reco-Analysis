#!/usr/bin/env python3
"""HCAL event displays: reconstructed clusters, member PFRecHits, and depth proxies.

Dependencies: numpy scipy shapely awkward uproot matplotlib. See README_multipion.md.
Selection uses shared detector cells and fractional assigned energies; it is an
illustrative event ranking, not a truth-matched reconstruction efficiency metric.
"""
from __future__ import annotations
import argparse
import json
import math
import tempfile
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
from shapely.geometry import Polygon
import awkward as ak
import uproot
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
from matplotlib.colors import to_rgb, to_hex
from matplotlib.legend_handler import HandlerTuple
from matplotlib.text import Text
from matplotlib.transforms import Bbox
from matplotlib.path import Path as MarkerPath


# Stable colors encode the seed depth of the entire timingPF cluster.
DEPTH_COLORS={1:'#0072B2',2:'#D55E00',3:'#009E73',4:'#6A3D9A',
              5:'#E69F00',6:'#56B4E9',7:'#80512C'}


def depth_color(depth):
    return DEPTH_COLORS.get(depth,'#555555')


def cluster_palette(clusters):
    """Stable depth hues, with distinct cluster shades ordered by stored ID."""
    palette={}
    for depth in sorted({c.seed[2] for c in clusters}):
        ids=sorted(c.cid for c in clusters if c.seed[2]==depth)
        base=np.asarray(to_rgb(depth_color(depth)))
        amounts=[0.] if len(ids)==1 else np.linspace(-.22,.28,len(ids))
        for cid,amount in zip(ids,amounts):
            rgb=base*(1+amount) if amount<0 else base+(1-base)*amount
            palette[cid]=to_hex(rgb)
    return palette


def colored_hit(ax,x,y,colors,seed=False):
    colors=colors or ['0.65']
    size=52 if seed else 26
    if len(colors)==1:
        ax.scatter(x,y,s=size,color=colors[0],alpha=1.,
                   edgecolors='black' if seed else 'white',linewidths=.55,zorder=6 if seed else 4)
    else:
        # Sectors are categorical memberships, not assigned-energy fractions.
        for i,color in enumerate(colors):
            theta=np.linspace(2*np.pi*i/len(colors),2*np.pi*(i+1)/len(colors),25)
            verts=np.vstack([[0.,0.],np.column_stack([np.cos(theta),np.sin(theta)]),[0.,0.]])
            ax.scatter(x,y,s=size,marker=MarkerPath(verts),color=color,
                       alpha=1.,linewidths=0,zorder=6 if seed else 4)
        ax.scatter(x,y,s=size,facecolors='none',edgecolors='black' if seed else 'white',
                   linewidths=.55,zorder=7 if seed else 5)


@dataclass
class Hit:
    key: tuple
    eta: float
    phi: float
    depth: int
    time: float
    energy: float
    fraction: float
    assigned: float

@dataclass
class Cluster:
    cid: int
    eta: float
    phi: float
    time: float
    energy: float
    hits: list
    seed: tuple

@dataclass
class Event:
    key: tuple
    eta0: float
    phi0: float
    gen: list
    clusters: list


def valid(t):
    return math.isfinite(t) and t != -999.0


def allowed_cluster_counts(args, generated_count):
    """Resolve display multiplicities; an explicit list takes precedence."""
    n = args.expected_clusters if args.expected_clusters is not None else generated_count
    if args.timing_clusters:
        return sorted(set(args.timing_clusters))
    if args.timing_count_offsets is not None:
        return sorted({n+d for d in args.timing_count_offsets if n+d > 0})
    if args.cluster_count_mode == 'exact':
        return [n]
    if args.cluster_count_mode == 'range':
        return list(range(2, n+3))
    return list(range(max(2, n-3), n+4))


def highest_energy_hit(c):
    """Histogram seed proxy: raw energy, before checking time validity.

    Exact energy ties are resolved by detector-cell key. No valid-time fallback.
    This proxy is independent of the stored algorithm seed used in event displays.
    """
    hits = [h for h in c.hits if math.isfinite(h.energy)]
    return min(hits, key=lambda h: (-h.energy, h.key)) if hits else None


def timing_spreads(events, source):
    """One entry per qualifying event, with independent cluster/seed selection."""
    entries = []
    for key, ev in sorted(events.items()):
        times = []
        ids = []
        for c in ev.clusters:
            if source == 'cluster':
                t = c.time
            elif source == 'seed':
                h = highest_energy_hit(c)
                t = h.time if h is not None else math.nan
            else:
                raise ValueError(f'Unknown time source {source!r}')
            if valid(t):
                times.append(float(t))
                ids.append(c.cid)
        if len(times) >= 2:
            entries.append(dict(run=key[0], lumi=key[1], event=key[2],
                                cluster_ids=ids, times_ns=times,
                                delta_t_ns=max(times)-min(times)))
    return entries


def histogram_data(std, timing, bin_width):
    result = dict(
        selection='All input events independently per algorithm; at least two valid times; no display cuts',
        valid_time='finite and not -999',
        seed_definition='Highest raw-energy associated rechit per cluster; ties by cell key; no valid-time fallback',
        input_events=dict(standardPF=len(std), seedTimingPF=len(timing)), pages={})
    for source in ('cluster', 'seed'):
        entries = {name: timing_spreads(events, source)
                   for name, events in (('standardPF', std), ('seedTimingPF', timing))}
        largest = max((r['delta_t_ns'] for rows in entries.values() for r in rows), default=0.)
        nbins = max(1, int(math.floor(largest/bin_width))+1)
        if nbins > 100000:
            raise ValueError('Histogram would exceed 100000 bins; increase --hist-bin-width.')
        edges = np.arange(nbins+1, dtype=float)*bin_width
        page = dict(bin_edges_ns=edges.tolist(), algorithms={})
        for name, rows in entries.items():
            counts, _ = np.histogram([r['delta_t_ns'] for r in rows], bins=edges)
            assert int(counts.sum()) == len(rows)
            page['algorithms'][name] = dict(event_count=len(rows), bin_counts=counts.tolist(), events=rows)
        result['pages'][source] = page
    return result


def draw_timing_histogram(pdf, data, source, page_number, total, label=None):
    page = data['pages'][source]
    with plt.rc_context({'font.size': 12, 'axes.linewidth': 1.2}):
        fig, ax = plt.subplots(figsize=(13, 7.3))
        fig.subplots_adjust(left=.12, right=.95, bottom=.19, top=.84)
        for name, color in [('standardPF', '#3399FF'), ('seedTimingPF', '#CC0000')]:
            stats = page['algorithms'][name]
            ax.stairs(stats['bin_counts'], page['bin_edges_ns'], color=color,
                      linewidth=2.2, label=f'{name}: {stats["event_count"]} events')
        title = 'Cluster-time spread' if source == 'cluster' else 'Highest-energy-rechit time spread'
        fig.text(.12, .94, title, fontsize=16, va='top')
        if label:
            fig.text(.12, .895, label, fontsize=11)
        ax.set_xlabel(r'$t_{\mathrm{latest}}-t_{\mathrm{earliest}}$ [ns]')
        ax.set_ylabel('Events')
        ax.set_xlim(page['bin_edges_ns'][0], page['bin_edges_ns'][-1])
        peak = max((max(s['bin_counts'], default=0) for s in page['algorithms'].values()), default=0)
        ax.set_ylim(0, max(1., peak*1.3))
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        ax.minorticks_on()
        ax.tick_params(which='both', direction='in', top=True, right=True)
        ax.legend(loc='upper right', frameon=False)
        fig.text(.12, .085, 'All input events; selection independent for each algorithm; at least two valid times.', fontsize=10)
        detail = ('Times: hcal_time. All valid cluster times enter the event maximum and minimum.'
                  if source == 'cluster' else
                  'Times: hbhe_pfrh_time of each highest-energy rechit. Invalid seed times have no fallback.')
        fig.text(.12, .055, detail, fontsize=9)
        fig.text(.98, .055, f'{page_number}/{total}', ha='right', fontsize=9)
        pdf.savefig(fig)
        plt.close(fig)


def angle(phi):
    return math.atan2(math.sin(phi), math.cos(phi))


def coordinate(eta, phi, ev, mode):
    dphi = angle(phi - ev.phi0)
    return math.hypot(eta - ev.eta0, dphi) if mode == 'dr' else dphi


def load_events(path, tree, validate_seeds=True):
    names = ['run', 'lumi', 'event', 'gen_eta', 'gen_phi', 'gen_pdgId',
             'gen_vx', 'gen_vy', 'gen_vz', 'clusterIdx', 'hcal_eta',
             'hcal_phi', 'hcal_time', 'hcal_energy',
             'hcal_seed_eta', 'hcal_seed_phi', 'hcal_seed_depth']
    names += ['hbhe_pfrh_' + x for x in ['eta','phi','depth','time','energy',
              'fracInCluster','energyFracInCluster','clusterIdx','ieta','iphi']]
    with uproot.open(path) as f:
        t = f[tree]
        missing = sorted(set(names) - set(t.keys()))
        if missing:
            raise ValueError(f'{path}: missing branches: {missing}')
        rows = ak.to_list(t.arrays(names))
    events = {}
    for r in rows:
        key = tuple(int(r[k]) for k in ['run','lumi','event'])
        if key in events:
            raise ValueError(f'Duplicate event key {key}')
        pion_indices = [i for i,pid in enumerate(r['gen_pdgId']) if abs(pid)==211]
        if not pion_indices:
            raise ValueError(f'{key}: no generated pion for the reference')
        ref = pion_indices[0]
        gen = [(r['gen_eta'][i],r['gen_phi'][i],r['gen_vx'][i],r['gen_vy'][i],r['gen_vz'][i]) for i in pion_indices]
        cs = [Cluster(int(cid), r['hcal_eta'][i],r['hcal_phi'][i],r['hcal_time'][i],r['hcal_energy'][i],[],
                      (r['hcal_seed_eta'][i],r['hcal_seed_phi'][i],int(r['hcal_seed_depth'][i]))) for i,cid in enumerate(r['clusterIdx'])]
        lookup = {c.cid:c for c in cs}
        if len(lookup) != len(cs):
            raise ValueError(f'{key}: duplicate cluster IDs')
        cols = ['eta','phi','depth','time','energy','fracInCluster','energyFracInCluster','clusterIdx','ieta','iphi']
        arrays = [r['hbhe_pfrh_'+c] for c in cols]
        if len({len(a) for a in arrays}) != 1:
            raise ValueError(f'{key}: mismatched PFRecHit branch lengths')
        for eta,phi,depth,time,energy,frac,assigned,cid,ieta,iphi in zip(*arrays):
            if cid < 0 or frac <= 0:
                continue
            if int(cid) not in lookup:
                raise ValueError(f'{key}: unknown associated cluster {cid}')
            # Signed ieta, iphi, depth identify HB/HE cells in these barrel samples.
            hit = Hit((int(ieta),int(iphi),int(depth)),eta,phi,int(depth),time,energy,frac,max(0.,assigned))
            lookup[int(cid)].hits.append(hit)
        for c in cs:
            # Depth merging can repeat one cell with complementary fractions.
            # Combine its assigned energies and fractions before plotting/counting.
            merged = {}
            for h in c.hits:
                if h.key in merged:
                    old = merged[h.key]
                    if not np.allclose([old.eta,old.phi,old.energy,old.time],
                                       [h.eta,h.phi,h.energy,h.time],rtol=1e-5,atol=1e-5,equal_nan=True):
                        raise ValueError(f'{key}: inconsistent duplicate cell {h.key}')
                    old.fraction += h.fraction
                    old.assigned += h.assigned
                else:
                    merged[h.key] = h
            c.hits = list(merged.values())
            if validate_seeds:
                seed_hit(c)  # Display marker uses the actual stored algorithm seed.
        events[key] = Event(key,r['gen_eta'][ref],r['gen_phi'][ref],gen,cs)
    return events


def assigned_map(c):
    return {h.key:h.assigned for h in c.hits}


def moments(c, min_fraction):
    hs = [h for h in c.hits if valid(h.time) and h.fraction >= min_fraction and h.assigned>0]
    if not hs:
        return None
    t = np.array([h.time for h in hs]); w = np.array([h.assigned for h in hs])
    mu = float(np.average(t, weights=w))
    return mu, float(np.average((t-mu)**2,weights=w)), float(w.sum()), len(hs)


def invalid_fractions(ev):
    # Merge shared cells so a cell counted in several clusters is counted once.
    cells = {}
    for c in ev.clusters:
        for h in c.hits:
            if h.key not in cells:
                cells[h.key] = [h.time, 0.]
            cells[h.key][1] += h.assigned
    if not cells:
        return 1.,1.
    bad = sum(not valid(t) for t,w in cells.values()) / len(cells)
    total = sum(w for t,w in cells.values())
    bad_energy = sum(w for t,w in cells.values() if not valid(t)) / max(total,1e-12)
    return bad,bad_energy


def evaluate(std, timing, args):
    """Require an allowed multiplicity and enough distinct cells in EVERY cluster.

    Timing contrast and invalid fractions rank candidates; there is no hidden
    energy/substantive-membership cut and no requirement for just one child pair.
    Angular assignment is a diagnostic, not a simulation truth association.
    """
    n = len(timing.gen)
    bad_s, bade_s = invalid_fractions(std)
    bad_t, bade_t = invalid_fractions(timing)
    bad, bad_e = max(bad_s, bad_t), max(bade_s, bade_t)
    counts = [sum(h.fraction >= args.min_hit_fraction for h in c.hits)
              for c in timing.clusters]
    valid_counts = [sum(valid(h.time) and h.fraction >= args.min_hit_fraction
                        for h in c.hits) for c in timing.clusters]
    record = dict(run=std.key[0], lumi=std.key[1], event=std.key[2],
                  generated_pion_count=n, standard_cluster_count=len(std.clusters),
                  timing_cluster_count=len(timing.clusters),
                  timing_cluster_ids=[c.cid for c in timing.clusters],
                  timing_rechit_counts=counts, timing_valid_rechit_counts=valid_counts,
                  timing_cluster_energies_GeV=[c.energy for c in timing.clusters],
                  timing_invalid_cluster_times=sum(not valid(c.time) for c in timing.clusters),
                  invalid_hit_fraction=bad, invalid_energy_fraction=bad_e,
                  qualified=False, score=0., rejection_reasons=[])

    # One-to-one direction assignment, minimizing total squared delta R.
    # No generated time is present in these files, so no time truth is invented.
    matches = []
    if timing.clusters:
        distances = np.array([[math.hypot(c.eta-g[0], angle(c.phi-g[1]))
                               for g in timing.gen] for c in timing.clusters])
        rows, cols = linear_sum_assignment(distances**2)
        matches = [dict(cluster_id=timing.clusters[i].cid, pion_index=int(j),
                        delta_r=float(distances[i,j])) for i,j in zip(rows,cols)]
    record['angular_assignment'] = matches
    record['max_assigned_delta_r'] = max((m['delta_r'] for m in matches), default=None)

    def pooled(ev):
        ms = [m for c in ev.clusters if (m := moments(c, args.min_hit_fraction)) is not None]
        energy = sum(m[2] for m in ms)
        return sum(m[1]*m[2] for m in ms)/energy if energy > 0 else None
    sv, tv = pooled(std), pooled(timing)
    improvement = (1-tv/sv) if sv is not None and sv > 0 and tv is not None else None
    record['standard_pooled_time_variance_ns2'] = sv
    record['timing_pooled_time_variance_ns2'] = tv
    record['spread_improvement'] = improvement
    record['timing_mean_hit_times_ns'] = [m[0] if (m := moments(c, 0.)) else None
                                         for c in timing.clusters]
    links = {str(c.cid): matched_children(c, timing, args) for c in std.clusters}
    record['standard_to_timing_shared_energy_links'] = links
    merge_excess = sum(max(0, len(children)-1) for children in links.values())
    reasons = record['rejection_reasons']
    allowed = allowed_cluster_counts(args, n)
    record['expected_cluster_count'] = args.expected_clusters if args.expected_clusters is not None else n
    record['allowed_timing_cluster_counts'] = sorted(set(allowed))
    if len(timing.clusters) not in allowed:
        reasons.append('timing cluster count is outside the requested counts')
    if any(count < args.min_rechits for count in counts):
        reasons.append('at least one timing cluster has too few associated rechits')
    if any(count < args.min_valid_hits for count in valid_counts):
        reasons.append('at least one timing cluster has too few valid timed rechits')
    if bad > args.max_invalid_fraction or bad_e > args.max_invalid_energy_fraction:
        reasons.append('event invalid-time fraction exceeds configured limit')
    if args.max_match_dr is not None and (len(timing.clusters) != n or len(matches) != n or
            record['max_assigned_delta_r'] > args.max_match_dr):
        reasons.append('one-to-one angular assignment exceeds max-match-dr')
    record['qualified'] = not reasons
    record['reason'] = '; '.join(reasons) if reasons else 'qualified by multiplicity and rechit cuts'
    if record['qualified']:
        quality = (1-bad)**2 * (1-bad_e)**2
        count_contrast = max(0, n-len(std.clusters))/n
        spread_gain = max(0., improvement) if improvement is not None else 0.
        record['score'] = quality * (1+count_contrast) * (1+spread_gain) * (1+merge_excess/n)
    return record


def parse_event_request(tokens, keys):
    selected=[]
    for token in tokens:
        for part in token.split(','):
            try:
                fields=tuple(map(int,part.split(':')))
            except ValueError as e:
                raise ValueError(f'Invalid event identifier {part!r}') from e
            matches=[k for k in keys if (k[2]==fields[0] if len(fields)==1 else k==fields)]
            if len(fields) not in (1,3) or len(matches)!=1:
                raise ValueError(f'{part!r}: expected one matched event; use run:lumi:event if ambiguous')
            if matches[0] not in selected:
                selected.append(matches[0])
    return selected


def representative_hits(c):
    return [max((h for h in c.hits if h.depth==depth),key=lambda h:(h.energy,h.assigned))
            for depth in sorted({h.depth for h in c.hits})]


def seed_hit(c):
    eta,phi,depth=c.seed
    matches=[h for h in c.hits if h.depth==depth and abs(h.eta-eta)<1e-5 and abs(angle(h.phi-phi))<1e-5]
    if len(matches)!=1:
        raise ValueError(f'Cluster {c.cid}: expected one cell matching the stored seed, found {len(matches)}')
    return matches[0]


def identical_clusters(a,b):
    ah={h.key:h for h in a.hits};bh={h.key:h for h in b.hits}
    return ah.keys()==bh.keys() and np.allclose([a.eta,a.phi,a.time,a.energy],
        [b.eta,b.phi,b.time,b.energy],atol=1e-6,rtol=1e-5,equal_nan=True) and all(
        np.allclose([ah[k].fraction,ah[k].assigned],[bh[k].fraction,bh[k].assigned],
                    atol=1e-7,rtol=1e-5) for k in ah)


def fit_envelope(points, yrange, trange, coverage=1., preferred=None):
    """Symbolic display ellipse, independent of physical cluster boundaries."""
    if not points:
        return None
    p=np.asarray(points,float)
    med=np.median(p/np.array([trange,yrange]),axis=0)
    order=np.argsort(np.sum((p/np.array([trange,yrange])-med)**2,axis=1))
    q=p[order[:max(1,math.ceil(coverage*len(p)))]]
    center=(q.min(axis=0)+q.max(axis=0))/2
    radii=np.maximum(np.ptp(q,axis=0)/2,[.025*trange,.04*yrange])
    radii*=max(1.,np.sqrt(np.max(np.sum(((q-center)/radii)**2,axis=1))))
    radii+=np.array([.012*trange,.015*yrange])
    if preferred is not None:
        # Prefer the stored center only when it does not inflate the footprint.
        # The center is not a member mark and no longer has to be enclosed.
        basearea=np.prod(radii)
        for blend in [1.,.75,.5,.25]:
            cc=(1-blend)*center+blend*np.asarray(preferred)
            rr=np.maximum(np.max(np.abs(q-cc),axis=0),[.015*trange,.025*yrange])
            rr*=max(1.,np.sqrt(np.max(np.sum(((q-cc)/rr)**2,axis=1))))
            rr+=np.array([.012*trange,.015*yrange])
            if np.prod(rr)<=1.20*basearea:
                center,radii=cc,rr
                break
    return center,radii


def ellipse_edge(spec):
    if spec is None: return []
    c,r=spec
    theta=np.linspace(0,2*np.pi,181)
    return list(c+np.column_stack([np.cos(theta),np.sin(theta)])*r)


def matched_children(source,timing,args):
    sm=assigned_map(source);se=sum(sm.values())
    matches=[]
    for child in timing.clusters:
        cm=assigned_map(child);ce=sum(cm.values())
        overlap=sum(min(w,cm.get(k,0.)) for k,w in sm.items())
        if (overlap > 0 and overlap/max(se,1e-12)>=args.min_child_fraction and
            overlap/max(ce,1e-12)>=args.min_child_purity):
            matches.append(child.cid)
    return matches


def draw_page(pdf,std,timing,record,page,total,args):
    # Invalid-time cluster centers get their own categorical lane and outline.
    fig,axes=plt.subplots(1,2,figsize=(13,7.3))
    fig.subplots_adjust(left=.075,right=.98,bottom=.16,top=.76,wspace=.20)
    palette=cluster_palette(timing.clusters)
    cell_clusters={}
    cell_depths={}
    for c in timing.clusters:
        for h in c.hits:
            cell_depths.setdefault(h.key,set()).add(c.seed[2])
            cell_clusters.setdefault(h.key,set()).add(c.cid)
    times=[o.time for ev in [std,timing] for c in ev.clusters for o in [c,*c.hits] if valid(o.time)]
    vals=[coordinate(o.eta,o.phi,std,args.coordinate) for ev in [std,timing] for c in ev.clusters for o in [c,*c.hits]]
    baselo=min(-.6,min(times,default=0)-1.)
    basehi=max(6.,max(times,default=5)+1.)
    baseylo=min(vals,default=0);baseyhi=max(vals,default=.3)
    yrange=max(.25,baseyhi-baseylo);trange=basehi-baselo
    baseylo-=.10*yrange;baseyhi+=.10*yrange
    # None signifies categorical invalid time, never a physical -999 coordinate.
    def rawxy(o):
        return (o.time if valid(o.time) else None,coordinate(o.eta,o.phi,std,args.coordinate))
    invalid_keys=[(ai,c.cid) for ai,ev in enumerate([std,timing]) for c in ev.clusters if not valid(c.time)]
    same={t.cid:s.cid for t in timing.clusters for s in std.clusters if identical_clusters(s,t)}
    invalid_group={key:((0,same[key[1]]) if key[0]==1 and key[1] in same else key) for key in invalid_keys}
    # One canonical mark per detector cell, shared by both algorithms. Prefer
    # its timingPF invalid-cluster lane, then a standardPF lane, then the general
    # invalid column. Membership colors retain every parent, regardless of lane.
    cells={};cell_parents={};seed_cells=set()
    for ai,ev in enumerate([std,timing]):
        for c in ev.clusters:
            seed_cells.add(seed_hit(c).key)
            for h in c.hits:
                cells.setdefault(h.key,h)
                cell_parents.setdefault(h.key,set()).add((ai,c.cid))
    cell_slots={};lane_conflicts=[]
    for cell,h in cells.items():
        candidates=[k for k in cell_parents[cell] if k in invalid_group] if not valid(h.time) else []
        ordered=sorted(candidates,key=lambda k:(-k[0],k[1]))
        cell_slots[cell]=invalid_group[ordered[0]] if ordered else None
        # A unique dot cannot be enclosed by two disjoint same-algorithm lanes.
        # Reject that layout instead of duplicating the cell or dropping a parent.
        if len({invalid_group[k] for k in candidates})>1:
            lane_conflicts.append(dict(cell=list(cell),reason='shared cell requires distinct invalid-cluster lanes'))
    centers={};points=[{},{}];glyphs={};glyph_clusters={};depth_positions={}
    for ai,ev in enumerate([std,timing]):
        for c in ev.clusters:
            key=(ai,c.cid);centers[key]=rawxy(c)
            points[0][key]=[rawxy(h) for h in c.hits];points[1][key]=[]
            for h in representative_hits(c):
                hx,hy=rawxy(h)
                slot=cell_slots[h.key] if hx is None else None
                pos=h.key
                if pos not in depth_positions:
                    # Separate touching 1 ns bars, including adjacent integer
                    # times. Reuse exactly coincident shared seeds across algos.
                    offsets=[0.]+[sign*n*.047*yrange for n in range(1,100) for sign in [1,-1]]
                    for offset in offsets:
                        yy=hy+offset
                        collision=any(
                            ((hx is None and ox is None and slot==oslot) or
                             (hx is not None and ox is not None and abs(hx-ox)<=1.05))
                            and abs(yy-oy)<.032*yrange
                            for ox,oy,depth,raw,oslot in glyphs.values())
                        if not collision:break
                    depth_positions[pos]=yy
                yy=depth_positions[pos]
                points[1][key].append((hx,yy))
                glyphs[h.key]=(hx,yy,h.depth,hy,slot)
                glyph_clusters.setdefault(h.key,set()).update(cell_clusters.get(h.key,set()))
    # Fit against marker/line footprints in display units. Iterate because padding
    # in points must be converted with the final axes scale. All member marks
    # participate, even if --contour-fraction was set below one in older commands.
    physical_left=baselo;right=basehi;ylow=baseylo;yhigh=baseyhi
    specs=[{},{}];stripwidth=1.7+2.5*len(set(invalid_group.values()))
    for iteration in range(30):
        separator=physical_left-.30;stripleft=separator-stripwidth
        for ax in axes:ax.set(xlim=(stripleft,right),ylim=(ylow,yhigh))
        ax=axes[0]
        dx=(right-stripleft)/ax.bbox.width*fig.dpi/72
        dy=(yhigh-ylow)/ax.bbox.height*fig.dpi/72
        # Only categorical x positions move: resolve coincident/nearby invalid
        # cells by 9 pt steps, preserving every measured angular coordinate.
        hit_offsets={};occupied=[]
        for cell,h in sorted(cells.items()):
            if valid(h.time):continue
            yy=coordinate(h.eta,h.phi,std,args.coordinate);slot=cell_slots[cell]
            for off in [0]+[sign*9*n for n in range(1,len(cells)+1) for sign in [-1,1]]:
                if all(oslot!=slot or math.hypot(off-ox,(yy-oy)/dy)>=8.5
                       for oslot,ox,oy in occupied):break
            hit_offsets[cell]=off*dx
            occupied.append((slot,off,yy))
        specs=[{},{}];supports=[{},{}]
        def box(x,y,halfx,halfy):
            return [(x-hafx,y-hafy) for hafx in [-halfx,halfx] for hafy in [-halfy,halfy]]
        for j in [0,1]:
            for ai,ev in enumerate([std,timing]):
                for c in ev.clusters:
                    key=(ai,c.cid);support=[]
                    for hx,hy in points[j][key]:
                        if hx is None:continue
                        # Seed endpoints ±0.5 ns, plus 5 pt clearance. Rechits 5 pt.
                        support+=box(hx,hy,(.5 if j==1 else 0)+5*dx,5*dy)
                        # Reserve space just above a depth seed for its D label.
                        if j==1:support+=box(hx,hy+8*dy,9*dx,6*dy)
                    # Both panels describe the full cluster extent: real valid
                    # member positions plus the angular span of invalid members.
                    valid_members=[rawxy(h) for h in c.hits if valid(h.time)]
                    all_y=[coordinate(h.eta,h.phi,std,args.coordinate) for h in c.hits]
                    if valid_members:
                        for hx,hy in valid_members:support+=box(hx,hy,5*dx,5*dy)
                        tx=(min(x for x,y in valid_members)+max(x for x,y in valid_members))/2
                        for yy in [min(all_y),max(all_y)]:support+=box(tx,yy,5*dx,5*dy)
                    cx,cy=centers[key]
                    supports[j][key]=support
                    specs[j][key]=fit_envelope(support,yrange,trange,1.,
                        preferred=(cx,cy) if cx is not None else None)
            for c in std.clusters:
                support=list(supports[j][(0,c.cid)])
                for cid in matched_children(c,timing,args):
                    for ex,ey in ellipse_edge(specs[j][(1,cid)]):
                        support+=box(ex,ey,4*dx,4*dy)
                specs[j][(0,c.cid)]=fit_envelope(support,yrange,trange,1.)
        invalid_specs=[{},{}];invalid_supports=[{},{}]
        for j in [0,1]:
            for key in invalid_keys:
                support=[]
                c=next(c for c in [std,timing][key[0]].clusters if c.cid==key[1])
                if j==0:
                    for h in c.hits:
                        if not valid(h.time):
                            support+=box(hit_offsets[h.key],rawxy(h)[1],5*dx,5*dy)
                else:
                    for hx,hy in points[j][key]:
                        if hx is None:
                            support+=box(0,hy,.5+5*dx,5*dy)
                            support+=box(0,hy+8*dy,9*dx,6*dy)
                # All-invalid clusters also retain their full angular range in
                # panel (b), even if only one depth representative is displayed.
                c=next(c for c in [std,timing][key[0]].clusters if c.cid==key[1])
                for h in c.hits:
                    if not valid(h.time):support+=box(0,coordinate(h.eta,h.phi,std,args.coordinate),5*dx,5*dy)
                invalid_supports[j][key]=support
                invalid_specs[j][key]=fit_envelope(support,yrange,2.,1.,preferred=(0,centers[key][1]))
            for tid,sid in same.items():
                if (1,tid) not in invalid_specs[j]:continue
                support=list(invalid_supports[j][(0,sid)])
                for ex,ey in ellipse_edge(invalid_specs[j][(1,tid)]):
                    support+=box(ex,ey,4*dx,4*dy)
                invalid_specs[j][(0,sid)]=fit_envelope(support,yrange,2.,1.)
        lane_radii={group:max(invalid_specs[j][key][1][0] for key in invalid_keys
                              if invalid_group[key]==group for j in [0,1]) for group in set(invalid_group.values())}
        newstrip=1.7+sum(2*r+.35 for r in lane_radii.values())
        allspec=[s for panel in specs for s in panel.values() if s is not None]
        yspec=allspec+[s for panel in invalid_specs for s in panel.values()]
        new=np.array([min([baselo]+[c[0]-r[0]-.15 for c,r in allspec]),
                      max([basehi]+[c[0]+r[0]+.15 for c,r in allspec]),
                      min([baseylo]+[c[1]-r[1]-.012*yrange for c,r in yspec]),
                      max([baseyhi]+[c[1]+r[1]+.012*yrange for c,r in yspec])])
        old=np.array([physical_left,right,ylow,yhigh])
        physical_left,right,ylow,yhigh=new
        converged=np.max(np.abs(new-old)/np.array([trange,trange,yrange,yrange]))<1e-5 and abs(newstrip-stripwidth)<1e-5
        stripwidth=newstrip
        if converged:break
    separator=physical_left-.30;stripleft=separator-stripwidth;invalid_x=separator-.85
    invalid_columns={};cursor=separator-1.7
    for key in sorted(lane_radii):
        radius=lane_radii[key];invalid_columns[key]=cursor-.175-radius
        cursor-=2*radius+.35
    def displayed(pt,key=None):return (invalid_columns.get(invalid_group.get(key,key),invalid_x) if pt[0] is None else pt[0],pt[1])
    yname=r'$\Delta R$ from first gen pion' if args.coordinate=='dr' else r'$\Delta\phi$ from first gen pion [rad]'
    stars=[];annotations=[]
    for j,ax in enumerate(axes):
        ax.set(xlim=(stripleft,right),ylim=(ylow,yhigh),xlabel=r'$t$ [ns]',ylabel=yname)
        ax.set_title(['(a) Clusters and rechits','(b) Clusters and depth seeds'][j],fontsize=12,pad=12)
        ax.axvspan(stripleft,separator,color='0.96',zorder=0)
        ax.axvline(separator,color='0.65',ls=':',lw=.7)
        ticks=np.arange(math.ceil(baselo/2)*2,right+.01,2.)
        ax.set_xticks([(stripleft+separator)/2,*ticks],['invalid\n(-999)',*[f'{v:g}' for v in ticks]])
        ax.tick_params(direction='in',top=True,right=True,labelsize=10)
        for spine in ax.spines.values():spine.set_linewidth(.8)
        # Each all-invalid cluster has one contour in its own categorical lane.
        for ai,ev in enumerate([std,timing]):
            for c in ev.clusters:
                spec=specs[j][(ai,c.cid)]
                if spec:
                    center,radii=spec
                    ax.add_patch(Ellipse(center,2*radii[0],2*radii[1],fill=False,edgecolor='black',
                                         ls='--' if ai==0 else '-',lw=1.1,zorder=2))
        for key,spec in invalid_specs[j].items():
            center,radii=spec;center=center+np.array([invalid_columns[invalid_group[key]],0.])
            ax.add_patch(Ellipse(center,2*radii[0],2*radii[1],fill=False,edgecolor='black',
                                 ls='--' if key[0]==0 else '-',lw=1.1,zorder=2))
        if j==0:
            for cell,h in cells.items():
                x,y=rawxy(h)
                if x is None:x=invalid_columns.get(cell_slots[cell],invalid_x)+hit_offsets[cell]
                colored_hit(ax,x,y,[palette[cid] for cid in sorted(cell_clusters.get(cell,set()))],cell in seed_cells)
    for glyphkey,(hx,hy,depth,raw,slot) in glyphs.items():
        x=invalid_columns.get(slot,invalid_x) if hx is None else hx
        # Exactly one ns on the physical axis; same visual width in the invalid
        # category (where width has no physical timing interpretation).
        axes[1].plot([x-.5,x+.5],[hy,hy],color='white',lw=3.5,zorder=4,solid_capstyle='butt')
        colors=[palette[cid] for cid in sorted(glyph_clusters[glyphkey])] or ['0.65']
        for i,color in enumerate(colors):
            axes[1].plot([x-.5+i/len(colors),x-.5+(i+1)/len(colors)],[hy,hy],
                         color=color,lw=2.,zorder=5,solid_capstyle='butt')
        if hy!=raw:axes[1].plot([x,x],[raw,hy],color='0.6',lw=.5,zorder=3)
        annotations.append(axes[1].annotate(f'D{depth}',(x,hy),xytext=(0,6),textcoords='offset points',
                                           ha='center',fontsize=8,zorder=8))
    fig.canvas.draw();renderer=fig.canvas.get_renderer();ax=axes[1]
    # Labels must avoid every depth bar, including bars in the invalid strip.
    placed=[]
    for hx,hy,depth,raw,slot in glyphs.values():
        x=invalid_columns.get(slot,invalid_x) if hx is None else hx
        left,yy=ax.transData.transform((x-.5,hy))
        right,_=ax.transData.transform((x+.5,hy))
        margin=2*fig.dpi/72
        placed.append(Bbox.from_extents(left-margin,yy-margin,right+margin,yy+margin))
    for pos in stars:
        px,py=ax.transData.transform(pos);placed.append(Bbox.from_bounds(px-7,py-7,14,14))
    for label in annotations:
        start=label.get_position();best=None
        for dxp,dyp in [(0,0),(0,8),(0,-19),(12,5),(-12,5),(0,18),(0,-27)]:
            label.set_position((start[0]+dxp,start[1]+dyp));label.update_positions(renderer)
            bounds=Text.get_window_extent(label,renderer).expanded(1.1,1.1)
            if not (bounds.x0>ax.bbox.x0 and bounds.x1<ax.bbox.x1 and bounds.y0>ax.bbox.y0 and bounds.y1<ax.bbox.y1):continue
            cost=10000*sum(bounds.overlaps(other) for other in placed)+dxp*dxp+dyp*dyp
            if best is None or cost<best[0]:best=(cost,label.get_position(),bounds)
        if best:label.set_position(best[1]);placed.append(best[2])
        else:label.set_position(start)
    fig.text(.075,.946,f'Event {std.key[2]}  (run {std.key[0]}, lumi {std.key[1]})',fontsize=13)
    nbad=sum(not valid(c.time) for c in timing.clusters)
    fig.text(.075,.908,f'{len(std.gen)} generated pions; {len(std.clusters)} standardPF; {len(timing.clusters)} seedTimingPF clusters; {nbad} invalid cluster times',fontsize=10)
    handles=[Line2D([],[],color='black',ls='--',label='standardPF'),
             Line2D([],[],color='black',ls='-',label='seedTimingPF'),
             Line2D([],[],color='black',ls='',marker='o',ms=5,label='Large dot: cluster seed (a)')]
    fig.legend(handles=handles,loc='upper right',bbox_to_anchor=(.98,.972),frameon=False,fontsize=10)
    depths=sorted({c.seed[2] for c in timing.clusters})
    colorhandles=[tuple(Line2D([],[],color=palette[c.cid],marker='o',ls='',ms=5)
                        for c in sorted(timing.clusters,key=lambda c:c.cid) if c.seed[2]==d) for d in depths]
    fig.text(.075,.848,'TimingPF cluster seed depth:',fontsize=10)
    fig.legend(handles=colorhandles,labels=[f'D{d}' for d in depths],
               handler_map={tuple:HandlerTuple(ndivide=None,pad=.3)},loc='center left',bbox_to_anchor=(.31,.854),ncol=len(depths),
               frameon=False,fontsize=10,handletextpad=.4,columnspacing=1.3)
    fig.text(.075,.815,'Shades distinguish clusters; split colors indicate shared membership, not energy fractions.',fontsize=8,color='0.3')
    fig.text(.075,.055,f'Invalid rechit times: {record["invalid_hit_fraction"]:.0%}.  Outlines: all-hit angular span, valid-hit time span.  Depth bars: 1 ns.',fontsize=8)
    fig.text(.98,.055,f'{page}/{total}',ha='right',fontsize=9)
    # Validate actual artist footprint containment, not just point centers.
    # The assertions also protect against accidentally reverting to center-only fits.
    for j in [0,1]:
        for key,support in supports[j].items():
            spec=specs[j][key]
            if spec and support:
                center,radii=spec
                assert np.max(np.sum(((np.asarray(support)-center)/radii)**2,axis=1))<1.00001
        for key,support in invalid_supports[j].items():
            center,radii=invalid_specs[j][key]
            assert np.max(np.sum(((np.asarray(support)-center)/radii)**2,axis=1))<1.00001
    # Require a real gap in display points. Same-algorithm regions must also
    # have disjoint interiors; cross-algorithm containment is intentional.
    collisions=list(lane_conflicts)
    for j,ax in enumerate(axes):
        outlines=[]
        for key,spec in specs[j].items():
            if spec:outlines.append((key,np.asarray(ellipse_edge(spec))))
        for key,spec in invalid_specs[j].items():
            edge=np.asarray(ellipse_edge(spec))+[invalid_columns[invalid_group[key]],0]
            outlines.append((key,edge))
        polygons=[(key,Polygon(ax.transData.transform(edge)*72/fig.dpi)) for key,edge in outlines]
        for i,(key,a) in enumerate(polygons):
            for other,b in polygons[i+1:]:
                same_algo=key[0]==other[0] and key!=other
                if (same_algo and a.intersects(b)) or a.boundary.distance(b.boundary)<args.outline_gap_pt:
                    collisions.append(dict(panel=j,clusters=[list(key),list(other)],reason='overlap or insufficient outline gap'))
    result=dict(clear=not collisions,collisions=collisions,
                timing_cluster_seed_depths={str(c.cid):c.seed[2] for c in timing.clusters},
                timing_cluster_colors={str(cid):color for cid,color in palette.items()},
                unique_rechit_mark_count=len(cells),unique_depth_seed_mark_count=len(glyphs),
                rechit_marks=[dict(cell=list(cell),timing_clusters=sorted(cell_clusters.get(cell,set())),
                                  invalid_lane=list(cell_slots[cell]) if cell_slots[cell] is not None else None,
                                  categorical_offset_pt=hit_offsets.get(cell,0)/dx)
                             for cell in sorted(cells)],
                shared_cells_with_multiple_timing_clusters=sum(len(cs)>1 for cs in cell_clusters.values()),
                angular_extent_uses_all_members=True,time_extent_uses_valid_members=True,
                shared_cells_with_multiple_seed_depths=sum(len(ds)>1 for ds in cell_depths.values()),
                identical_cluster_pairs=[dict(standard=s,timing=t) for t,s in same.items()])
    if pdf is not None:
        if collisions and not args.allow_overlapping_outlines:
            plt.close(fig)
            raise ValueError(f'Event {std.key[2]} cannot be displayed with separated outlines at unchanged coordinates. See layout audit.')
        pdf.savefig(fig)
    plt.close(fig)
    return result


def build_parser():
    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--input-dir',type=Path,default=Path('/eos/user/c/chtong/Public/Rereco/MultiPion_Oct_22/MultiPion_eta0.2_4ns_n1000'),help='Directory for automatically constructed ROOT filenames')
    p.add_argument('--standard',type=Path,help='Override the standardPF ROOT path; default is built from N and input-dir')
    p.add_argument('--timing',type=Path,help='Override the seedTimingPF ROOT path; default is built from N and input-dir')
    p.add_argument('--output',type=Path,help='Override the PDF filename; default is /eos/user/c/chtong/Public/Rereco/MultiPion_Oct_22/MultiPion_eta0.2_4ns_n1000/multipion_N{N}.pdf')
    p.add_argument('--sample-label',help='Optional sample description on the histogram pages')
    p.add_argument('--tree',default='pfObjectsNtupler/pfTree')
    p.add_argument('--mode',choices=['auto','manual'],default='auto')
    p.add_argument('--events',nargs='+',help='Event numbers or run:lumi:event; spaces or commas. Implies manual mode.')
    p.add_argument('--max-events',type=int,default=100,help='Auto mode only; never limits manual mode')
    p.add_argument('--coordinate',choices=['dr','dphi'],default='dr')
    p.add_argument('--min-rechits',type=int,default=2,help='Minimum unique associated cells in EVERY timing cluster')
    p.add_argument('--expected-clusters','--N',type=int,default=4,help='Reference count N for display selection and automatic input/output filenames')
    count_options=p.add_mutually_exclusive_group()
    count_options.add_argument('--cluster-count-mode',choices=['exact','nearby','range'],default='nearby',help='Display timingPF count: exact=N, nearby=max(2,N-2) through N+2, range=2 through N+1')
    count_options.add_argument('--timing-clusters',type=int,nargs='+',help='Explicit allowed timingPF counts, e.g. 2 3 4; replaces the count mode')
    count_options.add_argument('--timing-count-offsets',type=int,nargs='+',help='Custom offsets from expected count (or generated pion count); replaces the count mode')
    p.add_argument('--selection-order',choices=['balanced','score'],default='balanced',help='Balanced cycles through multiplicity groups in score order to show each available count')
    p.add_argument('--outline-gap-pt',type=float,default=3.,help='Minimum gap between outlines; automatic selection excludes events without a clear layout')
    overlap_options=p.add_mutually_exclusive_group()
    overlap_options.add_argument('--overlap',choices=['allow','forbid'],default='forbid',help='Allow intersecting outlines, or require separated outlines')
    overlap_options.add_argument('--allow-overlapping-outlines',action='store_true',help='Backward-compatible alias for --overlap allow')
    p.add_argument('--min-hit-fraction',type=float,default=0.,help='Optional minimum fractional assignment; zero counts all associated hits')
    p.add_argument('--min-valid-hits',type=int,default=0,help='Optional minimum valid timed hits in EVERY timing cluster; 0 permits invalid-time clusters')
    p.add_argument('--max-match-dr',type=float,help='Optional maximum angular distance for every one-to-one gen/cluster assignment; not a truth-purity test')
    p.add_argument('--min-child-fraction',type=float,default=0.,help='Minimum shared source-energy fraction for graphical source/child matching only')
    p.add_argument('--min-child-purity',type=float,default=.5,help='Minimum child-energy fraction shared with source for graphical matching only')
    p.add_argument('--max-invalid-fraction',type=float,default=.8,help='Maximum event fraction of unique associated cells with invalid time')
    p.add_argument('--max-invalid-energy-fraction',type=float,default=.5)
    p.add_argument('--contour-fraction',type=float,default=1.,help='Legacy option; padded outlines now include all displayed member marks')
    p.add_argument('--rank-only',action='store_true',help='Write selection audit without PDF')
    hist_options=p.add_mutually_exclusive_group()
    hist_options.add_argument('--histograms-only',action='store_true',help='Skip event selection/layout; output only the two all-event timing histograms')
    hist_options.add_argument('--no-timing-histograms',action='store_true',help='Omit the two appended timing histograms')
    p.add_argument('--hist-bin-width',type=float,default=1.,help='Shared bin width in ns for both algorithms; auto range includes all entries')
    return p


def parse_args(argv=None):
    p=build_parser()
    args=p.parse_args(argv)
    args.allow_overlapping_outlines = args.allow_overlapping_outlines or args.overlap=='allow'
    args.overlap = 'allow' if args.allow_overlapping_outlines else 'forbid'
    if args.events: args.mode='manual'
    if args.mode=='manual' and not args.events: p.error('manual mode requires --events')
    if min(args.max_events,args.min_rechits)<1: p.error('event and rechit counts must be positive')
    if args.timing_clusters and min(args.timing_clusters)<1:p.error('timing-clusters must be positive')
    if args.expected_clusters is not None and args.expected_clusters<1:p.error('expected-clusters must be positive')
    n=args.expected_clusters
    stem=f'MultiPi_N{n}_R177.7_eta0.2_E35.0_DT4.0_DR0.1_offsetFirst0.0.root'
    if args.standard is None:args.standard=args.input_dir/f'pfObjectsNtuple_standardPF_{stem}'
    if args.timing is None:args.timing=args.input_dir/f'pfObjectsNtuple_seedTimingPF_{stem}'
    if args.output is None:args.output=Path(f'/eos/user/c/chtong/Public/Rereco/MultiPion_Oct_22/MultiPion_eta0.2_4ns_n1000/multipion_N{n}.pdf')
    if not math.isfinite(args.hist_bin_width) or args.hist_bin_width<=0:p.error('hist-bin-width must be positive and finite')
    if args.histograms_only and (args.events or args.rank_only):p.error('histograms-only cannot be combined with events or rank-only')
    if args.min_valid_hits<0: p.error('min-valid-hits must be nonnegative')
    if args.max_match_dr is not None and args.max_match_dr<=0: p.error('max-match-dr must be positive')
    if args.outline_gap_pt<=0:p.error('outline-gap-pt must be positive')
    if not 0<=args.min_hit_fraction<=1: p.error('min-hit-fraction must be in [0,1]')
    for key in ['min_child_fraction','min_child_purity','max_invalid_fraction','max_invalid_energy_fraction','contour_fraction']:
        if not 0<=float(getattr(args,key))<=1: p.error(f'{key} must be in [0,1]')
    return args


def main(argv=None):
    args=parse_args(argv)
    print('Loading input events...', flush=True)
    std=load_events(args.standard,args.tree,validate_seeds=not args.histograms_only)
    timing=load_events(args.timing,args.tree,validate_seeds=not args.histograms_only)
    print(f'Loaded {len(std)} standardPF and {len(timing)} seedTimingPF events.', flush=True)
    histograms = None if args.no_timing_histograms else histogram_data(std, timing, args.hist_bin_width)
    common=sorted(std.keys() & timing.keys())
    if not common and not args.histograms_only:
        raise ValueError('No common run:lumi:event keys for event displays. Use --histograms-only for independent input samples.')
    if std.keys()!=timing.keys(): print('WARNING: event displays use matched keys; histograms still use ALL events in each file.')
    for k in ([] if args.histograms_only else common):
        if len(std[k].gen)!=len(timing[k].gen) or not np.allclose(std[k].gen,timing[k].gen,atol=1e-5,rtol=1e-6):
            raise ValueError(f'{k}: generator records differ between files')
        def cell_records(ev):
            return {h.key:(h.energy,h.time,h.eta,h.phi) for c in ev.clusters for h in c.hits}
        sc,tc=cell_records(std[k]),cell_records(timing[k])
        for cell in sc.keys() & tc.keys():
            if not np.allclose(sc[cell],tc[cell],rtol=1e-5,atol=1e-5,equal_nan=True):
                raise ValueError(f'{k}: same cell {cell} has inconsistent raw data across algorithms')
    records=[] if args.histograms_only else [evaluate(std[k],timing[k],args) for k in common]
    records.sort(key=lambda r:(-r['score'],r['run'],r['lumi'],r['event']))
    manual_keys=parse_event_request(args.events,common) if args.mode=='manual' else []
    checked=0
    for r in records:
        key=tuple(r[k] for k in ['run','lumi','event'])
        r['passes_physics_cuts']=r['qualified']
        if r['qualified'] or key in manual_keys:
            r['layout']=draw_page(None,std[key],timing[key],r,1,1,args)
            checked+=1
            if checked%25==0:print(f'Checked layouts for {checked} candidate events...', flush=True)
            if not r['layout']['clear'] and not args.allow_overlapping_outlines:
                r['qualified']=False
                r['rejection_reasons'].append('outlines overlap or touch at unchanged coordinates')
                r['reason']='; '.join(r['rejection_reasons'])
    qualified=[r for r in records if r['qualified']]
    if args.selection_order=='balanced':
        groups={n:[r for r in qualified if r['timing_cluster_count']==n]
                for n in sorted({r['timing_cluster_count'] for r in qualified})}
        qualified=[]
        while any(groups.values()):
            for group in groups.values():
                if group:qualified.append(group.pop(0))
    selected=(manual_keys if args.mode=='manual' else
              [tuple(r[k] for k in ['run','lumi','event']) for r in qualified][:args.max_events])
    args.output.parent.mkdir(parents=True,exist_ok=True)
    audit=args.output.with_suffix('.selection.json')
    payload=dict(selection_definition=('Allowed multiplicities and minimum unique cells in every timing cluster; overlapping outlines permitted' if args.allow_overlapping_outlines else 'Allowed multiplicities, minimum unique cells in every timing cluster, and separated schematic outlines at unchanged physical coordinates'),
                 parameters={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},
                 common_events=len(common),standard_only=sorted(std.keys()-timing.keys()),
                 timing_only=sorted(timing.keys()-std.keys()),selected=selected,ranking=records)
    payload['histograms_only']=args.histograms_only
    payload['input_event_counts']=dict(standardPF=len(std),seedTimingPF=len(timing))
    if histograms is not None:
        hist_audit=args.output.with_suffix('.timing_histograms.json')
        hist_audit.write_text(json.dumps(histograms,indent=2)+'\n')
        payload['timing_histograms_file']=str(hist_audit)
        for source,page in histograms['pages'].items():
            print(f'{source} histogram entries: '+', '.join(f'{name}={s["event_count"]}' for name,s in page['algorithms'].items()))
    audit.write_text(json.dumps(payload,indent=2)+'\n')
    print(f'{sum(r["qualified"] for r in records)} qualifying events / {len(common)} matched events')
    print('Selected:',', '.join(str(k[2]) for k in selected) or '(none)')
    if args.rank_only: return
    if not selected and histograms is None:
        raise SystemExit('No events pass the cuts. Audit written; use manual mode or explicitly adjust selection cuts.')
    if not selected and not args.histograms_only:
        print('No event displays pass the cuts; writing the two all-event histogram pages.')
    lookup={tuple(r[k] for k in ['run','lumi','event']):r for r in records}
    # Finish the PDF before writing it to a potentially synchronized output path.
    with tempfile.TemporaryDirectory(prefix='hcal-plots-') as work:
        rendered=Path(work)/'events.pdf'
        with PdfPages(rendered) as pdf:
            pdf.infodict().update(Title='HCAL cluster timing: event displays and time spreads',Subject='standardPF and seedTimingPF comparison')
            total=len(selected)+(2 if histograms is not None else 0)
            for i,k in enumerate(selected,1):
                draw_page(pdf,std[k],timing[k],lookup[k],i,total,args)
            if histograms is not None:
                for i,source in enumerate(('cluster','seed'),len(selected)+1):
                    draw_timing_histogram(pdf,histograms,source,i,total,args.sample_label)
        args.output.write_bytes(rendered.read_bytes())
    print('PDF:',args.output)
    print('Selection audit:',audit)

if __name__=='__main__':
    main()
