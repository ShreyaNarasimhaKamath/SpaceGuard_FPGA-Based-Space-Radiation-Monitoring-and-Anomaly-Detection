"""
SpaceGuard v2 - injection benchmark.

Ground truth: synthetic enhancements of known amplitude / duration / rise time are multiplied
onto the REAL GOES-16 flux (T1, bands 1-2). Natural variability stays in the data, so we also
measure how often each detector alarms on the un-injected record ("natural alarms/day").

Detectors are all causal / streaming:
  MED   : baseline = median of previous 1800 samples, alarm if x > (1+k)*baseline
  EMA   : asymmetric EMA baseline (alpha = 2^-n), alarm if x > (1+k)*baseline,
          optional "freeze" = baseline update gated off while alarm active (max tmax samples)
  LEGACY: old linear residual vs rolling median, threshold 6*1.4826*MAD from first 10 h (clean)
Persistence M samples to raise an alarm, R consecutive non-exceeding samples to release.

Dev days: Jan 1-2, held-out test day: Jan 3 (different random seeds).
"""
from pathlib import Path
import itertools, sys
import numpy as np, pandas as pd, numba as nb
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "raw" / "goes16_mpshi"
OUT = ROOT / "results" / "v2"; OUT.mkdir(parents=True, exist_ok=True)

DEV, TEST = ["20240101", "20240102"], ["20240103"]
BANDS = [0, 1]; TEL = 0
WARM = 1800; M = 10; R = 30
KS = [0.3, 0.5, 0.75, 1.0, 1.5, 2.0]
AMPS = [0.5, 1.0, 2.0, 4.0]; DURS = [60, 300, 1800]; RISES = [10, 120]
COMBOS = list(itertools.product(AMPS, DURS, RISES))
ALPHAS = [(10, 7), (9, 6), (11, 8)]            # (n_up, n_down): alpha = 2^-n
N_EVENTS, REPS, TMAX = 8, 25, 1800


def load(date, band):
    with Dataset(DATA / f"goes16_mpshi_{date}.nc") as nc:
        f = np.ma.filled(nc["DiffElectronFluxes"][:, TEL, band].astype(float), np.nan)
        q = np.ma.filled(nc["DiffElectronFluxDQFs"][:, TEL, band].astype(float), 255)
    ok = np.isfinite(f) & np.isin(q, (0, 16)) & (f >= 0)      # keep degraded-by-low-counts (16); drop 8/24/255
    return np.where(ok, f, np.nan), ok


@nb.njit(cache=True)
def alarm_from_cond(cond, valid, M, R):
    n = len(cond); a = np.zeros(n, np.bool_); run = 0; off = 0; on = False
    for i in range(n):
        if not valid[i]:
            a[i] = on; continue
        if cond[i]:
            run += 1; off = 0
        else:
            run = 0; off += 1
        if (not on) and run >= M: on = True
        elif on and off >= R: on = False
        a[i] = on
    return a


@nb.njit(cache=True)
def ema_alarm(x, valid, a_up, a_dn, k, M, R, freeze, tmax, warm):
    n = len(x); a = np.zeros(n, np.bool_); state = 0.0; init = False
    run = 0; off = 0; on = False; frz = 0; cnt = 0
    for i in range(n):
        if not valid[i]:
            a[i] = on; continue
        xi = x[i]
        if not init:
            state = xi; init = True; cnt = 1; continue
        cond = (cnt >= warm) and (xi > (1.0 + k) * state)
        if cond:
            run += 1; off = 0
        else:
            run = 0; off += 1
        if (not on) and run >= M:
            on = True; frz = 0
        elif on and off >= R:
            on = False
        a[i] = on
        if freeze and on and frz < tmax:
            frz += 1
        else:
            al = a_up if xi > state else a_dn
            state += al * (xi - state)
        cnt += 1
    return a


def med_baseline(x):
    return pd.Series(x).shift(1).rolling(WARM, min_periods=WARM // 2).median().to_numpy()


def profile(n, s, A, D, tr):
    p = np.zeros(n); e = s + tr + D + tr
    p[s:s + tr] = np.linspace(0, 1, tr, endpoint=False)
    p[s + tr:s + tr + D] = 1.0
    p[s + tr + D:e] = np.linspace(1, 0, tr, endpoint=False)
    return A * p, e


def make_events(rng, valid):
    n = len(valid); order = rng.permutation(len(COMBOS)); evs = []
    for k in range(N_EVENTS):
        A, D, tr = COMBOS[order[k % len(COMBOS)]]
        for _ in range(20):
            s = int(3600 + k * 9000 + rng.integers(0, 1800)); e = s + tr + D + tr
            if valid[s - 60:e].mean() > 0.9: break
        evs.append(dict(s=s, e=e, A=A, D=D, tr=tr))
    return evs


def inject(x, evs):
    g = np.zeros(len(x))
    for ev in evs:
        p, _ = profile(len(x), ev["s"], ev["A"], ev["D"], ev["tr"]); g = np.maximum(g, p)
    return x * (1 + g)


def count_intervals(a, start=WARM, gap=60):
    idx = np.flatnonzero(a[start:])
    if len(idx) == 0: return 0
    return int((np.diff(idx) > gap).sum() + 1)


def configs():
    for k in KS: yield f"MED_k{k}", ("med", k)
    for (nu, nd) in ALPHAS:
        for fr in (0, 1):
            for k in KS: yield f"EMA_{nu}_{nd}_{'frz' if fr else 'nof'}_k{k}", ("ema", nu, nd, fr, k)
    yield "LEGACY_lin6sig", ("legacy",)


def run_all(x, valid, thr_legacy, bmed):
    """return dict name -> alarm array"""
    out = {}
    cond_valid = valid
    for name, c in configs():
        if c[0] == "med":
            cond = np.isfinite(bmed) & (x > (1 + c[1]) * bmed) & valid
            cond[:WARM] = False
            out[name] = alarm_from_cond(cond, valid, M, R)
        elif c[0] == "ema":
            _, nu, nd, fr, k = c
            out[name] = ema_alarm(np.nan_to_num(x), valid, 2.0 ** -nu, 2.0 ** -nd, k, M, R, bool(fr), TMAX, WARM)
        else:
            cond = np.isfinite(bmed) & ((x - bmed) > thr_legacy) & valid
            cond[:WARM] = False
            out[name] = alarm_from_cond(cond, valid, M, R)
    return out


def main():
    rows_ev, rows_nat = [], []
    for split, dates, seed0 in (("dev", DEV, 100), ("test", TEST, 900)):
        for date in dates:
            for band in BANDS:
                x, valid = load(date, band)
                bm = med_baseline(x); res = x - bm
                cal = np.isfinite(res) & (np.arange(len(x)) < 36000)
                c = res[cal]; thr = 6 * 1.4826 * np.median(np.abs(c - np.median(c)))
                clean = run_all(x, valid, thr, bm)
                for name, a in clean.items():
                    rows_nat.append(dict(split=split, date=date, band=band + 1, config=name,
                                         natural_alarms=count_intervals(a), duty_pct=100 * a[WARM:].mean()))
                for rep in range(REPS):
                    rng = np.random.default_rng(seed0 + 1000 * band + 10 * int(date[-1]) + rep)
                    evs = make_events(rng, valid); xi = inject(x, evs)
                    bmi = med_baseline(xi)
                    inj = run_all(xi, valid, thr, bmi)
                    for name, a in inj.items():
                        for ev in evs:
                            s, e = ev["s"], ev["e"]
                            conf = bool(a[s - 1]); seg = a[s:e + 1]; det = bool(seg.any())
                            rows_ev.append(dict(split=split, date=date, band=band + 1, config=name, A=ev["A"], D=ev["D"], tr=ev["tr"],
                                                confounded=conf, detected=det, delay_s=(int(np.argmax(seg)) if det else np.nan),
                                                coverage=float(seg.mean())))
                print("done", split, date, band + 1, flush=True)
    ev = pd.DataFrame(rows_ev); nat = pd.DataFrame(rows_nat)
    ev.to_csv(OUT / "injection_events.csv", index=False); nat.to_csv(OUT / "natural_alarms.csv", index=False)

    e = ev[~ev.confounded]
    core = (e.A >= 1) & (e.D >= 300)
    summ = e.groupby(["split", "config"]).apply(lambda g: pd.Series(dict(
        Pd_all=g.detected.mean(),
        Pd_core=g[(g.A >= 1) & (g.D >= 300)].detected.mean(),
        Pd_short60=g[g.D == 60].detected.mean(),
        Pd_amp05=g[g.A == 0.5].detected.mean(),
        med_delay_s=g.delay_s.median(), p90_delay_s=g.delay_s.quantile(.9),
        coverage=g[g.detected].coverage.mean(), n=len(g))), include_groups=False).reset_index()
    n2 = nat.groupby(["split", "config"]).agg(natural_alarms_per_day=("natural_alarms", "mean"), duty_pct=("duty_pct", "mean")).reset_index()
    summ = summ.merge(n2, on=["split", "config"]); summ.to_csv(OUT / "benchmark_summary.csv", index=False)
    pd.set_option("display.width", 220)
    print(summ[summ.split == "dev"].sort_values("natural_alarms_per_day").round(3).to_string(index=False))


if __name__ == "__main__":
    main()
