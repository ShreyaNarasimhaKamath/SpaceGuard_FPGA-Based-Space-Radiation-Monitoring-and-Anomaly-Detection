"""
SpaceGuard v2 - hazard channel (>2 MeV integral electrons, deep-dielectric charging) as a COUNT detector.
Input: integer counts per second, c = round(flux/Q) with Q~1326 pfu per count/s (quantum measured in the L1b data).
Negative (background-subtracted, DQF 8/24) samples -> 0 counts but valid; NaN/fill/255 -> invalid (state held).
Detector (shift-only, one register): L <- L - (L>>SH) + (c<<FR), FR=8 fractional bits (leaky window ~2^SH s; steady state L = 2^(SH+FR) * mean rate)
Alarm: Schmitt trigger, set when L >= TH, clear when L < TL.  TH is an ABSOLUTE hazard level (SWPC >2 MeV alert = 1000 pfu).
Evaluation: (1) real days: alarms/day vs alert level; (2) Poisson count injection (trapezoid) for detect prob and delay.
"""
import numpy as np, pandas as pd, numba as nb
from netCDF4 import Dataset
from pathlib import Path
DATA = Path("../data/raw/goes16_mpshi"); OUT = Path("../results/v2"); OUT.mkdir(parents=True, exist_ok=True)
FR = 8   # fractional bits of the accumulator (without them the floor shift stops decay below 2^SH)
DAYS = ["20240101", "20240102", "20240103"]; Q = 1326.0; SH = 8


def load(day, tel):
    with Dataset(DATA / f"goes16_mpshi_{day}.nc") as nc:
        f = np.ma.filled(nc["IntgElectronFluxes"][:, tel].astype(float), np.nan)
        q = np.ma.filled(nc["IntgElectronFluxDQFs"][:, tel].astype(float), 255)
    valid = np.isfinite(f) & (q != 255)
    c = np.where(valid, np.maximum(np.round(np.nan_to_num(f) / Q), 0), 0).astype(np.int64)
    res = np.abs(np.nan_to_num(f)[valid & (f > 0)] / Q - np.round(np.nan_to_num(f)[valid & (f > 0)] / Q)).max()
    return c, valid, f, q, res


@nb.njit(cache=True)
def cdet(c, valid, sh, th, tl):
    n = len(c); a = np.zeros(n, np.bool_); L = np.int64(0); on = False; MAXL = (np.int64(1) << 32) - 1
    for i in range(n):
        if valid[i]:
            L = L - (L >> sh) + (c[i] << 8)
            if L > MAXL: L = MAXL
            if (not on) and L >= th: on = True
            elif on and L < tl: on = False
        a[i] = on
    return a


def thr(pfu): th = int(round(pfu / Q * (1 << (SH + FR)))); return th, int(round(th * 0.7))


def intervals(a, gap=60):
    idx = np.flatnonzero(a)
    if idx.size == 0: return 0
    return 1 + int(np.sum(np.diff(idx) > gap))


def profile(n, s, A, D, tr):
    g = np.zeros(n); e = s + tr + D + tr
    ramp = np.linspace(0, 1, tr + 1)[1:]
    g[s:s + tr] = ramp; g[s + tr:s + tr + D] = 1; g[s + tr + D:e] = ramp[::-1]
    return g


def main():
    # ---- data characterisation
    rows = []; cache = {}
    for d in DAYS:
        for tel in (0, 3):
            c, v, f, q, res = load(d, tel); cache[(d, tel)] = (c, v)
            neg = np.isfinite(f) & (f < 0); fill = ~v
            m5 = np.array([c[i:i + 300][v[i:i + 300]].mean() * Q if v[i:i + 300].sum() > 100 else np.nan for i in range(0, 86400, 300)])
            rows.append(dict(day=d, tel=tel + 1, invalid_pct=100 * fill.mean(), negative_as_zero_pct=100 * neg.mean(), mean_rate_cps=c[v].mean(),
                             max_5min_mean_pfu=np.nanmax(m5), bins_ge_1000pfu=int(np.nansum(m5 >= 1000)), quantum_residual=res))
    ch = pd.DataFrame(rows); ch.to_csv(OUT / "hazard_channel_characterisation.csv", index=False); print(ch.round(3).to_string(index=False))
    # ---- real-data operating curve
    rows = []
    for lvl in (400, 500, 600, 800, 1000, 2000):
        th, tl = thr(lvl); n = []
        for (d, tel), (c, v) in cache.items(): n.append(intervals(cdet(c, v, SH, th, tl)))
        rows.append(dict(alert_level_pfu=lvl, TH=th, TL=tl, alarms_over_6_channel_days=sum(n), per_channel_day_max=max(n)))
    oc = pd.DataFrame(rows); oc.to_csv(OUT / "hazard_real_operating_curve.csv", index=False); print("\nreal data alarms by alert level\n", oc.to_string(index=False))
    # ---- injection
    th, tl = thr(1000); BG = 350.0; rows = []; rng = np.random.default_rng(2024)
    for A in (0.3, 0.5, 0.8, 1.0, 2.0, 5.0, 10.0):
        for D in (300, 1800, 7200):
            hit = {1: 0, 2: 0}; tot = 0; delays = {1: [], 2: []}
            for (d, tel), (c, v) in [(k, cache[k]) for k in cache if k[1] == 0]:
                c4, v4 = cache[(d, 3)]
                for rep in range(8):
                    s = int(3600 + rep * 9000 + rng.integers(0, 1800)); tr = 60; e = s + 2 * tr + D
                    g = profile(86400, s, A, D, tr); lam = g * (A * 1000 / Q)
                    ci = c + rng.poisson(lam); cj = c4 + rng.poisson(lam)
                    a1 = cdet(ci.astype(np.int64), v, SH, th, tl); a4 = cdet(cj.astype(np.int64), v4, SH, th, tl)
                    if a1[s - 1] or a4[s - 1]: continue
                    tot += 1; va = a1 & a4
                    for k, a in ((1, a1), (2, va)):
                        w = a[s:e + 600]
                        if w.any(): hit[k] += 1; delays[k].append(int(np.argmax(w)))
            hz = (A * 1000 + BG) >= 1000
            rows.append(dict(A_x_1000pfu=A, total_mean_pfu=A * 1000 + BG, is_hazard=hz, dur_s=D, events=tot, P_alarm_1of1=hit[1] / tot, P_alarm_2of2=hit[2] / tot,
                             median_delay_s_1of1=np.median(delays[1]) if delays[1] else np.nan, median_delay_s_2of2=np.median(delays[2]) if delays[2] else np.nan))
    r = pd.DataFrame(rows); r.to_csv(OUT / "hazard_injection.csv", index=False); print("\ninjection (alert 1000 pfu, window 2^8 s)\n", r.round(3).to_string(index=False))


if __name__ == "__main__": main()
