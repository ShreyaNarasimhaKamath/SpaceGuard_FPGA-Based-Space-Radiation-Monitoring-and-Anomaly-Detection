"""
SpaceGuard v2 - controller timing study (ESC, HOLD) on the hazard-channel system (>2 MeV integral, telescopes 1+4, 2-of-2 vote).
Per day: 6 synthetic Poisson events injected coherently on both telescopes at (A x 1000 pfu, duration s):
  (2,600) (2,1800) (2,3600) (5,600) (1,1800) hazard-level;  (0.5,1800) sub-hazard (must NOT trigger).
Metrics per (ESC, HOLD): events reaching ALERT / PROTECT, median latency to each, hangover = samples non-NORMAL after the
event ends, mode transitions per event (flapping), non-NORMAL samples on the clean (un-injected) days, sub-hazard triggers.
"""
import numpy as np, pandas as pd
import v2_hazard_counts as H, v2_system_model as S
OUT = H.OUT
EV = [(2, 600), (2, 1800), (2, 3600), (5, 600), (1, 1800), (0.5, 1800)]; TR = 60


def build(day, rng, inject=True):
    ch = []; spans = []
    for tel in (0, 3):
        c, v, *_ = H.load(day, tel); ch.append([c.copy(), v])
    if inject:
        for k, (A, D) in enumerate(EV):
            s = 4000 + k * 12000; g = H.profile(86400, s, A, D, TR); lam = g * (A * 1000 / H.Q)
            for x in ch: x[0] += rng.poisson(lam)
            spans.append((A, D, s, s + 2 * TR + D))
    th, tl = H.thr(1000)
    al = [list(H.cdet(c.astype(np.int64), v, H.SH, th, tl).astype(int)) for c, v in ch]
    return al, [list(map(bool, v)) for _, v in ch], spans


def main():
    rng = np.random.default_rng(321); days = ["20240101", "20240102", "20240103"]
    inj = [build(d, rng) for d in days]; clean = [build(d, rng, False) for d in days]
    rows = []
    for ESC in (60, 120, 300, 600):
        for HOLD in (120, 300, 600, 1800):
            m = dict(ESC=ESC, HOLD=HOLD); n_h = n_a = n_p = 0; la = []; lp = []; hang = []; flap = []; sub = 0; n_sub = 0
            for al, vs, spans in inj:
                mo = np.array(S.ctrl(al, vs, 2, HOLD, ESC, 60))
                for A, D, s, e in spans:
                    w = mo[s:e + 1]
                    if A >= 0.8:
                        n_h += 1
                        if (w > 0).any(): n_a += 1; la.append(int(np.argmax(w > 0)))
                        if (w == 2).any(): n_p += 1; lp.append(int(np.argmax(w == 2)))
                        tail = mo[e:e + 20000]; z = np.flatnonzero(tail == 0); hang.append(int(z[0]) if z.size else len(tail))
                        seg = mo[s - 60:e + (hang[-1] if hang else 0) + 1]; flap.append(int(np.sum(np.diff(seg) != 0)))
                    else:
                        n_sub += 1; sub += int((mo[s:e + 1 + 3000] > 0).any())
            nn = sum(int((np.array(S.ctrl(al, vs, 2, HOLD, ESC, 60)) > 0).sum()) for al, vs, _ in clean)
            rows.append(dict(ESC=ESC, HOLD=HOLD, hazard_events=n_h, reach_ALERT=n_a, reach_PROTECT=n_p,
                             median_lat_ALERT_s=np.median(la) if la else np.nan, median_lat_PROTECT_s=np.median(lp) if lp else np.nan,
                             median_hangover_s=np.median(hang), mode_changes_per_event=np.mean(flap),
                             subhazard_triggers=f"{sub}/{n_sub}", clean_nonNORMAL_samples=nn))
    r = pd.DataFrame(rows); r.to_csv(OUT / "controller_timing_study.csv", index=False); print(r.round(1).to_string(index=False))


if __name__ == "__main__": main()
