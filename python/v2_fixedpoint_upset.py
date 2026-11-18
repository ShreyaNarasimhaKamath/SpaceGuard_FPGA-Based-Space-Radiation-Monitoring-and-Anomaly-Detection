"""
SpaceGuard v2 - (A) fixed-point quantisation and (B,C) single-event-upset study of the detector's own state.
Detector: asymmetric EMA, alpha_up=2^-10, alpha_down=2^-7, alarm if x > 1.5*baseline (k=0.5, shift-add),
persistence M=10, release R=30, optional freeze-on-alarm. 32-bit unsigned baseline register, F fractional bits.
Upset = one bit flipped in the baseline register (EMA) or in one stored sample of the 1800-sample median buffer.
Errors are counted as alarm-state mismatches vs the un-upset run over the following 3 h.
"""
import numpy as np, pandas as pd, numba as nb
import v2_injection_benchmark as B

OUT = B.OUT
NU, ND, KMODE = 10, 7, 0       # k = 0.5  -> thr = s + (s>>1)
WIN = 10800


@nb.njit(cache=True)
def ema_fx(xi, valid, nu, nd, kmode, M, R, freeze, tmax, warm, up_idx, up_bit):
    n = len(xi); a = np.zeros(n, np.bool_); state = np.int64(0); init = False
    run = 0; off = 0; on = False; frz = 0; cnt = 0; MAXS = (np.int64(1) << 32) - 1
    for i in range(n):
        if i == up_idx and init:
            state ^= (np.int64(1) << up_bit)
        if not valid[i]:
            a[i] = on; continue
        x = xi[i]
        if not init:
            state = x; init = True; cnt = 1; continue
        if kmode == 0: thr = state + (state >> 1)
        elif kmode == 1: thr = state + (state >> 1) + (state >> 2)
        else: thr = state << 1
        cond = (cnt >= warm) and (x > thr)
        if cond: run += 1; off = 0
        else: run = 0; off += 1
        if (not on) and run >= M: on = True; frz = 0
        elif on and off >= R: on = False
        a[i] = on
        if freeze and on and frz < tmax: frz += 1
        else:
            d = x - state
            if d > 0: state += d >> nu
            else: state += d >> nd
            if state < 0: state = 0
            if state > MAXS: state = MAXS
        cnt += 1
    return a


def toint(x, valid, F):
    return np.where(valid, np.round(np.nan_to_num(x) * (1 << F)), 0).astype(np.int64)


def series():
    for split, dates in (("dev", B.DEV), ("test", B.TEST)):
        for d in dates:
            for b in B.BANDS:
                x, v = B.load(d, b); yield split, d, b + 1, x, v


def part_a():
    rows = []
    for F in (0, 2, 4, 8, 12):
        agree, nat, pd_hits, pd_n = [], [], 0, 0
        for split, d, b, x, v in series():
            fa = B.ema_alarm(np.nan_to_num(x), v, 2.0 ** -NU, 2.0 ** -ND, 0.5, B.M, B.R, True, B.TMAX, B.WARM)
            qa = ema_fx(toint(x, v, F), v, NU, ND, KMODE, B.M, B.R, True, B.TMAX, B.WARM, -1, 0)
            agree.append((fa[B.WARM:] == qa[B.WARM:]).mean()); nat.append(B.count_intervals(qa))
            for rep in range(10):
                rng = np.random.default_rng(5000 + rep + 7 * b)
                evs = B.make_events(rng, v); xi = B.inject(x, evs)
                qi = ema_fx(toint(xi, v, F), v, NU, ND, KMODE, B.M, B.R, True, B.TMAX, B.WARM, -1, 0)
                for ev in evs:
                    if ev["A"] >= 1 and ev["D"] >= 300 and not qi[ev["s"] - 1]:
                        pd_n += 1; pd_hits += int(qi[ev["s"]:ev["e"] + 1].any())
        rows.append(dict(frac_bits=F, alarm_agreement_vs_float=np.mean(agree), natural_alarms_per_day=np.mean(nat), Pd_core=pd_hits / pd_n))
    r = pd.DataFrame(rows); r.to_csv(OUT / "fixedpoint_quantisation.csv", index=False); print("\n== A. Fixed-point quantisation (32-bit baseline register) ==\n", r.round(4).to_string(index=False))


def part_b(freeze):
    F = 12; rng = np.random.default_rng(1); recs = []
    for split, d, b, x, v in series():
        xi = toint(x, v, F)
        a0 = ema_fx(xi, v, NU, ND, KMODE, B.M, B.R, freeze, B.TMAX, B.WARM, -1, 0)
        for bit in range(32):
            for t in rng.integers(4000, 80000, 30):
                a = ema_fx(xi, v, NU, ND, KMODE, B.M, B.R, freeze, B.TMAX, B.WARM, int(t), bit)
                w0, w1 = a0[t:t + WIN], a[t:t + WIN]
                recs.append(dict(bit=bit, fa=int((w1 & ~w0).sum()), miss=int((~w1 & w0).sum())))
    df = pd.DataFrame(recs); df["err"] = df.fa + df.miss
    g = df.groupby("bit").agg(frac_with_error=("err", lambda s: (s > 0).mean()), median_err_s=("err", "median"),
                              p90_err_s=("err", lambda s: s.quantile(.9)), mean_false_alarm_s=("fa", "mean"), mean_missed_s=("miss", "mean")).reset_index()
    g.insert(1, "flux_weight_2^j", g.bit - F)
    g.to_csv(OUT / f"upset_ema_{'freeze' if freeze else 'nofreeze'}.csv", index=False)
    print(f"\n== B. EMA baseline-register upsets (freeze={freeze}); error seconds within 3 h of upset ==\n",
          g[g.bit.isin([0, 4, 8, 10, 12, 14, 16, 18, 20, 22, 24, 28, 31])].round(2).to_string(index=False))


def part_c():
    F = 12; rng = np.random.default_rng(2); recs = []
    bits = [0, 4, 8, 12, 14, 16, 18, 20, 24, 28, 31]
    for split, d, b, x, v in series():
        xi = toint(x, v, F)
        for bit in bits:
            for _ in range(8):
                t = int(rng.integers(6000, 78000)); j = t - int(rng.integers(0, 1800))
                lo, hi = t - 3600, t + 3600
                xs = x[lo:hi].copy(); vs = v[lo:hi]
                bad = xs.copy(); bad[j - lo] = float(int(xi[j]) ^ (1 << bit)) / (1 << F)
                def alarm(base_series):
                    bm = B.med_baseline(base_series); cond = np.isfinite(bm) & (xs > 1.5 * bm) & vs; cond[:1800] = False
                    return B.alarm_from_cond(cond, vs, B.M, B.R)
                a0, a1 = alarm(xs), alarm(bad)
                w0, w1 = a0[t - lo:], a1[t - lo:]
                recs.append(dict(bit=bit, fa=int((w1 & ~w0).sum()), miss=int((~w1 & w0).sum())))
    df = pd.DataFrame(recs); df["err"] = df.fa + df.miss
    g = df.groupby("bit").agg(frac_with_error=("err", lambda s: (s > 0).mean()), median_err_s=("err", "median"),
                              p90_err_s=("err", lambda s: s.quantile(.9)), max_err_s=("err", "max")).reset_index()
    g.to_csv(OUT / "upset_median_buffer.csv", index=False); print("\n== C. Median detector: one stored sample upset (24-31 bit word) ==\n", g.round(2).to_string(index=False))


if __name__ == "__main__":
    part_a(); part_b(True); part_b(False); part_c()
