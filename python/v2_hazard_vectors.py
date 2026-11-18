"""RTL vectors for the hazard (count) system: tel1 + tel4 >2 MeV integral -> count detector -> vote -> controller.
python v2_hazard_vectors.py <outdir>   writes hz_<set>_c{0,1}_{x,v,alarm}.hex and hz_<set>_mode_k{1,2}.hex"""
import sys
from pathlib import Path
import numpy as np
import v2_hazard_counts as H
import v2_system_model as S


def write(out, name, chans):
    als, vs = [], []
    th, tl = H.thr(1000)
    for c, (cnt, v) in enumerate(chans):
        a = H.cdet(cnt.astype(np.int64), v, H.SH, th, tl).astype(int)
        (out / f"{name}_c{c}_x.hex").write_text("\n".join(f"{int(t):08x}" for t in cnt) + "\n")
        (out / f"{name}_c{c}_v.hex").write_text("\n".join("1" if t else "0" for t in v) + "\n")
        (out / f"{name}_c{c}_alarm.hex").write_text("\n".join(str(t) for t in a) + "\n")
        als.append(list(a)); vs.append([bool(t) for t in v])
    for K in (1, 2):
        mo = S.ctrl(als, vs, K)
        (out / f"{name}_mode_k{K}.hex").write_text("\n".join(str(m) for m in mo) + "\n")
        print(f"{name} K={K}: alarm samples c0={sum(als[0])} c1={sum(als[1])}  NORMAL/ALERT/PROTECT={[mo.count(i) for i in range(3)]}")


if __name__ == "__main__":
    out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
    th, tl = H.thr(1000); rng = np.random.default_rng(99)
    c0, v0, *_ = H.load("20240101", 0); c1, v1, *_ = H.load("20240101", 3)
    write(out, "hz_real", [(c0, v0), (c1, v1)])                       # clean real day
    ch = []
    for tel in (0, 3):
        c, v, *_ = H.load("20240102", tel); ch.append([c.copy(), v.copy()])
    for k in range(6):                                                  # coherent injected events (A=2 x 1000 pfu, mixed durations)
        s = 4000 + k * 12000; D = (600, 1800, 3600)[k % 3]; g = H.profile(86400, s, 2.0, D, 60); lam = g * (2.0 * 1000 / H.Q)
        for ch_ in ch: ch_[0] += rng.poisson(lam)
    for _, v in ch:                                                     # dropout bursts + common blackout (watchdog)
        for s in rng.integers(0, 86000, 20): v[s:s + int(rng.integers(1, 80))] = False
        v[60000:60090] = False
    write(out, "hz_inj", [(a, b) for a, b in ch])
