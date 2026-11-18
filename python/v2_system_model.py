"""
SpaceGuard v2 - system-level golden model: 2 detector channels -> k-of-N vote -> mitigation controller.
Writes RTL vectors: <name>_c{0,1}_{x,v,alarm}.hex and <name>_mode_k{1,2}.hex   (python v2_system_model.py <outdir>)
Controller sees the detector alarms registered after the previous sample (1-sample latency) and valid of the current sample.
"""
import sys
from pathlib import Path
import numpy as np
import v2_rtl_refmodel as M

HOLD, ESC, WD = 600, 300, 60


def ctrl(alarms, valids, K, HOLD=HOLD, ESC=ESC, WD=WD):
    n = len(alarms[0]); mode = hold = esc = wd = 0; out = []
    for i in range(n):
        ev_in = sum(a[i - 1] if i else 0 for a in alarms) >= K
        wd = 0 if any(v[i] for v in valids) else min(wd + 1, 255)
        wdf = wd >= WD; ev = ev_in or wdf
        if mode == 0:
            if ev: mode, hold, esc = 1, HOLD, 0
        elif mode == 1:
            if ev:
                hold = HOLD; esc = min(esc + 1, 1023)
                if esc >= ESC: mode = 2
            elif hold <= 1: mode, hold, esc = 0, 0, 0
            else: hold -= 1
        else:
            if ev: hold = HOLD
            elif hold <= 1: mode, hold, esc = 1, HOLD, 0
            else: hold -= 1
        if wdf: mode, hold = 2, HOLD
        out.append(mode)
    return out


def write(out, name, chans):
    als, vs = [], []
    for c, (xf, v) in enumerate(chans):
        xi, al, _ = M.write_vectors(out, f"{name}_c{c}", xf, v)
        als.append(al); vs.append([bool(t) for t in v])
    for K in (1, 2):
        mo = ctrl(als, vs, K)
        (out / f"{name}_mode_k{K}.hex").write_text("\n".join(str(m) for m in mo) + "\n")
        print(f"{name} K={K}: samples in NORMAL/ALERT/PROTECT = {[mo.count(i) for i in range(3)]}")


if __name__ == "__main__":
    import v2_injection_benchmark as B
    out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
    # coherent event on both bands (same injected events), clean otherwise
    ch = []
    for b in (0, 1):
        x, v = B.load("20240103", b)
        evs = B.make_events(np.random.default_rng(123), v)
        ch.append((B.inject(x, evs), v))
    write(out, "sys_inj", ch)
    # dropouts: single-channel bursts plus one common blackout (triggers the watchdog fail-safe)
    rng = np.random.default_rng(11); ch = []
    for b in (0, 1):
        x, v = B.load("20240102", b)
        evs = B.make_events(np.random.default_rng(55), v)
        ch.append([B.inject(x, evs), v.copy()])
    for _, v in ch:
        for s in rng.integers(0, 86000, 30): v[s:s + int(rng.integers(1, 100))] = False
    for _, v in ch: v[40000:40090] = False; v[70000:70070] = False
    write(out, "sys_drop", [(xx, vv) for xx, vv in ch])
