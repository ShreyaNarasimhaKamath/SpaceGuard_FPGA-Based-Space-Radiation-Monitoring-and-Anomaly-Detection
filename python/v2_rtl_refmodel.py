"""
SpaceGuard v2 - bit-exact integer reference model of the single-channel detector (golden model for the RTL).

Parameters (board independent): 32-bit unsigned baseline register, input x is an unsigned 32-bit fixed-point
value (F=8 fractional bits). EMA: base += (x-base)>>NU if x>base else >>ND (arithmetic shift, floor).
Alarm condition: x > base + (base>>1) (k=0.5), armed after WARM valid samples.
Persistence M, release R, optional baseline freeze while alarm (at most TMAX samples per alarm).
Invalid samples hold all state. Alarm output = 'on' flag after the sample is processed.

Also generates RTL test vectors (hex files) from the GOES data:  python v2_rtl_refmodel.py <outdir>
"""
import sys
from pathlib import Path
import numpy as np

NU, ND, M, R, WARM, TMAX, FREEZE, F = 10, 7, 10, 30, 1800, 1800, 1, 8
MAXS = (1 << 32) - 1


def ref(xs, valid):
    """xs: list of ints in [0, 2^32-1]; returns (alarm list, baseline-after-sample list)."""
    base = 0; init = False; cnt = 0; run = 0; off = 0; on = 0; frz = 0
    al, bs = [], []
    for x, v in zip(xs, valid):
        if v:
            if not init:
                base, init, cnt = x, True, 1
            else:
                thr = base + (base >> 1)
                cond = cnt >= WARM and x > thr
                if cond: run = min(run + 1, M); off = 0
                else: run = 0; off = min(off + 1, R)
                on_set = (not on) and run >= M
                on_clr = bool(on) and off >= R
                if on_set: on = 1; frz = 0
                elif on_clr: on = 0
                hold = FREEZE and on and frz < TMAX
                if hold: frz += 1
                else:
                    d = x - base
                    base += (d >> NU) if d > 0 else (d >> ND)
                    base = min(max(base, 0), MAXS)
                cnt = min(cnt + 1, WARM)
        al.append(int(on)); bs.append(base)
    return al, bs


def write_vectors(out, name, xf, valid):
    xi = [int(round(v * (1 << F))) if ok else 0 for v, ok in zip(np.nan_to_num(xf), valid)]
    assert max(xi) <= MAXS
    al, bs = ref(xi, [bool(v) for v in valid])
    (out / f"{name}_x.hex").write_text("\n".join(f"{v:08x}" for v in xi) + "\n")
    (out / f"{name}_v.hex").write_text("\n".join("1" if v else "0" for v in valid) + "\n")
    (out / f"{name}_alarm.hex").write_text("\n".join(str(a) for a in al) + "\n")
    (out / f"{name}_base.hex").write_text("\n".join(f"{b:08x}" for b in bs) + "\n")
    print(f"{name}: n={len(xi)} valid={int(np.sum(valid))} alarm_samples={sum(al)} max_x={max(xi)}")
    return xi, al, bs


if __name__ == "__main__":
    import v2_injection_benchmark as B
    import v2_fixedpoint_upset as U
    out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
    sets = []
    x, v = B.load("20240103", 0)
    evs = B.make_events(np.random.default_rng(123), v)
    sets.append(("inj_d3_b1", B.inject(x, evs), v))
    x, v = B.load("20240101", 1); sets.append(("real_d1_b2", x, v))
    for name, xf, v in sets:
        xi, al, bs = write_vectors(out, name, xf, v)
        # cross-check against the numba model used for the study
        qa = U.ema_fx(np.array(xi, np.int64), v, NU, ND, 0, M, R, True, TMAX, WARM, -1, 0)
        print("   matches numba ema_fx:", bool(np.array_equal(qa.astype(int), np.array(al))))
