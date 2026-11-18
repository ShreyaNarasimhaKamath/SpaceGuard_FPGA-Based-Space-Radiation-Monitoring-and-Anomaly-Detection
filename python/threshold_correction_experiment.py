"""
Experiment 1: calibration-median correction of the detection threshold.
Asymmetric EMA Config A fixed (0.001 / 0.01). Event grouping unchanged.

  centred        : T = 6 * 1.4826 * MAD                      (current behaviour)
  median_offset  : T = median(r_cal) + 6 * 1.4826 * MAD

Two reference treatments:
  ref_fixed : rolling-median reference keeps its original (centred) threshold
  ref_same  : rolling-median reference uses the same correction as the candidate
"""
import numpy as np
import pandas as pd

import asymmetric_ema_sensitivity as base
import hysteresis_event_detection as hyst

RES = base.RESULTS_DIR
a_up, a_down = base.CONFIGURATIONS["A"]
THRS = [0.30, 0.50]


def threshold(residual, calib, offset):
    med, mad, sigma, t = base.calculate_threshold(residual, calib)
    return (med + t) if offset else t, med, sigma


rows, mrows = [], []
for date in base.DATES:
    for bi in base.BANDS:
        band = bi + 1
        sig, valid = base.load_channel(date, bi)
        ref_res = sig - base.calculate_reference(sig)
        _, res = base.calculate_asymmetric_ema(sig, valid, a_up, a_down)
        common = valid & np.isfinite(ref_res) & np.isfinite(res)
        calib = common & (np.arange(len(sig)) < base.CALIBRATION_SAMPLES)

        t_ref_c, ref_med, _ = threshold(ref_res, calib, False)
        t_ref_o, _, _ = threshold(ref_res, calib, True)
        ev_ref = {
            "ref_fixed": base.extract_events(common & (ref_res > t_ref_c)),
            "ref_same_offset": base.extract_events(common & (ref_res > t_ref_o)),
        }
        for mode, off in (("centred", False), ("median_offset", True)):
            T, med, sigma = threshold(res, calib, off)
            mask = common & (res > T)
            ev = base.extract_events(mask)
            for refname, evr in ev_ref.items():
                if refname == "ref_same_offset" and not off:
                    evr = ev_ref["ref_fixed"]          # centred vs centred
                scores = hyst.iou_matrix(evr, ev)
                rows.append(dict(date=date, band=band, threshold_mode=mode, reference=refname,
                                 cal_median=med, ref_cal_median=ref_med, sigma=sigma, threshold=T,
                                 effective_sigma_above_median=(T - med) / sigma,
                                 flagged_pct=100 * mask.sum() / common.sum(),
                                 ref_events=len(evr), cand_events=len(ev)))
                for thr in THRS:
                    for i, j, s in hyst.match(evr, ev, scores, thr, "rematched"):
                        mrows.append(dict(date=date, band=band, threshold_mode=mode, reference=refname,
                                          iou_threshold=thr, iou=s,
                                          onset_error_s=ev[j]["start"] - evr[i]["start"],
                                          duration_error_s=ev[j]["span_samples"] - evr[i]["span_samples"]))
        print("done", date, band, flush=True)

d = pd.DataFrame(rows); m = pd.DataFrame(mrows)
d.to_csv(RES / "threshold_correction_by_date.csv", index=False)
m.to_csv(RES / "threshold_correction_matches.csv", index=False)

q = []
for (band, mode, ref), g in d.groupby(["band", "threshold_mode", "reference"]):
    R, C = int(g.ref_events.sum()), int(g.cand_events.sum())
    for thr in THRS:
        x = m[(m.band == band) & (m.threshold_mode == mode) & (m.reference == ref) & (m.iou_threshold == thr)]
        tp = len(x); p = tp / C if C else 0; r = tp / R if R else 0
        q.append(dict(band=band, threshold_mode=mode, reference=ref, iou_threshold=thr, ref_events=R,
                      cand_events=C, matched=tp, precision=p, recall=r,
                      f1=2 * p * r / (p + r) if p + r else 0, mean_iou=x.iou.mean(),
                      med_abs_onset_s=x.onset_error_s.abs().median(),
                      mean_duration_err_s=x.duration_error_s.mean(),
                      mean_flagged_pct=g.flagged_pct.mean()))
q = pd.DataFrame(q); q.to_csv(RES / "threshold_correction_quality.csv", index=False)
pd.set_option("display.width", 220)
print(q.round(3).to_string(index=False))
print()
print(d[d.reference == "ref_fixed"][["date", "band", "threshold_mode", "cal_median", "threshold",
      "effective_sigma_above_median", "flagged_pct", "ref_events", "cand_events"]].round(2).to_string(index=False))
print("\nreference calibration medians:\n", d[d.threshold_mode == "centred"][["date", "band", "ref_cal_median"]].drop_duplicates().round(1).to_string(index=False))
