"""
SpaceGuard: hysteresis event-detection experiment (roadmap Phase 1, 20.1/20.2/20.3).

Asymmetric EMA Configuration A is held fixed (alpha_up=0.001, alpha_down=0.01).
Only the event-detection logic changes:
  H0 : single threshold (start = end = 1.0 x T)  -> must reproduce the existing Config A results
  H1-H4 : start at 1.0 x T, continue while residual > f x T, f = 0.8, 0.7, 0.6, 0.5

Everything else (calibration, threshold T, merge gap, min-active, common mask,
rolling-median reference) is reused unchanged from asymmetric_ema_sensitivity.py
so the numbers are directly comparable.

Matching is reported two ways:
  "filtered"  : existing method - Hungarian on IoU, then drop matches below the IoU threshold
  "rematched" : roadmap 20.2 - eligibility (IoU >= threshold) is applied BEFORE assignment,
                maximising the number of valid matches first, total IoU second.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment

import asymmetric_ema_sensitivity as base

ROOT = base.ROOT
RESULTS_DIR = base.RESULTS_DIR
PLOTS_DIR = base.PLOTS_DIR

ALPHA_UP, ALPHA_DOWN = base.CONFIGURATIONS["A"]
SETTINGS = {"H0": 1.0, "H1": 0.8, "H2": 0.7, "H3": 0.6, "H4": 0.5,
            "H1b": 0.8, "H2b": 0.7, "H3b": 0.6, "H4b": 0.5}
# "b" variants: minimum-duration filter counts only samples above the HIGH threshold
# (plain variants count every active sample, including low-threshold continuation).
IOU_THRESHOLDS = [0.30, 0.50]


def hysteresis_mask(residual, common_mask, t_high, t_low):
    """Sample is active from the first residual > t_high until residual <= t_low.
    A sample outside the common mask (invalid / not finite) is never active and
    ends any active run; the existing merge-gap logic can re-join the pieces."""
    n = len(residual)
    active = np.zeros(n, dtype=bool)
    on = False
    for i in range(n):
        if not common_mask[i]:
            on = False
            continue
        r = residual[i]
        if on:
            on = r > t_low
        else:
            on = r > t_high
        active[i] = on
    return active


def extract_events_high_count(active, high, min_high):
    """Group `active` with the base merge gap; keep groups with >= min_high samples above HIGH."""
    idx = np.flatnonzero(active)
    if len(idx) == 0:
        return []
    out = []
    for g in np.split(idx, np.where(np.diff(idx) > base.MERGE_GAP)[0] + 1):
        if high[g].sum() < min_high:
            continue
        out.append({"start": int(g[0]), "end": int(g[-1]),
                    "span_samples": int(g[-1] - g[0] + 1), "active_samples": len(g)})
    return out


def iou_matrix(ref, cand):
    m = np.zeros((len(ref), len(cand)))
    for i, r in enumerate(ref):
        for j, c in enumerate(cand):
            m[i, j] = base.calculate_event_iou(r, c)
    return m


def match(ref, cand, scores, thr, mode):
    """mode 'filtered' or 'rematched'. Returns list of (i, j, iou)."""
    if not ref or not cand:
        return []
    if mode == "filtered":
        rows, cols = linear_sum_assignment(-scores)
        return [(i, j, scores[i, j]) for i, j in zip(rows, cols)
                if scores[i, j] > 0 and scores[i, j] >= thr]
    eligible = scores >= thr
    cost = -(eligible * 1000.0 + scores * eligible)
    rows, cols = linear_sum_assignment(cost)
    return [(i, j, scores[i, j]) for i, j in zip(rows, cols) if eligible[i, j]]


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    date_rows, match_rows = [], []

    for date in base.DATES:
        for band_index in base.BANDS:
            band = band_index + 1
            signal, valid = base.load_channel(date, band_index)
            ref_baseline = base.calculate_reference(signal)
            ref_residual = signal - ref_baseline
            _, residual = base.calculate_asymmetric_ema(signal, valid, ALPHA_UP, ALPHA_DOWN)

            common = valid & np.isfinite(ref_residual) & np.isfinite(residual)
            calib = common & (np.arange(len(signal)) < base.CALIBRATION_SAMPLES)

            _, _, _, t_ref = base.calculate_threshold(ref_residual, calib)
            ref_events = base.extract_events(common & (ref_residual > t_ref))
            _, _, sigma, T = base.calculate_threshold(residual, calib)

            print(f"{date} band {band}: T={T:.1f} ref_events={len(ref_events)}")

            for name, frac in SETTINGS.items():
                mask = hysteresis_mask(residual, common, T, frac * T)
                if name.endswith("b"):
                    events = extract_events_high_count(mask, common & (residual > T), base.MIN_ACTIVE_SAMPLES)
                else:
                    events = base.extract_events(mask)
                scores = iou_matrix(ref_events, events)

                # merge / fragmentation diagnostics (positive overlap)
                if len(ref_events) and len(events):
                    ov = scores > 0
                    merged = int((ov.sum(axis=0) > 1).sum())       # candidates spanning >1 reference event
                    fragmented = int((ov.sum(axis=1) > 1).sum())   # reference events split across >1 candidates
                    unmatched_cand = int((ov.sum(axis=0) == 0).sum())
                else:
                    merged = fragmented = 0
                    unmatched_cand = len(events)

                date_rows.append({
                    "date": date, "band": band, "setting": name, "low_fraction": frac,
                    "threshold_high": T, "threshold_low": frac * T,
                    "flagged_percentage": 100 * mask.sum() / common.sum(),
                    "reference_events": len(ref_events), "candidate_events": len(events),
                    "candidates_spanning_multiple_refs": merged,
                    "refs_split_across_candidates": fragmented,
                    "candidates_with_no_overlap": unmatched_cand,
                })

                for mode in ("filtered", "rematched"):
                    for thr in IOU_THRESHOLDS:
                        for i, j, s in match(ref_events, events, scores, thr, mode):
                            match_rows.append({
                                "date": date, "band": band, "setting": name,
                                "low_fraction": frac, "mode": mode, "iou_threshold": thr,
                                "reference_event": i + 1, "candidate_event": j + 1, "iou": s,
                                "onset_error_s": events[j]["start"] - ref_events[i]["start"],
                                "duration_error_s": events[j]["span_samples"] - ref_events[i]["span_samples"],
                            })

    per_date = pd.DataFrame(date_rows)
    matches = pd.DataFrame(match_rows)

    quality = []
    for (setting, band), g in per_date.groupby(["setting", "band"]):
        R, C = int(g.reference_events.sum()), int(g.candidate_events.sum())
        for mode in ("filtered", "rematched"):
            for thr in IOU_THRESHOLDS:
                m = matches[(matches.setting == setting) & (matches.band == band)
                            & (matches["mode"] == mode) & (matches.iou_threshold == thr)]
                tp = len(m)
                p = tp / C if C else 0.0
                r = tp / R if R else 0.0
                f1 = 2 * p * r / (p + r) if p + r else 0.0
                quality.append({
                    "setting": setting, "low_fraction": SETTINGS[setting], "band": band,
                    "mode": mode, "iou_threshold": thr,
                    "reference_events": R, "candidate_events": C, "matched": tp,
                    "precision": p, "recall": r, "f1": f1,
                    "mean_iou": m.iou.mean() if tp else np.nan,
                    "mean_onset_error_s": m.onset_error_s.mean() if tp else np.nan,
                    "median_abs_onset_error_s": m.onset_error_s.abs().median() if tp else np.nan,
                    "mean_duration_error_s": m.duration_error_s.mean() if tp else np.nan,
                    "mean_abs_duration_error_s": m.duration_error_s.abs().mean() if tp else np.nan,
                    "mean_flagged_pct": g.flagged_percentage.mean(),
                    "merged_candidates": int(g.candidates_spanning_multiple_refs.sum()),
                    "split_refs": int(g.refs_split_across_candidates.sum()),
                })
    quality = pd.DataFrame(quality)

    per_date.to_csv(RESULTS_DIR / "hysteresis_by_date.csv", index=False)
    matches.to_csv(RESULTS_DIR / "hysteresis_event_matches.csv", index=False)
    quality.to_csv(RESULTS_DIR / "hysteresis_event_quality.csv", index=False)

    pd.set_option("display.width", 220)
    for thr in IOU_THRESHOLDS:
        for mode in ("filtered", "rematched"):
            print(f"\n=== IoU>={thr:.2f}  matching={mode} ===")
            print(quality[(quality.iou_threshold == thr) & (quality["mode"] == mode)][
                ["setting", "low_fraction", "band", "candidate_events", "matched", "precision",
                 "recall", "f1", "mean_iou", "median_abs_onset_error_s",
                 "mean_duration_error_s", "mean_flagged_pct", "merged_candidates", "split_refs"]
            ].round(3).to_string(index=False))

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
    q = quality[(quality.iou_threshold == 0.30) & (quality["mode"] == "rematched")]
    for band in (1, 2):
        s = q[q.band == band].sort_values("low_fraction")
        axes[0].plot(s.low_fraction, s.f1, marker="o", label=f"Band {band}")
        axes[1].plot(s.low_fraction, s.mean_duration_error_s, marker="o", label=f"Band {band}")
        axes[2].plot(s.low_fraction, s.mean_flagged_pct, marker="o", label=f"Band {band}")
    for ax, t in zip(axes, ["F1 (IoU>=0.30, rematched)", "Mean duration error (s)", "Mean flagged samples (%)"]):
        ax.set_title(t); ax.set_xlabel("Low threshold / high threshold"); ax.grid(alpha=0.3); ax.legend()
    fig.suptitle("Hysteresis experiment - asymmetric EMA Config A fixed")
    fig.tight_layout()
    fig.savefig(PLOTS_DIR / "hysteresis_experiment.png", dpi=150)
    print("\nSaved: hysteresis_by_date.csv, hysteresis_event_matches.csv, hysteresis_event_quality.csv, hysteresis_experiment.png")


if __name__ == "__main__":
    main()
