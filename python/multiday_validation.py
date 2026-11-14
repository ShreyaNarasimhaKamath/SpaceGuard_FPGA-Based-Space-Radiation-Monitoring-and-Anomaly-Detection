
import numpy as np
import pandas as pd
from netCDF4 import Dataset
from pathlib import Path

FILES = sorted(
    Path("data/raw/goes16_mpshi").glob("goes16_mpshi_2024010[1-3].nc")
)

BASELINE_WINDOW = 1800
CALIBRATION_HOURS = 10
THRESHOLD_MULTIPLIER = 6
MAX_GAP_SECONDS = 120
MIN_EVENT_SECONDS = 60

TELESCOPE = 0
BANDS = [0, 1]

OUTPUT = Path("results/evaluation")
OUTPUT.mkdir(parents=True, exist_ok=True)


def detect_events(candidate, flux, times):

    indices = np.where(candidate)[0]
    events = []

    if len(indices) == 0:
        return events

    start = previous = indices[0]
    active_count = 1

    groups = []

    for idx in indices[1:]:

        gap = idx - previous - 1

        if gap <= MAX_GAP_SECONDS:
            previous = idx
            active_count += 1

        else:
            groups.append((start, previous, active_count))
            start = previous = idx
            active_count = 1

    groups.append((start, previous, active_count))

    for start, end, active_count in groups:

        if active_count < MIN_EVENT_SECONDS:
            continue

        events.append({
            "start_hour": times[start],
            "end_hour": times[end],
            "span_seconds": end - start + 1,
            "active_samples": active_count,
            "peak_flux": float(np.nanmax(flux[start:end + 1]))
        })

    return events


all_summaries = []
all_events = []

for file in FILES:

    date = file.stem[-8:]

    print(f"\n========== {date} ==========")

    with Dataset(file, "r") as ds:

        timestamps = np.asarray(
            ds.variables["L1a_SciData_TimeStamp"][:],
            dtype=np.float64
        )

        time_hours = (
            timestamps - timestamps[0]
        ) / 3600

        for band in BANDS:

            flux = np.asarray(
                ds.variables["DiffElectronFluxes"][
                    :, TELESCOPE, band
                ],
                dtype=np.float64
            )

            dqf = np.asarray(
                ds.variables["DiffElectronFluxDQFs"][
                    :, TELESCOPE, band
                ]
            )

            valid = (
                (dqf == 0)
                & np.isfinite(flux)
                & (flux >= 0)
            )

            signal = pd.Series(
                np.where(valid, flux, np.nan)
            )

            baseline = signal.rolling(
                window=BASELINE_WINDOW,
                min_periods=BASELINE_WINDOW
            ).median()

            residual = signal - baseline

            calibration = residual[
                time_hours < CALIBRATION_HOURS
            ].dropna()

            if len(calibration) == 0:
                raise RuntimeError(
                    f"No calibration data: {date}, band {band + 1}"
                )

            calibration_values = calibration.to_numpy()

            residual_median = np.median(calibration_values)

            mad = np.median(
                np.abs(calibration_values - residual_median)
            )

            sigma = 1.4826 * mad
            threshold = THRESHOLD_MULTIPLIER * sigma

            if not np.isfinite(threshold) or threshold <= 0:
                raise RuntimeError(
                    f"Invalid threshold: {date}, band {band + 1}"
                )

            candidate = (
                (residual > threshold)
                & valid
                & baseline.notna()
            ).to_numpy()

            events = detect_events(
                candidate,
                flux,
                time_hours
            )

            print(f"\nElectron Band {band + 1}")
            print("Valid samples:", int(valid.sum()))
            print("Threshold:", round(threshold, 3))
            print("Candidate samples:", int(candidate.sum()))
            print("Candidate events:", len(events))

            all_summaries.append({
                "date": date,
                "energy_band": band + 1,
                "valid_samples": int(valid.sum()),
                "threshold": threshold,
                "candidate_samples": int(candidate.sum()),
                "candidate_percentage": (
                    candidate.sum() / len(candidate) * 100
                ),
                "event_count": len(events)
            })

            for number, event in enumerate(events, 1):

                all_events.append({
                    "date": date,
                    "energy_band": band + 1,
                    "event_number": number,
                    **event
                })

                print(
                    f"Event {number}: "
                    f"{event['start_hour']:.3f} - "
                    f"{event['end_hour']:.3f} h | "
                    f"Active samples: {event['active_samples']} | "
                    f"Peak: {event['peak_flux']:.2f}"
                )


summary_df = pd.DataFrame(all_summaries)
events_df = pd.DataFrame(all_events)

summary_df.to_csv(
    OUTPUT / "multiday_detector_summary.csv",
    index=False
)

events_df.to_csv(
    OUTPUT / "multiday_candidate_events.csv",
    index=False
)

print("\n========== MULTI-DAY SUMMARY ==========")
print(summary_df.to_string(index=False))

print("\nSaved:")
print(OUTPUT / "multiday_detector_summary.csv")
print(OUTPUT / "multiday_candidate_events.csv")
