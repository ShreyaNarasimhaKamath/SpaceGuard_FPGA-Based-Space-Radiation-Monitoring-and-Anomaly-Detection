
import pandas as pd

FILE = "results/evaluation/independent_channel_events.csv"

df = pd.read_csv(FILE)

band1 = df[df["energy_band"] == 1]
band2 = df[df["energy_band"] == 2]

matches = []

for _, a in band1.iterrows():

    for _, b in band2.iterrows():

        overlap_start = max(
            a["start_hour"],
            b["start_hour"]
        )

        overlap_end = min(
            a["end_hour"],
            b["end_hour"]
        )

        overlap = max(
            0,
            (overlap_end - overlap_start) * 3600
        )

        if overlap <= 0:
            continue

        duration1 = a["duration_seconds"]
        duration2 = b["duration_seconds"]

        union = duration1 + duration2 - overlap

        iou = overlap / union

        coverage1 = overlap / duration1
        coverage2 = overlap / duration2

        matches.append({
            "band1_event": int(a["event_number"]),
            "band2_event": int(b["event_number"]),
            "overlap_seconds": round(overlap, 1),
            "band1_coverage_pct": round(coverage1 * 100, 1),
            "band2_coverage_pct": round(coverage2 * 100, 1),
            "iou": round(iou, 3)
        })

result = pd.DataFrame(matches)

if result.empty:
    print("No overlapping events.")
else:
    result = result.sort_values(
        "iou",
        ascending=False
    )

    print("\n========== EVENT OVERLAP METRICS ==========")
    print(result.to_string(index=False))

    result.to_csv(
        "results/evaluation/event_overlap_metrics.csv",
        index=False
    )

    print("\nSaved to results/evaluation/event_overlap_metrics.csv")
