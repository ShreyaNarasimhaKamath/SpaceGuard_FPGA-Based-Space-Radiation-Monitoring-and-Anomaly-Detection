
import pandas as pd
import numpy as np
from scipy.optimize import linear_sum_assignment

FILE = "results/evaluation/multiday_candidate_events.csv"

df = pd.read_csv(FILE)

def iou(a, b):
    intersection = max(
        0,
        min(a["end_hour"], b["end_hour"])
        - max(a["start_hour"], b["start_hour"])
    )

    union = (
        (a["end_hour"] - a["start_hour"])
        + (b["end_hour"] - b["start_hour"])
        - intersection
    )

    return intersection / union if union > 0 else 0


summaries = []
all_matches = []

for date in sorted(df["date"].unique()):

    day = df[df["date"] == date]

    b1 = day[day["energy_band"] == 1].reset_index(drop=True)
    b2 = day[day["energy_band"] == 2].reset_index(drop=True)

    matrix = np.zeros((len(b1), len(b2)))

    for i, a in b1.iterrows():
        for j, b in b2.iterrows():
            matrix[i, j] = iou(a, b)

    if matrix.size:
        rows, cols = linear_sum_assignment(-matrix)

        matches = [
            (i, j, matrix[i, j])
            for i, j in zip(rows, cols)
            if matrix[i, j] > 0
        ]
    else:
        matches = []

    matched_b1 = {i for i, j, score in matches}
    matched_b2 = {j for i, j, score in matches}

    print(f"\n========== {date} ==========")
    print("Band 1 events:", len(b1))
    print("Band 2 events:", len(b2))
    print("One-to-one overlapping pairs:", len(matches))
    print("Band 1 unmatched:", len(b1) - len(matched_b1))
    print("Band 2 unmatched:", len(b2) - len(matched_b2))

    for i, j, score in matches:

        overlap_start = max(
            b1.loc[i, "start_hour"],
            b2.loc[j, "start_hour"]
        )

        overlap_end = min(
            b1.loc[i, "end_hour"],
            b2.loc[j, "end_hour"]
        )

        overlap_seconds = (overlap_end - overlap_start) * 3600

        all_matches.append({
            "date": date,
            "band1_event": int(b1.loc[i, "event_number"]),
            "band2_event": int(b2.loc[j, "event_number"]),
            "iou": score,
            "overlap_seconds": overlap_seconds
        })

        print(
            f"B1 Event {int(b1.loc[i, 'event_number'])} "
            f"<-> B2 Event {int(b2.loc[j, 'event_number'])} "
            f"| IoU: {score:.3f} "
            f"| Overlap: {overlap_seconds:.1f} s"
        )

    summaries.append({
        "date": date,
        "band1_events": len(b1),
        "band2_events": len(b2),
        "matched_pairs": len(matches),
        "band1_unmatched": len(b1) - len(matched_b1),
        "band2_unmatched": len(b2) - len(matched_b2),
        "mean_iou": (
            np.mean([m[2] for m in matches])
            if matches else 0
        )
    })

summary = pd.DataFrame(summaries)
matches_df = pd.DataFrame(all_matches)

summary.to_csv(
    "results/evaluation/multiday_matching_summary.csv",
    index=False
)

matches_df.to_csv(
    "results/evaluation/multiday_matching_details.csv",
    index=False
)

print("\n========== FINAL SUMMARY ==========")
print(summary.to_string(index=False))
