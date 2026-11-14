
import pandas as pd
import numpy as np
from scipy.optimize import linear_sum_assignment

FILE = "results/evaluation/independent_channel_events.csv"

GAP_SETTINGS = [60, 120, 180, 300, 600]

df = pd.read_csv(FILE)


def merge_events(channel, gap_seconds):
    """Merge same-channel intervals separated by short gaps."""

    intervals = sorted(
        zip(
            channel["start_hour"] * 3600,
            channel["end_hour"] * 3600
        )
    )

    if not intervals:
        return []

    merged = []
    start, end = intervals[0]

    for next_start, next_end in intervals[1:]:

        gap = next_start - end

        if gap <= gap_seconds:
            end = max(end, next_end)
        else:
            merged.append((start, end))
            start, end = next_start, next_end

    merged.append((start, end))

    return merged


def calculate_iou(a, b):
    """Intersection over Union for two time intervals."""

    intersection = max(
        0,
        min(a[1], b[1]) - max(a[0], b[0])
    )

    union = (
        (a[1] - a[0])
        + (b[1] - b[0])
        - intersection
    )

    if union <= 0:
        return 0.0

    return intersection / union


print("\n========== EVENT GAP SENSITIVITY ==========")

summary = []

for gap in GAP_SETTINGS:

    b1 = merge_events(
        df[df["energy_band"] == 1],
        gap
    )

    b2 = merge_events(
        df[df["energy_band"] == 2],
        gap
    )

    # Build IoU matrix
    matrix = np.zeros((len(b1), len(b2)))

    for i, a in enumerate(b1):
        for j, b in enumerate(b2):
            matrix[i, j] = calculate_iou(a, b)

    # Globally optimal one-to-one assignment
    if matrix.size:
        rows, cols = linear_sum_assignment(-matrix)

        matches = [
            matrix[i, j]
            for i, j in zip(rows, cols)
            if matrix[i, j] > 0
        ]
    else:
        matches = []

    matched_count = len(matches)

    mean_iou = (
        float(np.mean(matches))
        if matches else 0.0
    )

    summary.append({
        "merge_gap_seconds": gap,
        "band1_episodes": len(b1),
        "band2_episodes": len(b2),
        "matched_episode_pairs": matched_count,
        "mean_matched_iou": mean_iou
    })

    print(
        f"\nMerge gap: {gap} seconds"
        f"\n  Band 1 episodes: {len(b1)}"
        f"\n  Band 2 episodes: {len(b2)}"
        f"\n  Matched pairs: {matched_count}"
        f"\n  Mean matched IoU: {mean_iou:.3f}"
    )

result = pd.DataFrame(summary)

result.to_csv(
    "results/evaluation/event_gap_sensitivity.csv",
    index=False
)

print("\n========== SUMMARY TABLE ==========")
print(result.to_string(index=False))

print(
    "\nSaved: "
    "results/evaluation/event_gap_sensitivity.csv"
)
