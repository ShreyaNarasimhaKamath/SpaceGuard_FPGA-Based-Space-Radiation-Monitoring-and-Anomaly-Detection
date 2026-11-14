
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
path = ROOT / "results/evaluation/fair_event_comparison.csv"

df = pd.read_csv(path)

subset = df[
    (df["date"] == 20240103)
    & (df["band"] == 2)
]

for method in ["rolling_median", "moving_average", "ema_always"]:

    events = subset[subset["method"] == method]

    print(f"\n{'=' * 55}")
    print(method.upper())
    print(f"{'=' * 55}")

    if events.empty:
        print("No events detected")
        continue

    print(
        events[
            [
                "event_number",
                "start_sample",
                "end_sample",
                "duration_seconds",
                "active_samples",
            ]
        ].to_string(index=False)
    )

    print("\nEvent times (hours):")

    for _, event in events.iterrows():
        start = event["start_sample"] / 3600
        end = event["end_sample"] / 3600

        print(
            f"Event {int(event['event_number'])}: "
            f"{start:.3f} h – {end:.3f} h"
        )

