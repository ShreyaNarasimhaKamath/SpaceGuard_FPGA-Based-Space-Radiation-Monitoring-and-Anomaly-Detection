
import pandas as pd
import numpy as np

FILE = "results/evaluation/cross_channel_results.csv"

df = pd.read_csv(FILE)

event = df[df["event"] == 2].copy()

summary = (
    event.groupby("energy_band")
    .agg(
        telescope_count=("telescope", "nunique"),
        median_fold_change=("fold_change", "median"),
        minimum_fold_change=("fold_change", "min"),
        maximum_fold_change=("fold_change", "max"),
        median_percent_change=("percent_change", "median"),
        valid_measurements=("event_valid_samples", "sum")
    )
    .reset_index()
)

print("\n========== EVENT 2: ENERGY-BAND SUMMARY ==========")
print(summary.to_string(index=False))

print("\n========== PER-TELESCOPE FOLD CHANGE ==========")

table = event.pivot(
    index="telescope",
    columns="energy_band",
    values="fold_change"
)

print(table.to_string())
