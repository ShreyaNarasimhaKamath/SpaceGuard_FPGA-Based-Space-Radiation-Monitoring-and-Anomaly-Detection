
import pandas as pd
import matplotlib.pyplot as plt

FILE = "results/evaluation/cross_channel_results.csv"

df = pd.read_csv(FILE)

selected_events = [1, 2, 7, 8, 9, 12]

# Use Telescope 1 for both energy bands
subset = df[
    (df["event"].isin(selected_events))
    & (df["telescope"] == 1)
    & (df["energy_band"].isin([1, 2]))
]

table = subset.pivot(
    index="event",
    columns="energy_band",
    values="fold_change"
)

table = table.reindex(selected_events)
table.columns = ["Band 1", "Band 2"]

print("\n========== MULTI-EVENT CHANNEL COMPARISON ==========")
print(table.to_string())

print("\n========== BAND 2 / BAND 1 RESPONSE ==========")
table["Band 2 / Band 1"] = table["Band 2"] / table["Band 1"]
print(table.to_string())

# Plot comparison
ax = table[["Band 1", "Band 2"]].plot(
    kind="bar",
    figsize=(11, 6)
)

ax.set_xlabel("Candidate event")
ax.set_ylabel("Fold change relative to pre-event interval")
ax.set_title("Band 1 vs Band 2: Multi-Event Comparison")
ax.grid(axis="y")
plt.xticks(rotation=0)
plt.tight_layout()

plt.savefig(
    "results/plots/multi_event_channel_comparison.png",
    dpi=200
)

plt.show()
