
import pandas as pd

FILE = "results/evaluation/independent_channel_events.csv"

df = pd.read_csv(FILE)

band1 = df[df["energy_band"] == 1].copy()
band2 = df[df["energy_band"] == 2].copy()

print("\n========== BAND 1 EVENTS ==========")

for _, e in band1.iterrows():
    print(
        f"Event {int(e['event_number']):2d}: "
        f"{e['start_hour']:.3f} - {e['end_hour']:.3f} h | "
        f"{e['duration_seconds']} s"
    )

print("\n========== BAND 2 EVENTS ==========")

for _, e in band2.iterrows():
    print(
        f"Event {int(e['event_number']):2d}: "
        f"{e['start_hour']:.3f} - {e['end_hour']:.3f} h | "
        f"{e['duration_seconds']} s"
    )

print("\n========== TEMPORAL OVERLAPS ==========")

for _, a in band1.iterrows():
    for _, b in band2.iterrows():

        overlap = max(
            0,
            min(a["end_hour"], b["end_hour"])
            - max(a["start_hour"], b["start_hour"])
        )

        gap = max(
            0,
            a["start_hour"] - b["end_hour"],
            b["start_hour"] - a["end_hour"]
        )

        if overlap > 0:
            print(
                f"B1 Event {int(a['event_number'])} "
                f"<-> B2 Event {int(b['event_number'])} | "
                f"Overlap: {overlap * 3600:.1f} s"
            )
        elif gap <= 300 / 3600:
            print(
                f"B1 Event {int(a['event_number'])} "
                f"<-> B2 Event {int(b['event_number'])} | "
                f"Gap: {gap * 3600:.1f} s"
            )
