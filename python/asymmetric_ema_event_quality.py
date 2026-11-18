
from pathlib import Path

import numpy as np
import pandas as pd


# ============================================================
# CONFIGURATION
# ============================================================

ROOT = Path(__file__).resolve().parents[1]

RESULTS_DIR = ROOT / "results" / "evaluation"

MATCH_FILE = RESULTS_DIR / "asymmetric_ema_event_matches.csv"
SUMMARY_FILE = RESULTS_DIR / "asymmetric_ema_sensitivity.csv"

OUTPUT_FILE = RESULTS_DIR / "asymmetric_ema_event_quality.csv"

CONFIGURATIONS = ["A", "B", "C", "D", "E"]
IOU_THRESHOLDS = [0.30, 0.50]

# ============================================================
# LOAD DATA
# ============================================================

def load_data():
    if not MATCH_FILE.exists():
        raise FileNotFoundError(MATCH_FILE)

    if not SUMMARY_FILE.exists():
        raise FileNotFoundError(SUMMARY_FILE)

    matches = pd.read_csv(MATCH_FILE)
    summary = pd.read_csv(SUMMARY_FILE)

    required_match_columns = {
        "date",
        "band",
        "configuration",
        "iou",
        "onset_difference_seconds",
        "duration_difference_seconds",
    }

    required_summary_columns = {
        "date",
        "band",
        "configuration",
        "reference_events",
        "candidate_events",
    }

    if not required_match_columns.issubset(matches.columns):
        raise ValueError(
            f"Missing match columns: "
            f"{required_match_columns - set(matches.columns)}"
        )

    if not required_summary_columns.issubset(summary.columns):
        raise ValueError(
            f"Missing summary columns: "
            f"{required_summary_columns - set(summary.columns)}"
        )

    return matches, summary


# ============================================================
# EVENT QUALITY METRICS
# ============================================================

def calculate_metrics(matches, summary, iou_threshold):
    """
    Calculate event precision, recall and F1.

    The existing event matcher uses one-to-one Hungarian
    assignment and retains only positive-IoU matches.

    Therefore, filtering its matches by an IoU threshold
    gives the number of accepted one-to-one matches.
    """

    accepted = matches[
        matches["iou"] >= iou_threshold
    ]

    rows = []

    for configuration in CONFIGURATIONS:
        for band in [1, 2]:

            group_summary = summary[
                (summary["configuration"] == configuration)
                & (summary["band"] == band)
            ]

            group_matches = accepted[
                (accepted["configuration"] == configuration)
                & (accepted["band"] == band)
            ]

            reference_events = int(
                group_summary["reference_events"].sum()
            )

            candidate_events = int(
                group_summary["candidate_events"].sum()
            )

            true_matches = len(group_matches)

            false_positives = candidate_events - true_matches
            false_negatives = reference_events - true_matches

            precision = (
                true_matches / candidate_events
                if candidate_events else 0.0
            )

            recall = (
                true_matches / reference_events
                if reference_events else 0.0
            )

            f1 = (
                2 * precision * recall / (precision + recall)
                if precision + recall else 0.0
            )

            if true_matches:
                mean_iou = group_matches["iou"].mean()

                mean_onset_error = (
                    group_matches["onset_difference_seconds"].mean()
                )

                median_absolute_onset_error = (
                    group_matches["onset_difference_seconds"]
                    .abs()
                    .median()
                )

                mean_duration_error = (
                    group_matches["duration_difference_seconds"].mean()
                )

                mean_absolute_duration_error = (
                    group_matches["duration_difference_seconds"]
                    .abs()
                    .mean()
                )

            else:
                mean_iou = np.nan
                mean_onset_error = np.nan
                median_absolute_onset_error = np.nan
                mean_duration_error = np.nan
                mean_absolute_duration_error = np.nan

            rows.append({
                "configuration": configuration,
                "band": band,
                "iou_threshold": iou_threshold,
                "reference_events": reference_events,
                "candidate_events": candidate_events,
                "matched_events": true_matches,
                "false_positives": false_positives,
                "false_negatives": false_negatives,
                "precision": precision,
                "recall": recall,
                "f1_score": f1,
                "mean_matched_iou": mean_iou,
                "mean_onset_error_seconds": mean_onset_error,
                "median_absolute_onset_error_seconds":
                    median_absolute_onset_error,
                "mean_duration_error_seconds": mean_duration_error,
                "mean_absolute_duration_error_seconds":
                    mean_absolute_duration_error,
            })

    return pd.DataFrame(rows)


# ============================================================
# DATE-LEVEL QUALITY METRICS
# ============================================================

def calculate_date_metrics(matches, summary, iou_threshold):
    accepted = matches[
        matches["iou"] >= iou_threshold
    ]

    rows = []

    for _, summary_row in summary.iterrows():

        configuration = summary_row["configuration"]
        band = summary_row["band"]
        date = summary_row["date"]

        group_matches = accepted[
            (accepted["configuration"] == configuration)
            & (accepted["band"] == band)
            & (accepted["date"] == date)
        ]

        reference_events = int(summary_row["reference_events"])
        candidate_events = int(summary_row["candidate_events"])
        matched_events = len(group_matches)

        precision = (
            matched_events / candidate_events
            if candidate_events else 0.0
        )

        recall = (
            matched_events / reference_events
            if reference_events else 0.0
        )

        f1 = (
            2 * precision * recall / (precision + recall)
            if precision + recall else 0.0
        )

        rows.append({
            "date": date,
            "band": band,
            "configuration": configuration,
            "iou_threshold": iou_threshold,
            "reference_events": reference_events,
            "candidate_events": candidate_events,
            "matched_events": matched_events,
            "precision": precision,
            "recall": recall,
            "f1_score": f1,
        })

    return pd.DataFrame(rows)


# ============================================================
# MAIN
# ============================================================

def main():
    matches, summary = load_data()

    all_metrics = []
    date_metrics = []

    for threshold in IOU_THRESHOLDS:

        print("\n" + "=" * 65)
        print(f"EVENT QUALITY AT IoU >= {threshold:.2f}")
        print("=" * 65)

        metrics = calculate_metrics(
            matches, summary, threshold
        )

        per_date = calculate_date_metrics(
            matches, summary, threshold
        )

        all_metrics.append(metrics)
        date_metrics.append(per_date)

        print(
            metrics[
                [
                    "configuration",
                    "band",
                    "reference_events",
                    "candidate_events",
                    "matched_events",
                    "precision",
                    "recall",
                    "f1_score",
                    "mean_matched_iou",
                    "median_absolute_onset_error_seconds",
                    "mean_absolute_duration_error_seconds",
                ]
            ].to_string(index=False)
        )

    final_metrics = pd.concat(
        all_metrics, ignore_index=True
    )

    final_date_metrics = pd.concat(
        date_metrics, ignore_index=True
    )

    # Save aggregate and date-level results separately.
    output_file = OUTPUT_FILE

    final_metrics.to_csv(output_file, index=False)

    date_output = (
        RESULTS_DIR / "asymmetric_ema_event_quality_by_date.csv"
    )

    final_date_metrics.to_csv(date_output, index=False)

    print("\n" + "=" * 65)
    print("EVENT QUALITY EVALUATION COMPLETE")
    print("=" * 65)

    print(f"Aggregate metrics: {output_file}")
    print(f"Date-level metrics: {date_output}")


if __name__ == "__main__":
    main()
