"""d"""

from pandas import DataFrame, read_csv

# Snakemake
ACTIVITY = snakemake.input[0]  # type: ignore
WINDOW = snakemake.params.window  # type: ignore
THRESH = snakemake.params.threshold  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore


###
# Functions
###


def read_activity(
    filepath: str,
    fields: list = [
        "score",
        "unibind_count",
        "motif_count",
        "activity",
        "pval",
        "perc",
        "ci95_lbound",
        "ci95_rbound",
        "profile",
        "dataset",
    ],
) -> DataFrame:
    """Reads activity mapping into DF"""
    # Return matrix, set names
    return read_csv(
        filepath,
        sep="\t",
        engine="c",
        dtype=dict(
            zip(
                fields,
                [int, int, int, float, float, float, float, float, str, str],
            )
        ),
    )


def calculate_ci_width_and_rolling_avg(
    df, ci_lbound_col, ci_rbound_col, window_size, threshold
):
    # Calculate width of the 95ci
    df["ci_width"] = df[ci_rbound_col] - df[ci_lbound_col]

    # Calculate rolling average of ci_width
    df["ci_width_roll"] = df["ci_width"].rolling(window=window_size).mean()

    # Find last index where rolling average of ci_width is < threshold
    # if the map isn't long enough to get this, return the last index
    try:
        cap_index = df[df["ci_width_roll"] < threshold].index[-1]
    except IndexError:
        cap_index = df.index[-1]
    return df, cap_index


def main() -> None:
    """d"""
    # Read in activity
    activity = read_activity(ACTIVITY)

    # Calculate CI width and rolling average
    activity, cap_index = calculate_ci_width_and_rolling_avg(
        activity, "ci95_lbound", "ci95_rbound", WINDOW, THRESH
    )

    # Update to values at the cap index
    columns_to_update = ["activity", "ci95_lbound", "ci95_rbound"]
    for column in columns_to_update:
        activity.loc[cap_index:, column] = activity.loc[cap_index, column]

    # Write to file
    activity.to_csv(OUTPUT, sep="\t", index=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
