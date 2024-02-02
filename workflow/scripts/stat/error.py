"""d"""

from typing import List
from pandas import DataFrame, read_csv, merge
from statsmodels.stats.proportion import proportion_confint

# Snakemake
ACTIVITY = snakemake.input[0]  # type: ignore
PVALUES = snakemake.input[1]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore
PROFILE = snakemake.params[0]  # type: ignore
DATASET = snakemake.params[1]  # type: ignore

# ------------- #
# Functions     #
# ------------- #



def read_activity(filepath: str, fields: List[str] = None) -> DataFrame:
    """
    Reads activity mapping into DataFrame

    Parameters:
    filepath (str): Path to the file to read
    fields (List[str]): List of column names to use. Defaults to ["score", "unibind_count", "motif_count", "activity"].

    Returns:
    DataFrame: The read data
    """
    if fields is None:
        fields = ["score", "unibind_count", "motif_count", "activity"]

    return read_csv(
        filepath,
        header=None,
        sep="\t",
        names=fields,
        engine="c",
        dtype=dict(zip(fields, [int, int, int, float])),
    )


def read_profile_pvals(filepath: str, labels: List[str] = None) -> DataFrame:
    """
    Reads pval data for profile

    Parameters:
    filepath (str): Path to the file to read
    labels (List[str]): List of column names to use. Defaults to ["score", "pval", "perc"].

    Returns:
    DataFrame: The read data
    """
    if labels is None:
        labels = ["score", "pval", "perc"]

    return read_csv(
        filepath,
        delim_whitespace=True,
        header=None,
        names=labels,
        dtype=dict(zip(labels, [int, float, float])),
    )


def main(
    activity_file: str, pvalues_file: str, output_file: str, profile: str, dataset: str
) -> None:
    """
    Merges activity map with percentages and writes the result to a file.

    Parameters:
    activity_file (str): Path to the activity file
    pvalues_file (str): Path to the pvalues file
    output_file (str): Path to the output file
    profile (str): The profile name
    dataset (str): The dataset name
    """
    # Read inputs
    activity = read_activity(activity_file)
    pvalues = read_profile_pvals(pvalues_file)

    # Merge
    activity = merge(activity, pvalues, on="score", how="left")

    # Standard error on proportion
    activity["ci95_lbound"], activity["ci95_rbound"] = proportion_confint(
        activity["unibind_count"], activity["motif_count"]
    )

    # Add TF flag
    activity["profile"] = profile
    activity["dataset"] = dataset

    # Write out
    activity.to_csv(output_file, sep="\t", index=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main(ACTIVITY, PVALUES, OUTPUT, PROFILE, DATASET)
