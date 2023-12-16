"""D"""

from pandas import read_csv, DataFrame

# Snakemake
METADATA = snakemake.input[0]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore


def read_metadata(filepath: str) -> DataFrame:
    """Returns DF of metadata"""
    return read_csv(filepath, sep="\t", engine="c")


def main() -> None:
    """D"""
    # Read input
    metadata = read_metadata(METADATA)

    # group
    grouped = metadata.groupby(["percentage"])
    columns = ["score", "mean", "size", "sem"]

    # Aggregate
    summary = grouped.agg({"activity": ["mean", "size", "sem"]}).reset_index()

    # Set columns
    summary.columns = columns

    # 95% CI
    summary["ci95_lbound"] = summary["mean"] - (1.96 * summary["sem"])
    summary["ci95_rbound"] = summary["mean"] + (1.96 * summary["sem"])

    # Save results
    summary.to_csv(OUTPUT, sep="\t", index=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
