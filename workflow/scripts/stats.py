"""Utilities for the activity module"""

from re import findall
from subprocess import run, PIPE
from numpy import log10
from Bio.Seq import Seq
from pandas import DataFrame, Series, read_csv, merge
from Bio import motifs
from Bio.motifs.matrix import PositionSpecificScoringMatrix
from scipy.stats.mstats import winsorize
from statsmodels.stats.proportion import proportion_confint

# Snakemake
ACTIVITY = snakemake.input[0]  # type: ignore
PVALUES = snakemake.input[1]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore
PROFILE = snakemake.params.profile  # type: ignore


def read_activity(filepath: str) -> DataFrame:
    """Reads activity mapping into DF"""
    # Return matrix, set names
    fields = ["score", "unibind_count", "motif_count", "activity"]
    return read_csv(filepath, header=None, sep="\t", names=fields, engine="c")


def read_pvals(filepath: str) -> DataFrame:
    """Reads pval distribution into DF"""
    # Return matrix, set names
    fields = ["score", "p-value", "percentage"]
    return read_csv(
        filepath, header=None, delim_whitespace=True, names=fields, engine="c"
    )


def main() -> None:
    """Merges activity map with percentages"""
    # Read inputs
    activity = read_activity(ACTIVITY)
    pvalues = read_pvals(PVALUES)

    # Merge
    activity = merge(activity, pvalues, on="score", how="left")

    # Standard error on proportion
    activity["ci95_lbound"], activity["ci95_rbound"] = proportion_confint(
        activity["unibind_count"], activity["motif_count"]
    )

    # Winsorize activity and SE, 0% on left tail, 35% of values on the right tale which is about 90percentile PWM score
    activity["activity_winsor"] = winsorize(activity["activity"], limits=(0.00, 0.35))
    activity["ci95_lbound_winsor"] = winsorize(
        activity["ci95_lbound"], limits=(0.00, 0.35)
    )
    activity["ci95_rbound_winsor"] = winsorize(
        activity["ci95_rbound"], limits=(0.00, 0.35)
    )

    # Add TF flag
    activity["profile"] = PROFILE

    # Write out
    activity.to_csv(OUTPUT, sep="\t", index=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
