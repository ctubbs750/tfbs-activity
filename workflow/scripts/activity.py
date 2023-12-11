"""Utilities for the activity module"""

from re import findall
from subprocess import run, PIPE
from numpy import log10
from Bio.Seq import Seq
from pandas import DataFrame, Series, read_csv
from Bio import motifs
from Bio.motifs.matrix import PositionSpecificScoringMatrix
from scipy.stats.mstats import winsorize
from statsmodels.stats.proportion import proportion_confint

###
# Functions
###


def setup_pssm(pwm_filepath: str) -> PositionSpecificScoringMatrix:
    """D"""
    # read in, convert to dict
    pwm = read_csv(
        pwm_filepath, delim_whitespace=True, header=None, names=["A", "C", "G", "T"]
    ).to_dict(orient="list")

    # Convert to PSSM
    pssm = PositionSpecificScoringMatrix(alphabet="ACGT", values=pwm)
    return pssm


def read_damo(filepath: str) -> DataFrame:
    """Reads UniBind damo format file into pandas DF"""
    # Setup fast readin
    fields = ["chrm", "pos0", "pos1", "seq", "damo", "strand"]
    dtypes = {
        "chrm": str,
        "pos0": int,
        "pos1": int,
        "seq": str,
        "damo": str,
        "strand": str,
    }
    return read_csv(
        filepath, sep="\t", header=None, names=fields, dtype=dtypes, engine="pyarrow"
    )


def score_sequence(sequence: str, pssm: PositionSpecificScoringMatrix) -> float:
    """Scores input sequence on match to motif"""
    # Convert to sequence obj, assume sequence is strand corrected
    seq = Seq(sequence)

    # Sum to return raw score.
    return str(seq), pssm.calculate(seq)  # type: ignore


def run_subprocess(command: str) -> list:
    """Executes subprocess and returns decoded output"""
    # Execute process and capture stdout
    process = run([command], shell=True, stdout=PIPE, stderr=PIPE, check=False)

    # Return stdout. Last field is empty in these prcoesses
    return process.stdout.decode("utf-8").split("\n")[:-1]


def calculate_pwm(jaspar_filepath: str, output: str) -> None:
    """Converts dowloaded PFM to PWM"""
    # Connect to JASPAR motif matrix
    with open(jaspar_filepath) as handle:
        motif = motifs.read(handle, "jaspar")

        # Calculate pseudocounts
        motif.pseudocounts = motifs.jaspar.calculate_pseudocounts(motif)  # type: ignore

        # Caclulate PWM
        pwm = list(map(list, zip(*[motif.pssm[nt] for nt in "ACGT"])))  # type: ignore

        # Write out PWM
        with open(output, mode="w", encoding="utf-8") as op:
            for line in pwm:
                line = " ".join(
                    ["{:7d}".format(round(j * 100)) for j in line]
                )  # type: ignore
                op.write(f"{line}\n")


def pval_distribution(
    pwm_path: str,
    matrix_prob: str,
    rthresh: float,
    pvals_out: str,
    coeff_out: str,
) -> None:
    """Generates distribution of PWM scores for pvalue determination and sets CUTOFF"""
    # Calculate distribution of PWM scores
    cmd = f"{matrix_prob} {pwm_path}"
    out = run_subprocess(cmd)

    # Process output into dictionary: {raw_score[adjusted_pval]}
    probs = {}
    with open(pvals_out, mode="w", encoding="utf-8") as output:
        for line in out:
            # Strip and clean data
            data = findall("(\S+)", line)  # type: ignore

            # Process
            score = int(data[0])
            perc = float(data[2][:-1])
            pval = float(data[1])

            # Transform
            adj_pval = int(log10(pval) * 1000 / -10)

            # Update
            probs[score] = adj_pval

            # Get PWM score cutoff and write
            # if pval < pthresh and perc >= rthresh * 100:
            if perc >= rthresh * 100:
                with open(coeff_out, mode="w", encoding="utf-8") as cutoff:
                    cutoff.write(str(score))
            # Write out pvals
            output.write(
                "\t".join([str(score), str(adj_pval), str(pval), str(perc)]) + "\n"
            )


def proportion_ci(counts: Series, totals: Series) -> tuple:
    """Return proprtion 95 CI"""
    return proportion_confint(
        counts,
        totals,
        alpha=0.05,
        method="normal",
    )


def winsorize_array(array: Series, l_limit: int, r_limit: int) -> Series:
    """Winsorizes input array"""
    return winsorize(array, limits=(l_limit, r_limit))
