"""Scores matched sequence from UniBind DAMO files"""

from Bio.Seq import Seq
from pandas import DataFrame, read_csv
from Bio.motifs.matrix import PositionSpecificScoringMatrix


# Snakemake
BED = snakemake.input[0]  # type: ignore
PWM = snakemake.input[1]  # type: ignore
OUT = snakemake.output[0]  # type: ignore


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
    fields = ["chrm", "pos0", "pos1", "seq", "strand"]
    dtypes = {
        "chrm": str,
        "pos0": int,
        "pos1": int,
        "seq": str,
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
    return pssm.calculate(seq)  # type: ignore


def main():
    # Read in sites data and pwm
    bed = read_damo(BED)
    pssm = setup_pssm(PWM)
    # Score motif sequence and write out
    bed["."] = bed.apply(lambda row: score_sequence(row.seq, pssm), axis=1)
    bed.to_csv(OUT, sep="\t", index=False, header=False)


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
