"""d"""

from pypdf import PdfMerger

# Snakemake
PDFS = snakemake.input[0]  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore


def main():
    pdfs = PDFS

    merger = PdfMerger()

    for pdf in pdfs:
        merger.append(pdf)

    merger.write(OUTPUT)
    merger.close()


# ------------- #
# Main          #
# ------------- #

if __name__ == "__main__":
    main()
