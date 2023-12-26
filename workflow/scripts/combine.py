"""d"""

from pypdf import PdfMerger

# Snakemake
PDFS = snakemake.input  # type: ignore
OUTPUT = snakemake.output[0]  # type: ignore


def main():
    pdfs = PDFS

    print("BARK")
    print(pdfs)

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
