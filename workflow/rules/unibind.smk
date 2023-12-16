"""d"""


# Parameters
INSTALL_DIR = config["TFBS-ACTIVITY"]["INSTALL_DIR"]
PROCESS_DIR = config["TFBS-ACTIVITY"]["PROCESS_DIR"]
UNIBIND_URL = config["TFBS-ACTIVITY"]["UNIBIND_URL"]


rule all:
    input:
        INSTALL_DIR + "/damo_hg38_all_TFBS.unique_pwms.txt"


rule download_unibind:
    message:
        """
        - Downloads all DAMOS from UniBind database
        """
    output:
        temp(INSTALL_DIR + "/damo_hg39_all_TFBS.tar.gz"),
    params:
        url=UNIBIND_URL,
    shell:
        """
        curl -o {output} {params.url}
        """


checkpoint unpack_unibind:
    message:
        """
        - Unpacks UniBind tarbell
        """
    input:
        rules.download_unibind.output,
    output:
        directory(INSTALL_DIR + "/damo_hg38_all_TFBS"),
    params:
        outdir=INSTALL_DIR,
    shell:
        """
        tar -xzf {input} -C {params.outdir}
        """


checkpoint flatten_dir:
# TODO: Works but doesn't register with snakemake upon completion.
    message:
        """
        - Flattens unpacked tarbell. Not great - for the moment just check for creation of first in OP.
        Also removes all the empty dirs.
        """
    input:
        INSTALL_DIR + "/damo_hg38_all_TFBS",
    output:
        directory(INSTALL_DIR + "/damo_hg38_all_TFBS_unpacked"),
    shell:
        """
        find {input} -mindepth 2 -type f -exec mv -t {output} -i '{{}}' + &&
        find {input} -type d -empty -delete
        """

rule unique_pwms:
    message:
        """
        get list of unique profiles in unibind
        """
    input:
        rules.flatten_dir.output
    output:
        INSTALL_DIR + "/damo_hg38_all_TFBS.unique_pwms.txt"
    shell:
        """
        ls {input} | cut -d"." -f4-5 | sort -u > {output}
        """
##
