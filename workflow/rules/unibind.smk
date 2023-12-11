"""d"""


# Parameters
INSTALL_DIR = config["TFBS-ACTIVITY"]["INSTALL_DIR"]
PROCESS_DIR = config["TFBS-ACTIVITY"]["PROCESS_DIR"]
UNIBIND_URL = config["TFBS-ACTIVITY"]["UNIBIND_URL"]


rule all:
    input:
        INSTALL_DIR
        + "/damo_hg38_all_TFBS/ENCSR000AHD.MCF7_Invasive_ductal_breast_carcinoma.CTCF.MA0139.1.damo.bed",


rule download_unibind:
    message:
        """
        - Downloads all DAMOS from UniBind database
        """
    output:
        temp(INSTALL_DIR + "/data/unibind/damo_hg39_all_TFBS.tar.gz"),
    params:
        url=UNIBIND_URL,
    shell:
        """
        curl -o {output} {params.url}
        """


rule decompress_unibind:
    message:
        """
        - Decompresses UniBind tarbell
        """
    input:
        rules.download_unibind.output,
    output:
        temp(INSTALL_DIR + "/data/unibind/damo_hg39_all_TFBS.tar"),
    shell:
        """
        gunzip {input}
        """


rule unpack_unibind:
    message:
        """
        - Unpacks UniBind tarbell
        """
    input:
        rules.decompress_unibind.output,
    output:
        directory(INSTALL_DIR + "/data/unibind/damo_hg39_all_TFBS"),
    shell:
        """
        tar -xvf {input}
        """


rule flatten_dir:
    message:
        """
        - Flattens unpacked tarbell. Not great - for the moment just check for creation of first in OP.
        Also removes all the empty dirs.
        """
    input:
        rules.unpack_unibind.output,
    output:
        INSTALL_DIR
        + "/damo_hg38_all_TFBS/ENCSR000AHD.MCF7_Invasive_ductal_breast_carcinoma.CTCF.MA0139.1.damo.bed",
    shell:
        """
        find {input} -mindepth 2 -type f -exec mv -t {input} -i '{{}}' + &&
         ind {input} -type d -empty -delete
        """

rule organize_damos:
    message:
        """
        Organizes damos into dirs by TF
        """
    input:
        rules.flatten_dir.output,
    output:
        INSTALL_DIR
        + "/damo_hg38_all_TFBS/MA0139.1/ENCSR000AHD.MCF7_Invasive_ductal_breast_carcinoma.CTCF.MA0139.1.damo.bed",
    shell:
        """
        find {input} -mindepth 2 -type f -exec mv -t {input} -i '{}' + &&
        find {input} -type d -empty -delete
        """
