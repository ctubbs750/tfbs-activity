# Parameters
INSTALL_DIR = config["TFBS-ACTIVITY"]["INSTALL_DIR"]
PROCESS_DIR = config["TFBS-ACTIVITY"]["PROCESS_DIR"]
PROFILES = config["PWM-SCAN"]["TARGETS"]


# WC constraints - JASPAR matrix format
wildcard_constraints:
    PROFILE="[a-zA-Z\d]{6}.{1}\d{1}",


rule all:
    input:
        PROCESS_DIR + "/activity_metaplot.pdf",
        PROCESS_DIR + "/activity_all.pdf"


rule unibind_sites:
    message:
        """
        - Combines UniBind predictions into single sites file.
        1. awk seen removes duplicate seqs
        2. 
        """
    input:
        damos=INSTALL_DIR + "/damo_hg38_all_TFBS_unpacked"
    output:
        files=PROCESS_DIR + "/{PROFILE}/{PROFILE}-fileset.tsv",
        sites=PROCESS_DIR + "/{PROFILE}/{PROFILE}-unibind.tsv",
    params:
        profile=lambda wc: wc.PROFILE,
    shell:
        """
        set +o pipefail
        files="$(ls {input.damos}/*.bed | grep -w "{params.profile}")" &&
        for file in "${{files[@]}}"; do printf "${{file}}\n" >> {output.files}; done &&
        cat $files |
        vawk '!a[$1, $2, $3]++' | 
        #sed '/N/d' | TODO: Not sure what this does but it breaks some motifs...might be ok to just take out...
        awk -F "\t" 'BEGIN {{ OFS = FS }} {{ split ($4, motif, "_"); print $1, $2, $3, motif[2], $5, $6}}' |
        sort -k 1,1 -k2,2n > {output.sites}
        """

rule score_unibind:
    """
    - Updates UniBind sites file with motif PWM score.
    """
    input:
        sites=rules.unibind_sites.output.sites,
        pwm="results/pwmscan/{PROFILE}/{PROFILE}-pwm.txt",
    output:
        PROCESS_DIR + "/{PROFILE}/{PROFILE}-unibind_scored.tsv",
    params:
        mode="score_unibind"
    conda:
        "tfbs-activity"
    script:
        "../scripts/score.py"


rule intersect_motifs:
    """
    - Intersect UniBindm motifs with genome wide scan.
    """
    input:
        motifs="results/pwmscan/{PROFILE}/{PROFILE}-sites.bed.gz",
        unibind=rules.score_unibind.output,
    output:
        PROCESS_DIR + "/{PROFILE}/{PROFILE}-unibind_intersect.bed",
    shell:
        """
        bedtools intersect -a {input.motifs} -b {input.unibind} -c |
        vawk '{{ if ($7>=1) {{print $0, 1}} else {{print $0, 0}} }}' > {output}
        """


rule map_activity:
    """
    - Creates activity map between UniBind and reference motifs.
    - Sorts on PWM score for groupby operation.
    """
    input:
        rules.intersect_motifs.output,
    output:
        PROCESS_DIR + "/{PROFILE}/{PROFILE}-map.tsv",
    shell:
        """
        sort -k5,5n {input} |
        bedtools groupby -i stdin -g 5 -c 7  -o sum,count |
        vawk '{{ print $1, $2, $3, $2/$3 }}' > {output}
        """


rule activity_stats:
    """
    - Updates activity map with perecnetile score for each PWM score.
    - CI
    - winsorize
    """
    input:
        activity=rules.map_activity.output,
        pvals="results/pwmscan/{PROFILE}/{PROFILE}-pvals.txt",
    output:
        PROCESS_DIR + "/{PROFILE}/{PROFILE}-map_statistics.tsv",
    params:
        profile=lambda wc: wc.PROFILE,
    script:
        "../scripts/stats.py"

rule plot:
    """
    - d
    """
    input:
        activity=rules.activity_stats.output,
        nchip=rules.unibind_sites.output.files,
    output:
        PROCESS_DIR + "/{PROFILE}/{PROFILE}-activity.pdf",
    conda:
        "conda_R"
    params:
        profile=lambda wc: wc.PROFILE,
    script:
        "../scripts/plot.R"

rule combine_data:
    """
    - d
    """
    input:
        expand(PROCESS_DIR + "/{PROFILE}/{PROFILE}-map_statistics.tsv", PROFILE=PROFILES)
    output:
        PROCESS_DIR + "/activity_metadata.tsv",
    shell:
        """
        awk 'FNR==1 && NR!=1{{next;}}{{print}}' {input} > {output}
        """

rule activity_metadata:
    """
    - d
    """
    input:
        rules.combine_data.output
    output:
        PROCESS_DIR + "/activity_metadata-summary.tsv",
    script:
        "../scripts/metadata.py"

rule activity_metaplot:
    message:
        """
        d
        """
    input:
        metadata=rules.activity_metadata.output,
        combined=rules.combine_data.output,
    output:
        metaplot=PROCESS_DIR + "/activity_metaplot.pdf"
    conda:
        "conda_R"
    script:
        "../scripts/metaplot.R"

rule combine_plots:
    message:
        """
        d
        """
    input:
        expand(PROCESS_DIR + "/{PROFILE}/{PROFILE}-activity.pdf", PROFILE=PROFILES)
    output:
        all_plots=PROCESS_DIR + "/activity_all.pdf"
    conda:
        "tfbs-activity"
    script:
        "../scripts/combine.py"