from snakemake.utils import min_version


# Configuration
configfile: "config/config.yaml"


# Parameters TODO: think about how this working when doing the same in scan...
PROFILES = [i.split("|")[1] for i in config["TFBS-SCAN"]["targets"]]

# Settings
min_version("7.32.4")


# WC constraints - JASPAR matrix format
wildcard_constraints:
    PROFILE="[a-zA-Z\d]{6}.{1}\d{1}",


rule all:
    input:
        "results/activity/activity_data.combined.tsv",


rule profile_filesets:
    message:
        """
        Outputs a list of bed files covering a given TF
        """
    input:
        "resources/data/unibind/damo_hg38_all_TFBS_unpacked_flat",
    output:
        "results/activity/hg38/{PROFILE}/{PROFILE}-fileset.tsv",
    params:
        profile=lambda wc: wc.PROFILE,
    log:
        stdout="workflow/logs/profile_filesets_{PROFILE}.stdout",
        stderr="workflow/logs/uprofile_filesets_{PROFILE}.stderr",
    conda:
        "../envs/tfbs-activity.yaml"
    threads: 1
    shell:
        """
        find {input} -maxdepth 1 -type f -name "*{params.profile}*" > {output}
        """


rule unibind_sites:
    message:
        """
        Combines UniBind predictions into single sites file.
        1. vawk seen removes duplicate seqs across celltypes
        2. vawk split pulls out the matching sequence
        3. Sort by genomic coordinate
        """
    input:
        rules.profile_filesets.output,
    output:
        temp("results/activity/hg38/{PROFILE}/{PROFILE}-unibind.tsv"),
    params:
        profile=lambda wc: wc.PROFILE,
    conda:
        "../envs/tfbs-activity.yaml"
    log:
        stdout="workflow/logs/unibind_sites_{PROFILE}.stdout",
        stderr="workflow/logs/unibind_sites_{PROFILE}.stderr",
    threads: 1
    shell:
        """
        cat $(cat {input}) |
        vawk '!a[$1, $2, $3]++' | 
        vawk '{{split($4, motif, "_"); print $1, $2, $3, motif[2], $6}}' |
        sort -k 1,1 -k2,2n > {output}
        """


rule score_unibind:
    """
    Updates UniBind sites file with motif PWM score.
    """
    input:
        bed=rules.unibind_sites.output,
        pwm="results/tfbs-scan/hg38/{PROFILE}/{PROFILE}-pwm.txt",
    output:
        temp("results/activity/hg38/{PROFILE}/{PROFILE}-unibind_scored.tsv"),
    conda:
        "../envs/tfbs-activity.yaml"
    log:
        stdout="workflow/logs/score_unibind_{PROFILE}.stdout",
        stderr="workflow/logs/score_unibind_{PROFILE}.stderr",
    threads: 1
    script:
        "../scripts/pssm.py"


rule mask_unibind:
    """
    Filters UniBind sites out of custom masked regions used in scanning
    """
    input:
        sites=rules.score_unibind.output,
        exclude="results/tfbs-scan/hg38/hg38.masked_regions.bed",
    output:
        "results/activity/hg38/{PROFILE}/{PROFILE}-unibind_scored.masked.tsv",
    conda:
        "../envs/tfbs-activity.yaml"
    log:
        stdout="workflow/logs/mask_unibind_{PROFILE}.stdout",
        stderr="workflow/logs/mask_unibind_{PROFILE}.stderr",
    threads: 1
    shell:
        """
        bedtools intersect -a {input.sites} -b {input.exclude} > {output}
        """


rule intersect_motifs:
    """
    Intersect UniBind motifs with genome wide scanned motifs.
    - Intersect -c counts overlaps.
    - Vawk command makes overlaps binary (0/1), every so often more than one overlap.
    """
    input:
        motifs="results/tfbs-scan/hg38/{PROFILE}/{PROFILE}-sites.masked.bed.gz",
        unibind=rules.mask_unibind.output,
    output:
        "results/activity/hg38/{PROFILE}/{PROFILE}-unibind_intersect.bed",
    conda:
        "../envs/tfbs-activity.yaml"
    log:
        stdout="workflow/logs/intersect_motifs_{PROFILE}.stdout",
        stderr="workflow/logs/intersect_motifs_{PROFILE}.stderr",
    threads: 1
    shell:
        """
        bedtools intersect -a {input.motifs} -b {input.unibind} -c |
        vawk '{{ if ($7>=1) {{print $0, 1}} else {{print $0, 0}} }}' > {output}
        """


rule map_activity:
    """
    Creates activity map between UniBind and reference motifs.
    - Sorts on PWM score for groupby operation.
    - Note operate on 8th column which is binzrized unibind from above.
    """
    input:
        rules.intersect_motifs.output,
    output:
        temp("results/activity/hg38/{PROFILE}/{PROFILE}-map.tsv"),
    conda:
        "../envs/tfbs-activity.yaml"
    log:
        stdout="workflow/logs/map_activity_{PROFILE}.stdout",
        stderr="workflow/logs/map_activity_{PROFILE}.stderr",
    threads: 1
    shell:
        """
        sort -k5,5n {input} |
        bedtools groupby -i stdin -g 5 -c 8 -o sum,count |
        vawk '{{ print $1, $2, $3, $2/$3 }}' > {output}
        """


rule activity_stats:
    """
    Updates activity map with perecnetile score for each PWM score.
    - CI
    - winsorize
    """
    input:
        activity=rules.map_activity.output,
        pvals="results/tfbs-scan/hg38/{PROFILE}/{PROFILE}-pvals.txt",
    output:
        "results/activity/hg38/{PROFILE}/{PROFILE}-map_statistics.tsv",
    params:
        profile=lambda wc: wc.PROFILE,
    conda:
        "../envs/tfbs-activity.yaml"
    log:
        stdout="workflow/logs/activity_stats_{PROFILE}.stdout",
        stderr="workflow/logs/activity_stats_{PROFILE}.stderr",
    threads: 1
    script:
        "../scripts/stats.py"


rule combine_activity:
    """
    Combines all activity statistics into single file
    """
    input:
        expand(
            "results/activity/hg38/{PROFILE}/{PROFILE}-map_statistics.tsv",
            PROFILE=PROFILES,
        ),
    output:
        "results/activity/activity_data.combined.tsv",
    log:
        stdout="workflow/logs/combine_activity.stdout",
        stderr="workflow/logs/combine_activity.stderr",
    conda:
        "../envs/tfbs-activity.yaml"
    threads: 1
    shell:
        """
        awk 'FNR==1 && NR!=1{{next;}}{{print}}' {input} > {output}
        """
