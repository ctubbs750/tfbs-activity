from os import listdir, path
from snakemake.utils import min_version

# Settings
min_version("7.32.4")


# ------------- #
# Config        #
# ------------- #

GENOME_SCAN_DIR = config["genome_scan_dir"]
UNIBIND_SITES_DIR = config["unibind_sites_dir"]
OUTPUT_DIR = config["output_dir"]


# ------------- #
# I/O           #
# ------------- #

# Input Unibind hits file
UNIBIND_SITES = os.path.join(
    UNIBIND_SITES_DIR, "{tf_name}", "{profile}", "{dataset}.bed"
)

# Input genome-wide sites file
GENOME_SITES = os.path.join(
    GENOME_SCAN_DIR,
    "scan",
    "{tf_name}",
    "{profile}",
    "{dataset}",
    "sites.masked.genome.sorted.bed.starch",
)


# Gencode protein-coding exons
MASK_REGIONS = os.path.join(GENOME_SCAN_DIR, "mask", "hg38.masked_regions.bed")

# Masked unibind sites
UNIBIND_SITES_MASKED = os.path.join(
    OUTPUT_DIR,
    "indiv",
    "{tf_name}",
    "{profile}",
    "datasets",
    "{dataset}",
    "sites.masked.bed",
)

# Intersection between genome and unibind sites
SITES_IX = os.path.join(
    OUTPUT_DIR,
    "indiv",
    "{tf_name}",
    "{profile}",
    "datasets",
    "{dataset}",
    "unibind_intersect.bed",
)

# Non annotated activity map
ACTIVITY_MAP_RAW = os.path.join(
    OUTPUT_DIR,
    "indiv",
    "{tf_name}",
    "{profile}",
    "datasets",
    "{dataset}",
    "activity-map.raw.tsv",
)

# Pval dist from matrix_prob
PWM_PVALS = os.path.join(
    GENOME_SCAN_DIR, "scan", "{tf_name}", "{profile}", "{dataset}", "pvals.txt"
)


# AUC for binding
PROFILE_AUC = os.path.join(
    OUTPUT_DIR,
    "indiv",
    "{tf_name}",
    "{profile}",
    "datasets",
    "{dataset}",
    "profile_auc.tsv",
)

# Activity Error
ACTIVITY_CI = os.path.join(
    OUTPUT_DIR,
    "indiv",
    "{tf_name}",
    "{profile}",
    "datasets",
    "{dataset}",
    "activity.error.tsv",
)

# Final activity map with stats
ACTIVITY_MAP = os.path.join(
    OUTPUT_DIR,
    "indiv",
    "{tf_name}",
    "{profile}",
    "datasets",
    "{dataset}",
    "activity-map.final.tsv",
)

# Activity plot
ACTIVITY_PLOT = os.path.join(
    OUTPUT_DIR,
    "indiv",
    "{tf_name}",
    "{profile}",
    "datasets",
    "{dataset}",
    "activity.map.png",
)

# Meta-analyized activity over celltypes
# META_ACTIVITY = os.path.join(
#     OUTPUT_DIR,
#     "indiv",
#     "{tf_name}",
#     "{profile}",
#     "summary",
#     "{profile}-meta_activity.map.tsv",
# )

# Combined AUC
COMBINED_AUC = os.path.join(OUTPUT_DIR, "summary", "summary_aucs.tsv")

# ------------- #
# Params        #
# ------------- #

TF_NAMES, PROFILES, DATASETS = glob_wildcards(
    os.path.join(UNIBIND_SITES_DIR, "{tf_name}", "{profile}", "{dataset}.bed")
)

PROFILE_TO_DATASETS = {}
for profile in PROFILES:
    PROFILE_TO_DATASETS[profile] = [ds for ds in DATASETS if profile in ds]

# ------------- #
# Rules         #
# ------------- #


rule all:
    input:
        expand(
            PROFILE_AUC,
            zip,
            tf_name=TF_NAMES,
            profile=PROFILES,
            dataset=DATASETS,
        ),
        # expand(
        #     ACTIVITY_CI,
        #     zip,
        #     tf_name=TF_NAMES,
        #     profile=PROFILES,
        #     dataset=DATASETS,
        # ),
        # expand(
        #     ACTIVITY_PLOT,
        #     zip,
        #     tf_name=TF_NAMES,
        #     profile=PROFILES,
        #     dataset=DATASETS,
        # ),
        # expand(
        #     "results/activity/indiv/{tf_name}/{profile}/summary/{profile}-meta_activity.map.tsv",
        #     zip,
        #     tf_name=TF_NAMES,
        #     profile=PROFILES,
        # ),
        COMBINED_AUC,
    default_target: True


rule mask_unibind:
    """
    Filters UniBind sites out of custom masked regions used in scanning
    """
    input:
        unibind=UNIBIND_SITES,
        exclude=MASK_REGIONS,
    output:
        temp(UNIBIND_SITES_MASKED),
    conda:
        "../envs/tfbs-activity.yaml"
    log:
        stdout="workflow/logs/mask_unibind_{tf_name}_{profile}_{dataset}.stdout",
        stderr="workflow/logs/mask_unibind_{tf_name}_{profile}_{dataset}.stderr",
    shell:
        "bedtools intersect -a {input.unibind} -b {input.exclude} -v > {output}"


rule intersect_motifs:
    """
    Intersect UniBind motifs with genome wide scanned motifs.
    - Intersect -c counts overlaps.
    - Vawk command makes overlaps binary (0/1), every so often more than one overlap.
    """
    input:
        genome_sites=GENOME_SITES,
        unibind_sites=rules.mask_unibind.output,
    output:
        temp(SITES_IX),
    conda:
        "../envs/tfbs-activity.yaml"
    log:
        stdout="workflow/logs/intersect_motifs_{tf_name}_{profile}_{dataset}.stdout",
        stderr="workflow/logs/intersect_motifs_{tf_name}_{profile}_{dataset}.stderr",
    threads: 1
    shell:
        """
        unstarch {input.genome_sites} |
        bedtools intersect -a stdin -b {input.unibind_sites} -c |
        vawk '{{if ($7>0) {{print $1, $2, $3, $4, $5, $6, 1}} else {{print $1, $2, $3, $4, $5, $6, 0}} }}' > {output}
        """


rule map_activity:
    """
    Creates activity map between UniBind and reference motifs.
    - Sorts on PWM score for groupby operation.
    - Note operate on 7th column which is binzrized unibind from above.
    """
    input:
        rules.intersect_motifs.output,
    output:
        temp(ACTIVITY_MAP_RAW),
    conda:
        "../envs/tfbs-activity.yaml"
    log:
        stdout="workflow/logs/map_activity_{tf_name}_{profile}_{dataset}.stdout",
        stderr="workflow/logs/map_activity_{tf_name}_{profile}_{dataset}.stderr",
    shell:
        """
        sort -k5,5n {input} |
        bedtools groupby -i stdin -g 5 -c 7 -o sum,count |
        vawk '{{ print $1, $2, $3, $2/$3 }}' > {output}
        """


rule calculate_auc:
    """
    Calculate AUC of binding in unibind
    Thread count is to help memory issues...
    """
    input:
        unibind_ix=rules.intersect_motifs.output,
        profile_pvals=PWM_PVALS,
    output:
        PROFILE_AUC,
    params:
        profile=lambda wc: wc.profile,
        dataset=lambda wc: wc.dataset,
    conda:
        "../envs/tfbs-activity.yaml"
    log:
        stdout="workflow/logs/map_activity_{tf_name}_{profile}_{dataset}.stdout",
        stderr="workflow/logs/map_activity_{tf_name}_{profile}_{dataset}.stderr",
    threads: 12
    script:
        "../scripts/stat/auc.py"


rule calculate_error:
    """
    Calculate proportion error on activity scores
    """
    input:
        activity=rules.map_activity.output,
        pvals=PWM_PVALS,
    output:
        ACTIVITY_CI,
    params:
        profile=lambda wc: wc.profile,
        dataset=lambda wc: wc.dataset,
    conda:
        "../envs/tfbs-activity.yaml"
    log:
        stdout="workflow/logs/calculate_error_{tf_name}_{profile}_{dataset}.stdout",
        stderr="workflow/logs/calculate_error_{tf_name}_{profile}_{dataset}.stderr",
    script:
        "../scripts/stat/error.py"


rule confident_activity:
    """
    Adjusts activity scores based off of...
    """
    input:
        activity=rules.map_activity.output,
    output:
        ACTIVITY_MAP,
    params:
        window=50,
        threshold=0.10,
    conda:
        "../envs/tfbs-activity.yaml"
    log:
        stdout="workflow/logs/winsorize_activity_{tf_name}_{profile}_{dataset}.stdout",
        stderr="workflow/logs/winsorize_activity_{tf_name}_{profile}_{dataset}.stderr",
    threads: 1
    script:
        "../scripts/stat/adjust.py"


rule plot_activity:
    """
    Makes indvidual data plot based off activity map
    TODO: Update to plot off of winsorized activity scores
    """
    input:
        activity=rules.confident_activity.output,
    output:
        ACTIVITY_PLOT,
    params:
        profile=lambda wc: wc.profile,
        dataset=lambda wc: wc.dataset,
    conda:
        "../envs/report.yaml"
    log:
        stdout="workflow/logs/plot_activity_{tf_name}_{profile}_{dataset}.stdout",
        stderr="workflow/logs/plot_activity_{tf_name}_{profile}_{dataset}.stderr",
    script:
        "../scripts/plot/plot.R"


# rule meta_activity:
#     """
#     Combines each TF-profile cell specific activity maps to meta per profile.
#     """
#     input:
#         lambda wildcards: expand(
#             "results/activity/indiv/{{tf_name}}/{{profile}}/datasets/{dataset}/activity.error.tsv",
#             dataset=PROFILE_TO_DATASETS[wildcards.profile],
#         ),
#     output:
#         "results/activity/indiv/{tf_name}/{profile}/summary/{profile}-meta_activity.map.tsv",
#     conda:
#         "../envs/tfbs-activity.yaml"
#     log:
#         stdout="workflow/logs/meta_activity_{tf_name}_{profile}.stdout",
#         stderr="workflow/logs/meta_activity_{tf_name}_{profile}.stderr",
#     threads: 1
#     shell:
#         "echo {input} > {output}"


rule combine_aucs:
    """
    Combines indiv AUC estimation for each dataset into single file per TF-profile
    """
    input:
        expand(
            PROFILE_AUC,
            zip,
            tf_name=TF_NAMES,
            profile=PROFILES,
            dataset=DATASETS,
        ),
    output:
        COMBINED_AUC,
    conda:
        "../envs/tfbs-activity.yaml"
    log:
        stdout="workflow/logs/combine_aucs.stdout",
        stderr="workflow/logs/combine_aucs.stderr",
    shell:
        "cat {input} > {output}"
