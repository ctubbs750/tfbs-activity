# Parameters
INSTALL_DIR = config["TFBS-ACTIVITY"]["INSTALL_DIR"]
PROCESS_DIR = config["TFBS-ACTIVITY"]["PROCESS_DIR"]
PROFILES = config["PWM-SCAN"]["TARGETS"]


# WC constraints - JASPAR matrix format
wildcard_constraints:
    PROFILE="[a-zA-Z\d]{6}.{1}\d{1}",


rule all:
    input:
        expand(PROCESS_DIR + "/{PROFILE}/{PROFILE}-unibind.tsv", PROFILE=PROFILES),


rule unibind_sites:
    message:
        """
        - Combines UniBind predictions into single sites file.
        """
    output:
        PROCESS_DIR + "/{PROFILE}/{PROFILE}-unibind.tsv",
    params:
        damodir=INSTALL_DIR + "/damo_hg38_all_TFBS",
        profile=lambda wc: wc.PROFILE,
    shell:
        """
        set +o pipefail
        damos="$(ls {params.damodir}/*.bed | grep -w {params.profile})" &&
        cat $damos |
        vawk '!a[$1, $2, $3]++' | 
        sed '/N/d' |
        awk -F "\t" 'BEGIN {{ OFS = FS }} {{ split ($4, motif, "_"); print $1, $2, $3, motif[2], $5, $6}}' |
        sort -k 1,1 -k2,2n > {output}
        """


# rule score_unibind:
#     """
#     - Updates UniBind sites file with motif PWM score.
#     """
#     input:
#         sites=rules.unibind_sites.output,
#         pwm=rules.calculate_pwm.output,
#     output:
#         f"{JASPAR}-unibind_scored.tsv",
#     run:
#         # Read in sites data and pwm
#         sites = utilities.read_damo(input[0])
#         pssm = utilities.setup_pssm(input[1])
#         # Score motif sequence and write out
#         sites["s"] = sites.apply(
#             lambda row: utilities.score_sequence(row.seq, pssm), axis=1
#         )
#         sites.to_csv(output[0], sep="\t", index=False, header=False)


# rule intersect_motifs:
#     """
#     - Intersect UniBindm motifs with genome wide scan.
#     """
#     input:
#         motifs=rules.scan_genome.output,
#         unibind=rules.score_unibind.output,
#     output:
#         f"{JASPAR}-unibind_intersect.bed",
#     shell:
#         """
#         bedtools intersect -a {input.motifs} -b {input.unibind} -c |
#         vawk '{{ if ($7>=1) {{print $0, 1}} else {{print $0, 0}} }}' > {output}
#         """


# rule map_acitvity:
#     """
#     - Creates activity map between UniBind and reference motifs.
#     - Sorts on PWM score for groupby operation.
#     """
#     input:
#         rules.intersect_motifs.output,
#     output:
#         temp(f"{JASPAR}-map.tsv"),
#     shell:
#         """
#         sort -k5,5n {input} |
#         bedtools groupby -i stdin -g 5 -c 7  -o sum,count |
#         vawk '{{ print $1, $2, $3, $2/$3 }}' > {output}
#         """


# rule merge_percentiles:
#     """
#     - Updates activity map with perecnetile score for each PWM score.
#     """
#     input:
#         activity=rules.map_acitvity.output,
#         pvals=rules.calculate_pvals.output,
#     output:
#         temp(f"{JASPAR}-map_percentiles.tsv"),
#     shell:
#         """
#         join --check-order -t $'\t' -j1 <(sort -k1n {input.activity}) <(sort -k1n {input.pvals})
#         """

# rule activity_ci:
#     input:
#         rules.merge_percentiles.output,
#     output:
#         temp(f"{JASPAR}-propse.tsv"),
#     run:
#         # Read in
#         activity = read_csv(input[0], sep="\t", header=None)
#         # Standard error on proportion
#         activity["ci95_lbound"], activity["ci95_rbound"] = activity.proportion_ci(
#             activity["unibind_count"], activity["motif_count"]
#         )
#         # Write final
#         activity.to_csv(output[0], sep="\t", index=False)
# rule winsorize_activity:
#     input:
#         rules.proportion_se.output,
#     output:
#         f"{JASPAR}-activity.tsv",
#     run:
#         # Read activity scores
#         data = read_csv(input[0], sep="\t", header=None)
#         # Winsorize activity and SE, 0% on left tail, 35% of values on the right tale which is about 90percentile PWM score
#         data["activity_winsor"] = activity.winsorize_array(
#             data["activity"], limits=(0.00, 0.35)
#         )
#         data["ci95_lbound_winsor"] = activity.winsorize_array(
#             data["ci95_lbound"], limits=(0.00, 0.35)
#         )
#         data["ci95_rbound_winsor"] = activity.winsorize_array(
#             data["ci95_rbound"], limits=(0.00, 0.35)
#         )
#         # Write out
#         data.to_csv(output[0], index=False, sep="\t")
