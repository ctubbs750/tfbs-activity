from snakemake.utils import min_version


# Configuration
configfile: "config/config.yaml"


# Parameters TODO: think about how this working when doing the same in scan...
UNIBIND_URL = config["TFBS-ACTIVITY"]["unibind_url"]

# Settings
min_version("7.32.4")


rule all:
    input:
        "resources/data/unibind/targets.psv",
        "resources/data/unibind/damo_hg38_all_TFBS_unpacked_flat",


rule download_unibind:
    message:
        """
        Downloads all DAMOS from UniBind database
        """
    output:
        temp("resources/data/unibind/damo_hg39_all_TFBS.tar.gz"),
    params:
        url=UNIBIND_URL,
    log:
        stdout="workflow/logs/download_unibind.stdout",
        stderr="workflow/logs/download_unibind.stderr",
    conda:
        "../envs/tfbs-activity.yaml"
    threads: 1
    shell:
        """
        curl -o {output} {params.url}
        """


rule unpack_unibind:
    message:
        """
        Unpacks UniBind data download
        """
    input:
        rules.download_unibind.output,
    output:
        temp(directory("resources/data/unibind/damo_hg38_all_TFBS_unpacked")),
    log:
        stdout="workflow/logs/unpack_unibind.stdout",
        stderr="workflow/logs/unpack_unibind.stderr",
    conda:
        "../envs/tfbs-activity.yaml"
    threads: 1
    shell:
        """
        mkdir -p {output} && tar -xzf {input} -C {output}
        """


rule flatten_dir:
    message:
        """
        Flattens unpacked tarbell. Not great - for the moment just check for creation of first in OP.
        Also removes all the empty dirs.
        """
    input:
        rules.unpack_unibind.output,
    output:
        directory("resources/data/unibind/damo_hg38_all_TFBS_unpacked_flat"),
    log:
        stdout="workflow/logs/flatten_dir.stdout",
        stderr="workflow/logs/flatten_dir.stderr",
    conda:
        "../envs/tfbs-activity.yaml"
    threads: 1
    shell:
        """
        mkdir {output} &&
        find {input} -mindepth 2 -type f -exec mv -t {output} -i '{{}}' + &&
        find {input} -type d -empty -delete
        """


rule fetch_targets:
    message:
        """
        Returns list of all TFs and their JASPAR profiles from UniBind
        """
    input:
        rules.flatten_dir.output,
    output:
        temp("resources/data/unibind/make_targets.txt"),
    log:
        stdout="workflow/logs/make_targets.stdout",
        stderr="workflow/logs/make_targets.stderr",
    conda:
        "../envs/tfbs-activity.yaml"
    threads: 1
    shell:
        """
        ls {input} | cut -d"." -f3,4-5 | sort -u > {output}
        """


rule failed_targets:
    message:
        """
        Test each matrix ID for membership in JASPAR database.
        """
    input:
        rules.fetch_targets.output,
    output:
        temp("resources/data/unibind/failed_targets.txt"),
    conda:
        "../envs/tfbs-activity.yaml"
    log:
        stdout="workflow/logs/failed_targets.stdout",
        stderr="workflow/logs/failed_targets.stderr",
    threads: 1
    shell:
        """
        coreapi get https://jaspar.elixir.no/api/v1/docs/
        readarray -t matrices < <(cut -d"." -f2-3 {input})
        for matrix in "${{matrices[@]}}";
        do
            coreapi action matrix read -p matrix_id=$matrix || echo $matrix >> {output}
        done
        """


rule filter_targets:
    message:
        """
        Filter profiles not in JASPAR.
        """
    input:
        target=rules.fetch_targets.output,
        failed=rules.failed_targets.output,
    output:
        temp("resources/data/unibind/filter_targets.txt"),
    log:
        stdout="workflow/logs/filter_targets.stdout",
        stderr="workflow/logs/filter_targets.stderr",
    conda:
        "../envs/tfbs-activity.yaml"
    threads: 1
    shell:
        """
        cat {input.target} | grep -v -f {input.failed} > {output}
        """


rule fetch_lengths:
    message:
        """
        Returns list of all matrix lengths from JASPAR API
        """
    input:
        rules.filter_targets.output,
    output:
        temp("resources/data/unibind/target_lengths.txt"),
    conda:
        "../envs/tfbs-activity.yaml"
    log:
        stdout="workflow/logs/fetch_lengths.stdout",
        stderr="workflow/logs/fetch_lengths.stderr",
    threads: 1
    shell:
        """
        coreapi get https://jaspar.elixir.no/api/v1/docs/
        readarray -t matrices < <(cut -d"." -f2-3 {input})
        for matrix in "${{matrices[@]}}";
        do
            echo $matrix
            coreapi action matrix read -p matrix_id=$matrix | jq '.pfm.A[]' | wc -l >> {output}
        done
        """


rule combine_targets:
    message:
        """
        Combines targets and their profile lengths.
        """
    input:
        targets=rules.filter_targets.output,
        lengths=rules.fetch_lengths.output,
    output:
        temp("resources/data/unibind/combine_targets.txt"),
    log:
        stdout="workflow/logs/combine_targets.stdout",
        stderr="workflow/logs/combine_targets.stderr",
    conda:
        "../envs/tfbs-activity.yaml"
    threads: 1
    shell:
        """
        paste -d"." {input.targets} {input.lengths} > {output}
        """


rule format_targets:
    message:
        """
        Makes final target list.
        """
    input:
        rules.combine_targets.output,
    output:
        "resources/data/unibind/targets.psv",
    log:
        stdout="workflow/logs/format_targets.stdout",
        stderr="workflow/logs/format_targets.stderr",
    conda:
        "../envs/tfbs-activity.yaml"
    threads: 1
    shell:
        """
        cat {input} |
        cut -d"." -f1-4 --output-delimiter=$'\t' |
        vawk '{{print $1"|"$2"."$3"|"$4}}' > {output}
        """
