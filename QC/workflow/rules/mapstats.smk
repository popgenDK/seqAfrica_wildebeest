rule collect_jsons:
    input:
        jsons = config["jsons"]
    output:
        directory( OUTMAIN / "jsons_statistics")
    conda:
        "../envs/python.yaml"
    log:
         OUTMAIN / "jsons_statistics.log"
    shell:
        "python3 {SCRIPT_DIR}/collect_jsons.py --jsons {input} --output_dir {output} > {log}"

rule plot_whatever:
    input:
         OUTMAIN / "jsons_statistics"
    output:
         OUTMAIN / "plots" / "matches_pct.passed.png",
         OUTMAIN / "plots" / "insert_sizes.passed.png",
         OUTMAIN / "plots" / "mapping_quality.passed.png",
         OUTMAIN / "plots" / "query_lengths.passed.png",
    conda:
        "../envs/r_stuff.yaml"
    shell: """
        Rscript {SCRIPT_DIR}/plotting/whatever.R {input}/matches_pct.passed.txt {output[0]} matches_pct bar
        Rscript {SCRIPT_DIR}/plotting/whatever.R {input}/insert_sizes.passed.txt {output[1]} insert_size line
        Rscript {SCRIPT_DIR}/plotting/whatever.R {input}/mapping_quality.passed.txt {output[2]} mapping_quality bar
        Rscript {SCRIPT_DIR}/plotting/whatever.R {input}/query_lengths.passed.txt {output[3]} query_length bar
    """

rule plot_whatever_stats_excl:
    input:
         OUTMAIN / "jsons_statistics"
    output:
         OUTMAIN / "plots" / "stats_exclusion_freq.pdf",
         OUTMAIN / "plots" / "stats_exclusion_count.pdf",
    conda:
        "../envs/r_stuff.yaml"
    shell:
        "Rscript {SCRIPT_DIR}/plotting/whatever_stats_exclusioncrit.R {input}/stats.txt {output}"

rule plot_whatever_stats_pass:
    input:
         OUTMAIN / "jsons_statistics"
    output:
         OUTMAIN / "plots" / "stats_passed_failed.png"
    conda:
        "../envs/r_stuff.yaml"
    shell:
        "Rscript {SCRIPT_DIR}/plotting/whatever_stats_passed.R {input}/stats.txt {output}"


