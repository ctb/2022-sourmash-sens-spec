rule all:
    input:
        "basic-detection-100bp.csv"

rule make_curve:
    input:
        "single.fa.gz"
    output:
        "basic-detection-100bp.csv"
    shell: """
       ./scripts/make-detection-curve.py single.fa.gz \
          -o basic-detection-100bp.csv
    """

rule generate_single:
    output: "single.fa.gz"
    shell: """
       ./scripts/generate-genomes.py -n 1 -o single.fa.gz
    """
