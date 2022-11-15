rule all:
    input:
        "basic-detection-100bp.csv",
        "basic-detection-10000bp.csv",

rule make_curve_wc:
    input:
        genome="single.fa.gz"
    output:
        csv="basic-detection-{bp}bp.csv"
    shell: """
       ./scripts/make-detection-curve.py single.fa.gz \
          -o {output.csv} -r {wildcards.bp}
    """

rule generate_single:
    output: "single.fa.gz"
    shell: """
       ./scripts/generate-genomes.py -n 1 -o single.fa.gz
    """
