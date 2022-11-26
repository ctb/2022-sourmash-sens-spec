rule all:
    input:
        "basic-detection-100bp.csv",
        "basic-detection-10000bp.csv",
        expand("genomes-100k.sketch.k{k}.sqldb", k=range(9, 21, 2)),
        "gtdb.family_confused.sig.gz",

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
       ./scripts/generate-genomes.py -n 1 -o {output}
    """

rule genomes_100k:
    output: protected("genomes-100k.fa.gz"),
    shell: """
       ./scripts/generate-genomes.py -n 100000 -o {output}
    """

rule sketch_100k:
    input: "genomes-100k.fa.gz",
    output: protected("genomes-100k.sketch.k{ksize}.zip"),
    shell: """
       sourmash sketch dna {input} -o {output} -p k={wildcards.ksize} \
          --singleton
    """

rule make_sqldb:
    input: "genomes-100k.{x}.zip",
    output: protected("genomes-100k.{x}.sqldb"),
    shell: """
       sourmash sig cat {input} -o {output}
    """

rule count_coherent:
    input: "gtdb-rs207.genomic.k31.sqldb",
    output: "gtdb.family_confused.sig.gz"
    shell: """
       scripts/count-coherent-hashvals.py {input} -o {output}
    """
