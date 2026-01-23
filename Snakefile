SAMPLES = ["SRR2584857_1"]
GENOME = ["ecoli-rel606"]

rule make_vcf:
    input:
        expand("outputs/{sample}.x.{genome}.vcf",
               sample=SAMPLES, genome=GENOME),

rule uncompress_genome:
    input: "{xyz}.fa.gz"
    output: "outputs/{xyz}.fa"
    shell: """
        gunzip -c {input} > {output}
    """

rule map_reads:
    input:
        genome="outputs/{genome}.fa",
        sample="{sample}.fastq.gz",
    output: "outputs/{sample}.x.{genome}.sam"
    conda: "mapping"
    shell: """
        minimap2 -ax sr {input.genome} {input.sample} > {output}
    """

rule sam_to_bam:
    input: "outputs/{sample}.x.{genome}.sam"
    output: "outputs/{sample}.x.{genome}.bam"
    conda: "mapping"
    shell: """
        samtools view -b {input} > {output}
     """

rule sort_bam:
    input: "outputs/{sample}.x.{genome}.bam"
    output: "outputs/{sample}.x.{genome}.sorted.bam"
    conda: "mapping"
    shell: """
        samtools sort {input} > {output}
    """

rule index_bam:
    input: "outputs/{sample}.x.{genome}.sorted.bam"
    output: "outputs/{sample}.x.{genome}.sorted.bam.bai"
    conda: "mapping"
    shell: """
        samtools index {input}
    """

rule call_variants:
    input:
        genome="outputs/{genome}.fa",
        bam="outputs/{sample}.x.{genome}.sorted.bam",
        bai="outputs/{sample}.x.{genome}.sorted.bam.bai",
    output:
        pileup="outputs/{sample}.x.{genome}.pileup",
        bcf="outputs/{sample}.x.{genome}.bcf",
        vcf="outputs/{sample}.x.{genome}.vcf"
    conda: "mapping"
    shell: """
        bcftools mpileup -Ou -f {input.genome} {input.bam}> {output.pileup}
        bcftools call -mv -Ob {output.pileup} -o {output.bcf}
        bcftools view {output.bcf} -o {output.vcf}
    """
