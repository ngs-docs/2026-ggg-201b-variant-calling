rule uncompress_genome:
    input: "ecoli-rel606.fa.gz"
    output: "outputs/ecoli-rel606.fa"
    shell: """
        gunzip -c ecoli-rel606.fa.gz > outputs/ecoli-rel606.fa
    """

rule map_reads:
    input: "outputs/ecoli-rel606.fa", "SRR2584857_1.fastq.gz"
    output: "outputs/SRR2584857_1.x.ecoli-rel606.sam"
    conda: "mapping"
    shell: """
        minimap2 -ax sr outputs/ecoli-rel606.fa SRR2584857_1.fastq.gz > \
            outputs/SRR2584857_1.x.ecoli-rel606.sam
    """

rule sam_to_bam:
    input: "outputs/SRR2584857_1.x.ecoli-rel606.sam"
    output: "outputs/SRR2584857_1.x.ecoli-rel606.bam"
    conda: "mapping"
    shell: """
        samtools view -b outputs/SRR2584857_1.x.ecoli-rel606.sam > \
            outputs/SRR2584857_1.x.ecoli-rel606.bam
     """

rule sort_bam:
    input: "outputs/SRR2584857_1.x.ecoli-rel606.bam"
    output: "outputs/SRR2584857_1.x.ecoli-rel606.sorted.bam"
    conda: "mapping"
    shell: """
        samtools sort outputs/SRR2584857_1.x.ecoli-rel606.bam > \
            outputs/SRR2584857_1.x.ecoli-rel606.sorted.bam
    """

rule index_bam:
    input: "outputs/SRR2584857_1.x.ecoli-rel606.sorted.bam"
    output: "outputs/SRR2584857_1.x.ecoli-rel606.sorted.bam.bai"
    conda: "mapping"
    shell: """
        samtools index outputs/SRR2584857_1.x.ecoli-rel606.sorted.bam
    """

rule call_variants:
    input: "outputs/ecoli-rel606.fa", "outputs/SRR2584857_1.x.ecoli-rel606.sorted.bam"
    output: "outputs/SRR2584857_1.x.ecoli-rel606.pileup", "outputs/SRR2584857_1.x.ecoli-rel606.bcf", "outputs/SRR2584857_1.x.ecoli-rel606.vcf"
    conda: "mapping"
    shell: """
        bcftools mpileup -Ou -f outputs/ecoli-rel606.fa \
            outputs/SRR2584857_1.x.ecoli-rel606.sorted.bam  > \
            outputs/SRR2584857_1.x.ecoli-rel606.pileup
        bcftools call -mv -Ob outputs/SRR2584857_1.x.ecoli-rel606.pileup \
            -o outputs/SRR2584857_1.x.ecoli-rel606.bcf
        bcftools view outputs/SRR2584857_1.x.ecoli-rel606.bcf -o \
            outputs/SRR2584857_1.x.ecoli-rel606.vcf
    """
