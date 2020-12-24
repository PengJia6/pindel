localrules: all, pindel_conf
sample_list="/home/DATA/1000G/1000G_Deep/adding_698_data/working_peng_pindel_filter/final2"
pindel="/opt/pindel/pindel"
python="/home/pengjia/miniconda3/envs/default/bin/python"
pindel_filter_script=python+" /home/DATA/1000G/1000G_Deep/adding_698_data/scripts/pindelfilter/pindel_filter.py"
samtools="/home/pengjia/miniconda3/envs/default/bin/samtools"
bcftools="/home/pengjia/miniconda3/envs/default/bin/bcftools"
bgzip="/home/pengjia/miniconda3/envs/default/bin/bgzip"
tabix="/home/pengjia/miniconda3/envs/default/bin/tabix"
ref="/home/DATA/REFGENOMEDB/human/GRCh38_full_analysis_set_plus_decoy_hla/genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
crams="/home/DATA/1000G/1000G_Deep/adding_698_data/rawCram/"
bams="/home/DATA/1000G/1000G_Deep/adding_698_data/bams/"
pindel_output="/home/DATA/1000G/1000G_Deep/adding_698_data/pindel/"
pindel_output_filter="/home/DATA/1000G/1000G_Deep/adding_698_data/pindel_filter/"
mako_output="/home/DATA/1000G/1000G_Deep/adding_698_data/mako/"
pindel_output_vcf="/home/DATA/1000G/1000G_Deep/adding_698_data/pindelvcf/"
pindel_filter_vcf="/home/DATA/1000G/1000G_Deep/adding_698_data/pindelvcf_filter/"
logs="/home/DATA/1000G/1000G_Deep/adding_698_data/working_peng/logs/"
cram2bam_threads=8
bamindex_threads=8
pindel_threads=2
pindel2vcf_threads=1
chroms=["chr"+str(i) for i in range(1,23)]+["chrX","chrY","chrM"]
samples=[]
for line in open(sample_list):
    if len(line)>2:
        samples.append(line[:-1].rstrip(".final.cram").rstrip(".final.bam"))

#samples=["NA18913"]
wildcard_constraints:
    sample="|".join(samples),
    chrom="|".join(chroms)

rule all: 
    input: 
#      expand(bams+"{sample}.final.bam.bai",sample=samples),
      expand(pindel_filter_vcf+"{sample}.pindel.illumina_high_coverage.20201009.sites.vcf.gz.tbi",sample=samples),
#      expand(pindel_output+"{sample}.final/{sample}.final_{chrom}.vcf",sample=samples,chrom=chroms),
#   expand(mako_output+"{sample}.mako.cfg",sample=samples),
#     expand(mako_output+"{sample}.mako.finished.tag",sample=samples),
rule cram2bam:
    input: 
      cram=crams+"{sample}.final.cram",
      ref=ref
    output:
      bam=bams+"{sample}.final.bam"
    threads:
      cram2bam_threads
    log:
      logs+"cram2bam/{sample}.log"
    benchmark:
      logs+"cram2bam/{sample}.bm"
    shell:
      samtools+" view -T {input.ref} -b -@ {threads} -o {output.bam} {input.cram} 2>{log} 1>{log}"
rule bamindex:
    input: 
      bam=bams+"{sample}.final.bam"
    output:
      bam_index=bams+"{sample}.final.bam.bai"
    threads:
      bamindex_threads
    log:
      logs+"bamindex/{sample}.log"
    benchmark:
      logs+"bamindex/{sample}.bm"
    shell: 
      samtools+" index -@ {threads} {input.bam} 2>{log} 1>{log}"

rule pindel_conf:
    input:
      bam=bams+"{sample}.final.bam",
      bai=bams+"{sample}.final.bam.bai",
    output: 
      conf=pindel_output+"{sample}.final/{sample}.final.conf" 
    shell: 
      "echo {input.bam} 450 {wildcards.sample} > {output.conf}" 

rule pindel_chrom: 
    input: 
      conf=rules.pindel_conf.output.conf,
      ref=ref
    output: 
      out=pindel_output+"{sample}.final/{sample}.final_{chrom}" 
    threads:
      pindel_threads
    log:
      logs+"pindel/{sample}_{chrom}.log"
    benchmark:
      logs+"pindel/{sample}_{chrom}.bm"
    run: 
      shell(pindel+" -f {input.ref} -i {input.conf} -T {threads} -c {wildcards.chrom} -x 4 -M 3 -B 100000 -A 20 -o {output.out}")
      shell("touch {output.out}")
      
rule pindel_filter:
    input:
      out=rules.pindel_chrom.output.out,
    output:
      out=pindel_output_filter+"{sample}.final/{sample}.final_{chrom}"
    run:
        shell(pindel_filter_script+" -i {input.out} -o {output.out} ")
        shell("touch {output.out}")

rule pindel2vcf: 
    input: 
      out=rules.pindel_filter.output.out,
      ref=ref
    output: 
      out=pindel_output_filter+"{sample}.final/{sample}.final_{chrom}.vcf" 
    threads:
      pindel2vcf_threads
    log:
      logs+"pindel2vcf2/{sample}_{chrom}.log"
    benchmark:
      logs+"pindel2vcf2/{sample}_{chrom}.bm"
    run: 
      shell("/opt/pindel/pindel2vcf -P {input.out} -r {input.ref} -c {wildcards.chrom} -d 20200828 -v {output} -R  GRCh38 2>{log} 1>{log}")

rule pindel_vcf_compress: 
    input: 
      vcf=rules.pindel2vcf.output.out,
    output: 
      out=pindel_output_filter+"{sample}.final/{sample}.final_{chrom}.vcf.gz" 
    threads:
      2
    log:
      logs+"pindel_vcf_compress2/{sample}_{chrom}.log"
    benchmark:
      logs+"pindel_vcf_compress2/{sample}_{chrom}.bm"
    run: 
      shell(bgzip+" -c  {input.vcf} >{output.out}")


rule pindel_vcf_index: 
    input: 
      vcf="{prefix}.vcf.gz",
    output: 
      out="{prefix}.vcf.gz.tbi" 
    threads:
      1
    run: 
      shell(tabix+" {input.vcf} ")

rule pindel_vcf_merge: 
    input: 
      vcfs=expand(pindel_output_filter+"{{sample}}.final/{{sample}}.final_{chrom}.vcf.gz",chrom=chroms),
      index=expand(pindel_output_filter+"{{sample}}.final/{{sample}}.final_{chrom}.vcf.gz.tbi",chrom=chroms),
    output: 
      out=pindel_filter_vcf+"{sample}.pindel.illumina_high_coverage.20201009.sites.vcf.gz" 
    threads:
      1
    log:
      logs+"pindel_vcf_merge/{sample}.log"
    benchmark:
      logs+"pindel_vcf_merge/{sample}.bm"
    run: 
      inputs= " ".join(["{}".format(f) for f in input.vcfs])
      shell(bcftools+" concat -Oz -o {output.out}  {inputs} 2>{log} 1>{log} ")


# $PINDEL2VCF -P $OUTDIR/$SAMPLE.pindel -r $REF -R GRCH38 -d 20160321
# $PINDEL -f $REF -i $OUTDIR/.$SAMPLE.pindel_config -T $NO_THREAD -x 4 -M 3 -B 100000 -A 20 -o $OUTDIR/$SAMPLE.pindel
