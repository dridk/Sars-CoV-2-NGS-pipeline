configfile: "config.yml"
report: "report/workflow.rst"

from glob import glob 
import re
import os


# Global variable 
GENOM = config["GENOM"]
GENOM_BWA_INDEX = GENOM + ".bwt"
GENOM_FAI_INDEX = GENOM + ".fai"
GENOM_NAME="NC_045512.2"

# Get samples names from fastq 
FASTQ_DIR= config["FASTQ_DIR"]
SAMPLES = [re.findall(f"{FASTQ_DIR}/([^_]+)_1.fastq", i)[0] for i in glob(os.path.join(FASTQ_DIR,"/*_1.fastq.gz"))]


print(SAMPLES)

rule index_genom:
	input:
		GENOM
	output:
		GENOM_FAI_INDEX,
		GENOM_BWA_INDEX
	shell:
		"bwa index {input}; samtools faidx {input}"

rule align:
	input:
		R1=FASTQ_DIR + "/{sample}_1.fastq.gz",
		R2=FASTQ_DIR + "/{sample}_2.fastq.gz",
		index=GENOM_BWA_INDEX

	output:
		temporary("{sample}.sam")
	log:
		"{sample}.align.log"
	params:
		genom = GENOM
	shell:
		"bwa mem -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}' {params.genom} {input.R1} {input.R2} > {output} 2> {log}"


rule sam2bam:
	input:
		"{sample}.sam"
	output:
		"{sample}.bam"
	shell:
		"samtools sort -O BAM {input} > {output}"

rule samIndex:
	input:
		"{sample}.bam"
	output:
		"{sample}.bam.bai"
	shell:
		"samtools index {input}"


rule vcf:
	input:
		"{sample}.bam", "{sample}.bam.bai", GENOM_FAI_INDEX
	output:
		"{sample}.vcf"
	params:
		genom = GENOM
	shell:
		"freebayes -f {GENOM} -p 1 -C10 {input[0]} > {output} "


rule bgzip:
	input:
		"{sample}.vcf"
	output:
		"{sample}.vcf.gz"
	shell:
		"bgzip {input} ;  tabix -p vcf {output}"


rule merge_vcf:
	input: 
		[f'{sample}.vcf.gz' for sample in SAMPLES]
	output:
		"final.vcf.gz"
	shell:
		"bcftools merge {input} -Oz -o {output}"

rule merge_annotation:
	input:
		"final.vcf.gz"
	output:
		"final.ann.vcf"
	log:
		"final.snpeff.log"
	shell:
		"snpEff -Xmx10G -v {GENOM_NAME} {input}> {output} 2> {log}"

rule single_annotation:
	input:
		"{sample}.vcf.gz"
	output:
		"{sample}.ann.vcf"
	log:
		"{sample}.snpeff.log"
	shell:
		"snpEff -Xmx10G -v {GENOM_NAME} {input}> {output} 2> {log}"


rule ann_to_csv:
	input:
		"{sample}.ann.vcf"
	output:
		"{sample}.results.csv"
	params:
		qual = config["VCF_QUAL_FILTER"]
	shell:
		"""
		SnpSift filter 'QUAL > {params.qual}' {input} |
		./scripts/vcfEffOnePerLine.pl|
		SnpSift extractFields - 'ANN[*].GENE' 'ANN[*].FEATUREID' 'POS' 'REF' 'ALT' 'ANN[*].HGVS_C' 'ANN[*].HGVS_P' 'ANN[*].IMPACT' 'ANN[*].EFFECT' > {output}
		"""



rule consensus:
	input:
		"{sample}.vcf.gz"
	output:
		report("{sample}.fa",caption="report/sample.rst",category="sample")
	params:
		genom = GENOM
	log:
		"{sample}.consensus.log"
	shell:
		"bcftools consensus {input} -f {params.genom} --sample {wildcards.sample} > {output} 2> {log}"





# rule pangolin:
# 	input:
# 		"{sample}.fa"
# 	output:
# 		"{sample}.lineage_report.csv"
# 	conda:
# 		"pangolin"
# 	shell:
# 		"pangolin {input} --outfile {output}"