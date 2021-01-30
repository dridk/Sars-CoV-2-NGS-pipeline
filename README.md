# Sars-CoV-2-NGS-pipeline
A simple snakemake pipeline to call variant from NGS data of Sars-CoV-2 genome. 
It required bwa, freebayes, samtools, snpEff , SnpSift , bcftools

## Installation 

1. ``git clone https://github.com/dridk/Sars-CoV-2-NGS-pipeline.git``
2. ``conda env create -f environment.yml``
3. ``conda activate covid``
4: ``snakemake final.ann.vcf``

