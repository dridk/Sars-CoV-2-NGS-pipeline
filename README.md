# Sars-CoV-2-NGS-pipeline
A simple snakemake pipeline to call variant from NGS data of Sars-CoV-2 genome. 
It depends on [bwa](http://bio-bwa.sourceforge.net/), [freebayes](https://github.com/freebayes/freebayes), [SnpEff/SnpSift](https://pcingola.github.io/SnpEff/) and [samtools](http://www.htslib.org/)



## Installation 
You can use conda to install dependencies. 
1. ``git clone https://github.com/dridk/Sars-CoV-2-NGS-pipeline.git``
2. ``conda env create -f environment.yml``
3. ``conda activate covid``       

You can test the pipeline with our toys dataset : 

``snakemake -p A.results.csv B.results.csv  -j4``

## Configuration and execution 

From ```config.yml``` set the FASTQ_DIR to the folder containing fastq files. 
These files must follow the follwing pattern : 

- SAMPLENAME_1.fastq.gz
- SAMPLENAME_2.fastq.gz

To get result of a specific SAMPLENAME:

    snakemake -p SAMPLENAME.results.csv
    
To get fasta genom  of a specific SAMPLENAME:

    snakemake -p SAMPLENAME.fa
    
You can pass this consensus sequence  to [Pangolin](https://github.com/cov-lineages/pangolin) to get the lineage. 


## Results 

Each sample comes with a csv file with the following columns : 

1. Gene Name 
2. Feature ID 
3. Variant position
4. Reference bases 
5. Alternative bases
6. HGVS coding name 
7. HGVS protein name 
8. Impact 
9. effect 

```
ANN[*].GENE     ANN[*].FEATUREID        POS     REF     ALT     ANN[*].HGVS_C   ANN[*].HGVS_P   ANN[*].IMPACT   ANN[*].EFFECT
ORF1ab  GU280_gp01      490     T       A       c.225T>A        p.Asp75Glu      MODERATE        missense_variant
ORF1ab  YP_009725297.1  490     T       A       c.225T>A        p.Asp75Glu      MODERATE        missense_variant
ORF1ab  YP_009742608.1  490     T       A       c.225T>A        p.Asp75Glu      MODERATE        missense_variant
ORF1ab  GU280_gp01.2    490     T       A       c.225T>A        p.Asp75Glu      MODERATE        missense_variant
ORF1ab  YP_009725298.1  490     T       A       c.-316T>A               MODIFIER        upstream_gene_variant
ORF1ab  YP_009742609.1  490     T       A       c.-316T>A               MODIFIER        upstream_gene_variant
ORF1ab  YP_009725299.1  490     T       A       c.-2230T>A              MODIFIER        upstream_gene_variant
ORF1ab  YP_009742610.1  490     T       A       c.-2230T>A              MODIFIER        upstream_gene_variant
```
