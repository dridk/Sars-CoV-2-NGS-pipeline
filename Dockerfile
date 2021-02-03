
FROM continuumio/miniconda3

WORKDIR /usr/app 

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge


ADD Snakefile /usr/app/Sars-CoV-2-NGS-pipeline/
ADD config.yml /usr/app/Sars-CoV-2-NGS-pipeline/
ADD report/ /usr/app/Sars-CoV-2-NGS-pipeline/report/

ADD data/fastq /fastq/
ADD data/genome /genome/


# Install env 
#RUN git clone https://github.com/dridk/Sars-CoV-2-NGS-pipeline


WORKDIR /usr/app/Sars-CoV-2-NGS-pipeline

RUN conda install bwa freebayes samtools snpEff

RUN python -m pip install snakemake


# # Install pangolin
# RUN git clone https://github.com/cov-lineages/pangolin.git
# WORKDIR /usr/app/pangolin
# RUN conda env create -f environment.yml
# RUN conda activate pangolin
# RUN python setup.py install

RUN ls /fastq
RUN ls /genome
ENTRYPOINT ["snakemake","-d","/fastq/", "-pFj4","--configfile","/usr/app/Sars-CoV-2-NGS-pipeline/config.yml","-C","FASTQ_DIR=/fastq", "GENOME=/genome/whuan.fasta"] 

#ENTRYPOINT ["snakemake", "A.bam", "-j4"]