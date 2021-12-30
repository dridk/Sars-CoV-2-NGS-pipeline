
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

RUN conda install mamba -c conda-forge

RUN mamba create -n pangolin -y -c bioconda -c conda-forge -c defaults "python==3.8" "pangolin>=3.1.17" bwa freebayes samtools snpEff snpsift bcftools tabix snakemake


#RUN conda install bwa freeibayes samtools snpEff bcftools tabix snakemake
#RUN conda install -c bioconda -c conda-forge -c defaults pangolin

#RUN python -m pip install snakemake

ENV PATH /opt/conda/envs/pangolin/bin:$PATH
ENV CONDA_DEFAULT_ENV pangolin


#RUN mamba activate pangolin

RUN snakemake --version
RUN python --version
RUN pangolin --version
# # Install pangolin
# RUN git clone https://github.com/cov-lineages/pangolin.git
# WORKDIR /usr/app/pangolin
# RUN conda env create -f environment.yml
# RUN conda activate pangolin
# RUN python setup.py install

RUN ls /fastq
RUN ls /genome
ENTRYPOINT ["snakemake","-d","/fastq/", "-pj4","--configfile","/usr/app/Sars-CoV-2-NGS-pipeline/config.yml","-C","FASTQ_DIR=/fastq", "GENOME=/genome/whuan.fasta"] 

#ENTRYPOINT ["snakemake", "A.bam", "-j4"]
