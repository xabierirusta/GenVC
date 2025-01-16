# **GenVC Introduction**

**GenVC** is a Nextflow 24.10.2 pipeline that processes raw FASTQ files and calls SV, SNV and CNV variants using [Delly](https://github.com/dellytools/delly), GATK [HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)/[Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851-Mutect2) and [CNVKit](https://github.com/etal/cnvkit).

## Requirements

- Nextflow 24.10.2 version.
- All the processes are containerized in Docker images, so there is no need to install all the required tools.
  - However, in order to annotate the VCF files, **Ensembl Variant Effect Predictor (VEP)** should be installed in order to get the **homo_sapiens_merged_113_GRCh38 cache** that will be manually inputed using --vep_cache.
- Valid FASTQ files that can be checked using [fastq_check.py](fastq_check.py), paired-ended files have to be differentiated as ...{1,2}.fastq.gz
- targets.bed file adecuated to the analysis and its respective .gz and .gz.tbi files that have to follow the next structure:
    - 3 columns (chr_number target_start target_end) separated by tabs and the chromosome numbers have to be chr1,chr2,chX...
    - If all the exons are to be analyzed, [format_bedtools.sh](bed/format_bedtools.sh) script can be useful to format the [all_exons.txt](bed/all_exons.txt) file.
    - The new targets file HAS to be named **targets.bed**.

## Installation

Follow these steps to set up and run the pipeline locally:

### 1. Clone the repository

First, clone the repository to your local machine:

```bash
git clone https://github.com/xabierirusta/GenVC.git
cd GenVC
```

### 2. Preparation of targets.bed file

There is already a targets.bed file that includes all human exons ready to use in the repository but if other regions of interest were needed, a custom targets.bed and all the files needed could be created using [format_bedtools.sh](bed/format_bedtools.sh) script taking into account that the .txt used has to have the same structure as [all_exons.txt](bed/all_exons.txt) that can be obtained in [USCS Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables).

```bash
cd GenCV/bed
format_bedtools.sh <exons_file.txt>
```

### 3. Download required VCF files

Download the following files and save them inside a directory called **vcf** in the GenVC directory.

Use the [GATK](https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38;tab=objects?prefix=&forceOnObjectsSortingFiltering=false&inv=1&invt=Abm_VA) best practises google cloud bucket to download the following files.
- 1000g_pon.hg38.vcf.gz
- 1000g_pon.hg38.vcf.gz.tbi
- af-only-gnomad.hg38.vcf.gz
- af-only-gnomad.hg38.vcf.gz.tbi

Also use [ClinVar](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/) to get the ClinVar VCF files.
- clinvar.vcf.gz   
- clinvar.vcf.gz.tbi

Use the next [link](https://mega.nz/folder/MGFklDrA#CNUzj2MqgWE03UGFYA3YZA) to download the dbSNP database VCF file.

### 4. Download required Reference files

Download the following files and save them inside a directory called **ref** in the GenVC directory

Use [MEGA] to access the reference genome files and download the next directory.
- ref

But in order to obtain the primary assembly file go to the [Ensembl](https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/) page and download and decompress it:
- Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

Index the reference primary assembly FASTA file using:
```bash
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa

```

Use the next [link](https://mega.nz/folder/pX1HEIrK#AVwUR8ozUKaTqjwNngdh7Q) to download the rest of necessary reference files.
## Input

This pipeline only allows inputs in FASTQ format that have been obtained by Illumina sequencing.

## Output

This pipeline generates vcf files and their indexes for SV, SNV and CNV variants using Delly, GATK HaplotypeCaller/Mutect2 and CNVKit. Also some plots from the CNVKit tool that would require installing [CNVkit](https://github.com/etal/cnvkit) for further personalization.

## Pipeline Usage

```bash
cd GenVC
nextflow run GenVC --input_type <WGS/WES> --input path/to/fastq_files_dir --seq <s/p> --v <germline/somatic> --vep_cache HOME/.vep --outDir path/to/output_dir
```

**Optional arguments:**
- --help -->          Shows a help message and exits the pipeline.
- --ref -->           Directory where reference files are stored, default "./ref/".
- --vcf -->           Directory where vcf files are stored, default "./vcf/".
- --tools -->         Directory where files by tools are stored, default "./tools/".
- --t -->             Number of threads used in the pipeline.
- --normal -->        Path to the location directory of the fastq file that is going to be used as normal control (Compulsory when using --v somatic or --input_type WGS).
- --bed -->           Path to the directory where the BED files are located (For WES, target Bed file HAS to be named targets.bed).

## Contact
If you have any questions, feel free to reach out to me at:

- **Email:** xabi.irusta@gmail.com
- **GitHub:** [xabierirusta](https://github.com/xabierirusta)
