# Dockerfiles

This repository contains Dockerfiles used to build Docker images for bioinformatics tools, with a focus on evolutionary biology and phylogenetics. All images are published on [DockerHub](https://hub.docker.com/r/evolbioinfo/) under the `evolbioinfo` organization.

A tentative of matching between the containerized tools and bio.tools is described in [BIOTOOLS.md](BIOTOOLS.md).

## Table of Contents

- [Usage](#usage)
  - [Pulling an image](#pulling-an-image)
  - [Running a container](#running-a-container)
  - [Mounting local data](#mounting-local-data)
- [Available Images](#available-images)
  - [Phylogenetics](#phylogenetics)
  - [Sequence Alignment](#sequence-alignment)
  - [Sequence Mapping](#sequence-mapping)
  - [Genomic Tools](#genomic-tools)
  - [Variant Calling](#variant-calling)
  - [Quality Control](#quality-control)
  - [Read Trimming](#read-trimming)
  - [Assembly](#assembly)
  - [Molecular Evolution](#molecular-evolution)
  - [Population Genetics](#population-genetics)
  - [Gene Expression](#gene-expression)
  - [Taxonomy & Metagenomics](#taxonomy--metagenomics)
  - [Viral & Pandemic Analysis](#viral--pandemic-analysis)
  - [Ancient DNA](#ancient-dna)
  - [Sequence Analysis](#sequence-analysis)
  - [Haplotyping](#haplotyping)
  - [Simulation](#simulation)
  - [Read Correction](#read-correction)
  - [Visualization](#visualization)
  - [Workflow Management](#workflow-management)
  - [Base Images](#base-images)
  - [Utilities](#utilities)
- [Building Images Locally](#building-images-locally)
- [Contributing](#contributing)

## Usage

### Pulling an image

Each image is hosted on DockerHub under `evolbioinfo/<tool-name>:<version>`. To pull an image:

```bash
docker pull evolbioinfo/<tool-name>:<version>
```

For example, to pull the RAxML-NG image:

```bash
docker pull evolbioinfo/raxml-ng:v1.2.2
```

### Running a container

Most images define an `ENTRYPOINT` pointing directly to the tool executable. You can run a tool as follows:

```bash
docker run evolbioinfo/<tool-name>:<version> [tool options]
```

For example, to display the RAxML-NG help:

```bash
docker run evolbioinfo/raxml-ng:v1.2.2 --help
```

Or to run FastTree:

```bash
docker run evolbioinfo/fasttree:v2.1.9 -help
```

### Mounting local data

To use your local files inside the container, mount your working directory using the `-v` flag:

```bash
docker run -v /path/to/your/data:/data evolbioinfo/<tool-name>:<version> [tool options]
```

For example, to run RAxML-NG on a local alignment file:

```bash
docker run -v $(pwd):/data evolbioinfo/raxml-ng:v1.2.2 --msa /data/alignment.fasta --model GTR+G --prefix /data/output
```

## Available Images

The versions listed represent the latest containerized builds in this repository, not necessarily the most recent versions of the underlying tools.

### Phylogenetics

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [Bio++](biopp/) | v3.0.0 |  The Bio++ Libraries for phylogenetic and sequence analysis |
| [ebg](ebg/) | v0.13.3 | [Educated Bootstrap Guesser](https://github.com/wiegertj/EBG/wiki) |
| [epa-ng](epa-ng/) | v0.3.8 | Evolutionary placement algorithm for short reads into reference trees |
| [fastme](fastme/) | v2.1.6.4 | Fast and accurate distance-based phylogenetic tree construction |
| [fasttree](fasttree/) | v2.1.9 | Approximate maximum-likelihood phylogenetic trees for large alignments |
| [goalign](goalign/) | v0.4.0 | Multiple sequence alignment analysis toolkit |
| [gotree](gotree/) | v0.5.1 | Phylogenetic tree manipulation toolkit |
| [guppy](guppy/) | v3.1.5 | Tools for working with pplacer phylogenetic placement files |
| [iqtree](iqtree/) | v3.1.0 | Fast and accurate maximum-likelihood phylogenetic tree inference |
| [lsd](lsd/) | v0.3.3 | Least-squares dating of phylogenetic trees |
| [lsd2](lsd2/) | v2.4.1 | Least-squares dating of phylogenetic trees (version 2) |
| [maple](maple/) | v0.7.5 | Maximum likelihood phylogenetic estimation with reduced memory |
| [ml_bootstrap](ml_bootstrap/) | fec985c | Machine learning based support (see [this article](https://academic.oup.com/bioinformatics/article/40/Supplement_1/i208/7700891)) |
| [mrbayes](mrbayes/) | v3.2.7 | Bayesian inference of phylogenetic trees |
| [newick_utilities](newick_utilities/) | v1.6 | Utilities for manipulating Newick format trees |
| [ngphylogeny_multitools](ngphylogeny_multitools/) | seqtype_detect | Multi-tool image for phylogenetic workflows |
| [phylodeep](phylodeep/) | v0.9 | Deep learning for phylodynamic parameter estimation |
| [phyml](phyml/) | v3.3.20250515 | Maximum-likelihood phylogenetic tree estimation |
| [phyml-sms](phyml-sms/) | v1.8.1.1 | PhyML with Smart Model Selection |
| [ptp](ptp/) | v4bb2daf | Species delimitation from phylogenetic trees |
| [rappas](rappas/) | v1.21 | Rapid alignment-free phylogenetic identification via statistical hypothesis testing |
| [raxml](raxml/) | v8.2.13 | Randomized accelerated maximum likelihood phylogenetic inference |
| [raxml-ng](raxml-ng/) | v2.0.0 | RAxML next-generation |
| [table2itol](table2itol/) | latest | Converts annotation tables to iTOL dataset files |
| [tqdist](tqdist/) | v1.0.2 | Computing quartet and triplet distances between trees |
| [treedater](treedater/) | 89a0df0 | Scalable relaxed clock phylogenetic dating |
| [treemmer](treemmer/) | v0.3 | Reduce tree size while preserving phylogenetic diversity |
| [treesimulator](treesimulator/) | v0.2.27 | Simulating rooted phylogenetic trees under various models |
| [treestructure](treestructure/) | a831a66 | [Identification of hidden population structure in time-scaled phylogenies](https://academic.oup.com/sysbio/article/69/5/884/5734655) |
| [treetime](treetime/) | v0.11.4 | Maximum-likelihood phylogenetic time-trees inference |
| [treewas](treewas/) | v1.0 | Genome-wide association studies on phylogenetic trees |

### Sequence Alignment

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [bmge](bmge/) | v1.12 | Block Mapping and Gathering with Entropy for alignment trimming |
| [clustal_omega](clustal_omega/) | v1.2.4 | Fast and scalable multiple sequence alignment |
| [gblocks](gblocks/) | v0.91b | Alignment trimming by selecting conserved blocks |
| [lastal](lastal/) | v980 | Local alignment of biological sequences |
| [mafft](mafft/) | v7.526 | Multiple sequence alignment using fast Fourier transform |
| [muscle](muscle/) | v3.8.31 | Multiple sequence alignment |
| [noisy](noisy/) | v1.5.12 | Identify homoplastic characters in multiple sequence alignments |
| [tcoffee](tcoffee/) | Version_11.00.18778a8 | Multiple sequence alignment using T-Coffee |
| [trimal](trimal/) | v1.5.1 | Automated removal of spurious sequences or poorly aligned regions |

### Sequence Mapping

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [bbmap](bbmap/) | v39.01 | Short read aligner for DNA and RNA-seq data |
| [bowtie](bowtie/) | v1.3.1 | Aligning sequencing reads to references |
| [bowtie2](bowtie2/) | v2.5.5 | Mapping DNA sequences against a large reference genome |
| [bwa](bwa/) | v0.7.19 | Burrows-Wheeler aligner for short DNA sequences |
| [hisat2](hisat2/) | v2.2.2 | Alignment program for mapping next-generation sequencing reads to a population of human genomes |
| [mash](mash/) | v2.3 | Fast genome and metagenome distance estimation using MinHash |
| [minimap2](minimap2/) | v2.30 | Versatile pairwise aligner for genomic and spliced nucleotide sequences |
| [papara](papara/) | v2.5 | Phylogeny-aware short-read alignment |
| [star](star/) | v2.7.11b | Spliced Transcripts Alignment to a Reference (RNA-seq) |

### Genomic Tools

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [alfred](alfred/) | v0.5.3 | BAM alignment statistics, feature counting and feature annotation |
| [bam-readcount](bam-readcount/) | v1.0.1 | Per-position read counts from BAM files |
| [bamUtil](bamUtil/) | v1.0.15 | Programs for working on SAM/BAM files |
| [bedtools](bedtools/) | v2.31.1 | Genome arithmetic and interval manipulation |
| [dsrc](dsrc/) | v2.0.2 | DNA sequence compression tool |
| [fastqutils](fastqutils/) | v0.1.7 | Utilities for manipulating FASTQ files |
| [fastxtoolkit](fastxtoolkit/) | v0.0.14 | FASTX toolkit for preprocessing FASTQ/FASTA files |
| [gofasta](gofasta/) | v1.2.3 | Command-line utilities for working with genomic alignments |
| [picard](picard/) | v3.4.0 | Command-line tools for manipulating high-throughput sequencing data |
| [samtools](samtools/) | v1.23 | Reading, writing, and manipulating SAM/BAM/CRAM files |
| [seqkit](seqkit/) | v2.4.0 | Ultrafast toolkit for FASTA/Q file manipulation |
| [seqtk](seqtk/) | v1.5 | Toolkit for processing sequences in FASTA/Q formats |
| [sra-tools](sra-tools/) | v3.0.1 | NCBI SRA toolkit for downloading and processing sequencing data |
| [sratoolkit](sratoolkit/) | v3.0.1 | Alternate NCBI SRA toolkit image |
| [vcftools](vcftools/) | v0.1.16 | Tools for working with VCF files |

### Variant Calling

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [bcftools](bcftools/) | v1.23 | VCF/BCF variant manipulation and calling |
| [freebayes](freebayes/) | v1.3.10 | Bayesian genetic variant detector |
| [ivar](ivar/) | v1.4.4 | Tools for viral amplicon-based sequencing |

### Quality Control

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [catch](catch/) | v1.5.2 | Compact Aggregation of Targets for Comprehensive Hybridization |
| [fastqc](fastqc/) | v0.12.1 | Quality control analysis of high-throughput sequencing data |
| [minionqc](minionqc/) | v1.4.2 | Quality control for Oxford Nanopore sequencing data |
| [multiqc](multiqc/) | v1.9 | Aggregate bioinformatics results across samples into a report |
| [nanoplot](nanoplot/) | v1.46.2 | Plotting tools for long-read sequencing data |
| [rna-seqc](rna-seqc/) | v1.1.9 | Quality control metrics for RNA-seq data |

### Read Trimming

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [adapterremoval](adapterremoval/) | v2.3.3 | Trimming of adapters and low-quality bases from NGS reads |
| [alien_trimmer](alien_trimmer/) | v2.1 | Adapter trimming for sequencing reads |
| [trimgalore](trimgalore/) | v0.6.11 | Wrapper for Cutadapt and FastQC for adapter trimming |

### Assembly

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [canu](canu/) | v2.3 | Long-read assembler |
| [savage](savage/) | v0.4.1 | Sequence assembly for viral genomes |
| [spades](spades/) | v3.15.4 | Assembly and analysis of sequencing data |
| [velvet](velvet/) | v1.2.10 | De novo genomic assembler |

### Molecular Evolution

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [bayestraits](bayestraits/) | v5.0.0 | Bayesian analysis of trait evolution on phylogenies |
| [fastcodeml](fastcodeml/) | v1.1.0 | Accelerated codeml for detecting positive selection |
| [hyphy](hyphy/) | v2.5.96 | Hypothesis testing using phylogenies |
| [paml](paml/) | v4.8a | Phylogenetic analysis by maximum likelihood |
| [pcoc](pcoc/) | v898c138 | Detection of Convergent Amino-Acid Evolution |
| [pastml](pastml/) | v1.9.51 | Ancestral state reconstruction and phylogeographic inference |


### Population Genetics

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [admixture](admixture/) | v1.3.1 | Maximum-likelihood estimation of individual ancestries |
| [finestructure](finestructure/) | v4.1.1 | Population structure inference using haplotypes |

### Gene Expression

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [deseq](deseq/) | v1.39.0 | Differential expression analysis from RNA-seq count data |
| [stringtie](stringtie/) | v2.2.1 | Transcript assembly and quantification for RNA-seq |
| [subread](subread/) | v2.0.3 | Subread/featureCounts read summarization for RNA-seq |

### Taxonomy & Metagenomics

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [checkm](checkm/) | v1.2.5 | Quality assessment of genome bins from metagenomes |
| [fastani](fastani/) | v1.34 | Fast and accurate whole-genome ANI estimation |
| [khmer](khmer/) | v2.1.2 | Probabilistic k-mer counting data structure |
| [kraken](kraken/) | v2.17.1 | Taxonomic classification of metagenomic sequences |
| [krakenuniq](krakenuniq/) | v1.0.4 | Metagenomics classification using unique k-mer counts |
| [vamb](vamb/) | v3.0.2 | Variational autoencoders for metagenomic binning |

### Viral & Pandemic Analysis

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [artic-ncov2019](artic-ncov2019/) | e814ed4 | ARTIC network bioinformatics tools for SARS-CoV-2 |
| [civet](civet/) | v2.1.2 | Cluster investigation and virus epidemiology tool |
| [irma](irma/) | v1.0.3 | Iterative refinement meta-assembler for viral genomics |
| [label](label/) | v0.6.4 | Sequence labeling and annotation tool |
| [nextstrain-base](nextstrain-base/) | build-20251119T000157Z | Nextstrain base environment for viral phylodynamics |
| [pangolin](pangolin/) | v4.3.1 | Phylogenetic assignment of named global outbreak LINeages |
| [polecat](polecat/) | b4a36f3 | Phylogenetic Overview & Local Epidemiological Cluster Analysis Tool |
| [vivan](vivan/) | v0.43 | Virus variation analyzer |

### Ancient DNA

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [damageprofiler](damageprofiler/) | v1.1 | Profiling damage patterns in ancient DNA reads |
| [mapdamage](mapdamage/) | v2.2.3 | Identifying and quantifying DNA damage in ancient DNA |
| [pathphynder](pathphynder/) | v1.2.4 | Ancient DNA placement into reference phylogenies |
| [schmutzi](schmutzi/) | v1.5.6 | Estimation of ancient DNA contamination |

### Sequence Analysis

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [cd-hit](cd-hit/) | v4.8.1 | Sequence clustering|
| [gubbins](gubbins/) | v3.4 | Rapid detection of recombination in bacterial genomes |
| [hmmer](hmmer/) | v3.3 | Biosequence analysis using profile hidden Markov models |
| [jphmm](jphmm/) | v03.2015 | Jumping profile hidden Markov model for HIV subtyping |
| [jphmm_tools](jphmm_tools/) | v0.1.4 | Tools for working with jpHMM output |
| [sdrmhunter](sdrmhunter/) | v0.2.1.6 | HIV surveillance drug resistance mutation identification |

### Haplotyping

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [haploconduct](haploconduct/) | v0.2.1 | Haplotype-aware genome assembly toolkit |
| [haplogrep](haplogrep/) | v2.4.0 | Mitochondrial haplogroup classification |
| [predicthaplo](predicthaplo/) | v1.0 | Predicting HIV haplotypes from next-generation sequencing |
| [shorah](shorah/) | v1.99.2 | Short Reads Assembly into Haplotypes |
| [strainline](strainline/) | commit-8af032906e | Full-length de novo viral haplotype reconstruction |

### Simulation

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [indelible](indelible/) | v1.03 | Flexible evolutionary sequence simulator |
| [nanosim](nanosim/) | v3.2.3 | Nanopore sequence read simulator |
| [seq-gen](seq-gen/) | v1.3.4 | Simulation of molecular sequence data along phylogenetic trees |
| [snag](snag/) | master | Sequence simulation along a tree |
| [reseq](reseq/) | 053b8d1 | Realistic simulation of Illumina sequencing data |

### Read Correction

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [khmer](khmer/) | v2.1.2 | Probabilistic k-mer counting data structure |
| [musket](musket/) | v1.1 | k-spectrum based short read error correction |

### Visualization

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [igv](igv/) | v2.9.0 | Integrative Genomics Viewer for alignment and variant data |
| [inkscape](inkscape/) | latest | Vector graphics editor |

### Workflow Management

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [snakemake](snakemake/) | v5.6.0 | Workflow management system for reproducible bioinformatics |

### Base Images

These images serve as base environments for building other images or running custom analyses:

| Image | Latest Version | Description |
|-------|---------------|-------------|
| [perl](perl/) | v5.32.1 | Perl with BioPerl modules |
| [python](python/) | v3.8.2 | Python base image |
| [python-dl](python-dl/) | v3.13 | Python with deep learning packages (TensorFlow, PyTorch) |
| [python-evol](python-evol/) | v3.8.2 | Python with evolutionary biology packages |
| [python-ml](python-ml/) | v3.8.2 | Python with machine learning packages |
| [r-base](r-base/) | v4.0.2 | R statistical computing base image |
| [r-evol](r-evol/) | v4.2.2 | R with evolutionary biology packages |
| [r-extended](r-extended/) | v4.3.3 | R with extended bioinformatics packages |
| [r-gisaid](r-gisaid/) | v4.1.2 | R environment for GISAID data analysis |
| [r-sra](r-sra/) | v3.6.1 | R environment for SRA data analysis |
| [ubuntu](ubuntu/) | v24.04 | Ubuntu base image |

### Utilities

| Tool | Latest Version | Description |
|------|---------------|-------------|
| [jq](jq/) | v1.8.1 | Lightweight and flexible command-line JSON processor |
| [s3cmd](s3cmd/) | v2.4.0 | Command Line S3 Client and Backup for Linux and Mac |
| [s3utils](s3utils/) | v0.6.1 | Utilities for interacting with Amazon S3 |
| [sphinx](sphinx/) | v1.8.5 | Python documentation generator |
| [wget](wget/) | v1.17.1 | Network utility to retrieve files from the Web |

## Building Images Locally

To build a Docker image from this repository, navigate to the tool's version directory and run `docker build`:

```bash
cd <tool-name>/<version>/
docker build -t evolbioinfo/<tool-name>:<version> .
```

For example, to build the RAxML-NG image:

```bash
cd raxml-ng/v1.2.2/
docker build -t evolbioinfo/raxml-ng:v1.2.2 .
```

## Contributing

Contributions of new Dockerfiles or updates to existing ones are welcome. Please follow the conventions used in this repository:

- Place each tool's Dockerfile in a subdirectory named `<tool-name>/<version>/`.
- Use a `LABEL maintainer=` instruction to identify the image author.
- Use `ENV VERSION=<version>` to specify the tool version.
- Clean up build dependencies and package manager caches at the end of the `RUN` step.
- Define an `ENTRYPOINT` pointing to the tool's main executable where appropriate.
- Create a `/pasteur` directory at the end of the build (this is the shared mount point convention used in these images).

