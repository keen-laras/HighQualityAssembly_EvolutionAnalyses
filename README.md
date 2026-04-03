# High Quality Genome Assembly & Annotation + Evolution Analyses
A guide on genome assembly study. The whole pipeline contains assembly - polishing - scaffolding - repeat annotation - genome annotation and how to use the data for evolution analyses

## 1. Genome Assembly [🔗](https://github.com/keen-laras/HighQualityAssembly_EvolutionAnalyses/tree/main/1.Genome%20Assembly)
### High-quality genome assembly from long-read sequencing data requires efficient error correction and accurate graph-based assembly.

  This project uses **NextDenovo** to assemble a draft genome from raw **CycloneSEQ** reads by:

  a. Performing read correction to reduce sequencing errors

  b. Constructing a consensus assembly from corrected reads

  c. Running in parallel for improved performance and scalability

## 2. Genome Polish [🔗](https://github.com/keen-laras/HighQualityAssembly_EvolutionAnalyses/tree/main/2.Genome%20Polishing)
### Genome assemblies generated from long-read sequencing technologies (e.g., CycloneSEQ or PacBio) often contain base-level errors such as insertions, deletions, and mismatches.

  This pipeline applies **NextPolish** to iteratively refine a draft genome by:

  a. Correcting errors using short-read data (high accuracy)

  b. Refining structure using long-read data (long-range continuity)

  c. Running multiple polishing rounds for optimal results

## 3. Scaffolding [🔗](https://github.com/keen-laras/HighQualityAssembly_EvolutionAnalyses/tree/main/3.Scaffolding)
### Draft genome assemblies are often fragmented into contigs. To achieve chromosome-level assemblies, long-range interaction data such as Hi-C or Pore-C can be used to order and orient these contigs.

  This pipeline performs:

  a. Hi-C reads alignment **(Bowtie2)** for contact data by **HiC-Pro or Chromap**

  b. Scaffolding using **Hi-C/Pore-C** contact data

  c. Contig ordering and orientation with **YaHS**

  d. Contact map generation for visualization

  e. Manual curation support via **Juicebox** (.hic format)

## 4. Repeat Annotation [🔗](https://github.com/keen-laras/HighQualityAssembly_EvolutionAnalyses/tree/main/4.Repeat%20Annotation) 
### Repetitive elements are a major component of eukaryotic genomes and can interfere with downstream analyses such as gene prediction and alignment.

  This pipeline performs:

  a. Repeat annotation using **RepeatMasker** (multiple libraries)

  b. Soft-masked genome processing (lowercase regions)

  c. Extraction of repeat intervals (BED format)

  d. Merging repeat annotations from multiple sources

  e. Final hard-masked genome generation (N-masked FASTA)

## 5. Genome Annotation [🔗](https://github.com/keen-laras/HighQualityAssembly_EvolutionAnalyses/tree/main/5.Genome%20Annotation) 
### Accurate genome annotation requires combining multiple lines of evidence.

  This pipeline integrates:

  a. **BRAKER** → predicts genes directly from genomic sequence

  b. **EviAnn** → homology-based annotation using protein evidence from related species

## 6. Phylogenetic Tree [🔗](https://github.com/keen-laras/HighQualityAssembly_EvolutionAnalyses/tree/main/6.Phylogenetic%20Tree)
### This directory contains the pipeline used to construct gene-based phylogenetic trees for comparative genomic analyses.

  The workflow is designed for:

  a. Orthologous gene sets

  b. Coding sequence (CDS) alignment

## 7. QUIBL-ABBA Analyses [🔗](https://github.com/keen-laras/HighQualityAssembly_EvolutionAnalyses/tree/main/7.QUIBL-ABBA%20Analyses)
### This module focuses on detecting and characterizing introgression (gene flow between species) using two complementary approaches:

  a. **ABBA-BABA (D-statistics)** — tests for the presence of introgression

  b. **QuIBL (Quartet-based Inference of Branch Lengths)** — distinguishes introgression from incomplete lineage sorting (ILS)

  Together, these methods provide a robust framework for identifying discordant phylogenetic signals and inferring their evolutionary causes.

## 8. Synteny Analyses [🔗](https://github.com/keen-laras/HighQualityAssembly_EvolutionAnalyses/tree/main/8.Synteny%20Analyses)

  This module performs genome-wide synteny analysis and visualization to investigate chromosomal conservation and structural evolution  across species.

  It integrates:

  a. **GeneSpace** → detection of syntenic blocks and conserved genomic regions

  b. **GenomeSyn** → visualization of synteny relationships across genomes
  
