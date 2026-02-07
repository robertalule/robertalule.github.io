---
layout: page
title: "RNA-Seq Pipeline: From Raw Reads to Gene Counts"
permalink: /rnaseq-pipeline/
order: 3
---

# RNA-Seq Pipeline Walkthrough: From Raw Reads to Gene Counts

This document provides a detailed, step-by-step walkthrough of a standard RNA-Seq analysis pipeline, from raw sequencing data (FASTQ files) to a gene-level count matrix. The pipeline is designed for a single-end, unstranded experiment using a yeast (*Saccharomyces cerevisiae*) reference. This workflow is intended for educational use and serves as a technical reference for reproducible analysis.

**Learning Objectives**
*   Understand the purpose and interpretation of quality control (QC) for RNA-Seq data.
*   Execute and justify read trimming and filtering.
*   Perform splice-aware alignment of RNA-Seq reads to a reference genome.
*   Generate a gene-level count matrix suitable for downstream differential expression analysis.

**Prerequisites**
*   A Linux-based command-line environment.
*   Installed software: FastQC, Trim Galore (which uses Cutadapt), HISAT2, samtools, featureCounts.
*   A reference genome (FASTA) and its gene annotation (GTF file).
*   Raw RNA-Seq reads in FASTQ format (`.fastq.gz`).

---

## 1. Data Organization and Initial Inspection

**Concept:** Before any analysis, establish a clean, documented directory structure. This prevents file mix-ups and ensures reproducibility.

**Biological Rationale:** RNA-Seq data is the primary measurement. Its initial quality dictates the validity of all downstream biological conclusions.

**Computational Logic:** A structured workspace allows for systematic processing and clear tracking of input, intermediate, and output files.

**Action:** Create a project directory and navigate to it.
```bash
mkdir -p rnaseq_project/{00_raw_data,01_fastqc,02_trimmed,03_aligned,04_counts}
cd rnaseq_project
```
Place your raw FASTQ files (e.g., `sample1.fastq.gz`) in `00_raw_data/`. For this guide, we assume a single sample file.

---

## 2. Pre-trimming Quality Control with FastQC

**Concept:** FastQC performs a series of diagnostic tests on raw sequence data to identify potential technical artifacts.

**Biological Rationale:** We must assess if the sequencing data is trustworthy. Problems like adapter contamination, severe sequence bias, or pervasive low-quality bases can invalidate biological interpretation.

**Computational Logic:** FastQC runs multiple independent modules, each producing a metric or plot. The goal is not a single "pass/fail" score but a holistic interpretation of all modules in the context of RNA-Seq.

**Action:** Run FastQC on the raw FASTQ file.
```bash
fastqc 00_raw_data/sample1.fastq.gz -o 01_fastqc/
```

**Interpretation & Decision:**
Open the generated HTML report (`sample1_fastqc.html`). Key modules and their interpretation for RNA-Seq are:

*   **Per base sequence quality:** A gradual decline in quality towards the 3' end of reads is normal for Illumina sequencing. A sharp drop early in the read indicates a problem. **Decision:** Proceed if quality is generally high (Q > 28).
*   **Per base sequence content:** The proportion of A/T/G/C at each position. RNA-seq libraries prepared with random hexamer priming **typically show a biased sequence composition at the first ~10-12 bases**. This is a technical artifact of the protocol, not a reason to discard data.
*   **Adapter content:** This detects the presence of sequencing adapter sequences within the reads. **Any non-zero adapter contamination is a problem that must be fixed by trimming**, as adapters can interfere with alignment.
*   **Sequence duplication levels:** In RNA-seq, a high level of duplication can reflect truly highly expressed genes (biology), not just PCR bias (artifact). This requires careful interpretation based on the organism's transcriptome size.

**Explicit Decision (based on attached log):**
The report for `SRR453566` showed:
*   **PASS** for base quality and GC content.
*   **FAIL/WARN** for per-base sequence content (expected RNA-seq bias).
*   **FAIL** for adapter content and overrepresented sequences.
**Conclusion:** The data is of good sequencing quality but contains adapter contamination. **Proceed to trimming.**

**Alternative:** While FastQC is the standard diagnostic tool, `multiQC` can aggregate reports from many samples into a single summary.

---

## 3. Read Trimming with Trim Galore

**Concept:** Trim Galore is a wrapper for `Cutadapt` that automates the removal of adapter sequences and low-quality bases from read ends.

**Biological Rationale:** Adapters are non-biological sequences that must be removed to allow the true biological sequence to align correctly to the genome. Trimming low-quality bases can improve alignment accuracy.

**Computational Logic:** The tool scans reads for matches to a known adapter sequence and removes it. It can also trim bases from the 3' end that fall below a specified quality threshold.

**Action:** Run Trim Galore. It will automatically detect common adapter sequences.
```bash
trim_galore \
  --output_dir 02_trimmed \
  --quality 20 \
  --length 25 \
  00_raw_data/sample1.fastq.gz
```

**Interpretation & Decision:**
Check the trimming report (`sample1.fastq.gz_trimming_report.txt`). It details the percentage of reads with adapter content and the length distribution after trimming.
Run FastQC again on the trimmed file (`02_trimmed/sample1_trimmed.fq.gz`).

**Explicit Decision (based on attached log):**
The post-trimming FastQC report showed:
*   **PASS** for **Adapter Content**. This is the critical check – the primary issue has been resolved.
*   **FAIL** for Per base sequence content and Sequence duplication levels.
*   **WARN** for Sequence length distribution.
**Interpretation:** The remaining "FAILs" are expected consequences of trimming and RNA-seq biology. The per-base content bias persists due to library prep. Duplication levels are expected in a compact yeast transcriptome. The length distribution warning confirms reads were variably trimmed, which is correct. **Conclusion: Trimming was successful. Proceed to alignment.**

**Alternative:** `fastp` is a modern, all-in-one alternative for QC, filtering, and trimming, offering superior speed.

---

## 4. Splice-Aware Alignment with HISAT2

**Concept:** HISAT2 aligns RNA-seq reads to a reference genome, accounting for reads that span exon-exon junctions (spliced alignments).

**Biological Rationale:** Eukaryotic genes contain introns that are removed during splicing. An RNA-seq read derived from a mature mRNA may align discontinuously across the junction where two exons join.

**Computational Logic:** HISAT2 uses a hierarchical graph FM-index to efficiently map reads across the genome, including to known and novel splice sites. It requires a pre-built genome index.

**Action:**
1.  **Build the genome index (once per reference).**
    ```bash
    hisat2-build reference_genome.fa hisat2_index_name
    ```
2.  **Align the trimmed reads.**
    ```bash
    hisat2 \
      -x path/to/hisat2_index_name \
      -U 02_trimmed/sample1_trimmed.fq.gz \
      --known-splicesite-infile known_splice_sites.txt \
      -p 2 \
      | samtools view -bS \
      | samtools sort -o 03_aligned/sample1.sorted.bam

    samtools index 03_aligned/sample1.sorted.bam
    ```
    *   `-x`: Path to the HISAT2 index.
    *   `-U`: Input file (single-end, unpaired reads).
    *   `--known-splicesite-infile`: File of known splice sites (optional, can improve alignment).
    *   The output is piped to `samtools` to convert SAM to BAM, sort by genomic coordinate, and create an index.

**Interpretation & Decision:** Check the alignment summary printed to the terminal by HISAT2. The key metric is the **overall alignment rate**.

**Explicit Decision (based on attached log):**
HISAT2 reported an overall alignment rate of **97.13%** for the yeast sample. This is an excellent mapping rate, indicating the reference and data are compatible. **Proceed to alignment QC.**

**Alternative:** **STAR** is a widely used, highly accurate aligner that is faster but requires more memory (RAM). For standard analyses, both HISAT2 and STAR are excellent choices.

---

## 5. Alignment Quality Control

**Concept:** Verify that the alignment process itself behaved as expected and that spliced alignments are present.

**Biological Rationale:** We need to confirm that reads are mapping to the genome and that the aligner correctly identified splicing events, confirming we have processed mRNA data.

**Computational Logic:** Use `samtools` to extract statistics and inspect individual alignments.

**Action:**
1.  **Get alignment statistics.**
    ```bash
    samtools flagstat 03_aligned/sample1.sorted.bam
    ```
2.  **Check for spliced alignments.**
    ```bash
    samtools view 03_aligned/sample1.sorted.bam | head -20
    ```
    Look for the `CIGAR` string containing an `N` (e.g., `100M` is unspliced, `50M1000N50M` is spliced).

**Interpretation & Decision:**
*   `flagstat`: Confirm the primary mapped percentage is high (e.g., >70-80% for yeast/mouse/human) and consistent with the HISAT2 summary.
*   Spliced reads: The presence of reads with `N` in the CIGAR string is a positive indicator of correct splice-aware alignment.

**Explicit Decision (based on attached log):**
`flagstat` confirmed **96.79% of primary reads mapped**. Inspection of alignments showed multiple reads with `N` in the CIGAR string (e.g., `12S60M`). **Conclusion: Alignment is successful and splice-aware. The BAM file is ready for quantification.**

---

## 6. Read Quantification with featureCounts

**Concept:** `featureCounts` counts the number of reads assigned to genomic features, such as genes, based on their alignment coordinates.

**Biological Rationale:** The fundamental output of RNA-Seq is a measure of gene expression abundance. This is achieved by counting how many reads originate from each gene.

**Computational Logic:** The tool overlaps the genomic coordinates of each aligned read with the coordinates of features (exons) defined in a GTF annotation file. Reads overlapping a unique gene are counted for that gene.

**Action:** Run featureCounts. The most critical parameter is strandedness (`-s`).
*   `-s 0`: unstranded (most common historic protocol)
*   `-s 1`: stranded, read same strand as gene
*   `-s 2`: reverse-stranded

For our unstranded yeast data:
```bash
featureCounts \
  -T 2 \
  -s 0 \
  -t exon \
  -g gene_id \
  -a reference_annotation.gtf \
  -o 04_counts/gene_counts.txt \
  03_aligned/sample1.sorted.bam
```

**Interpretation & Decision:**
1.  Examine the summary file: `cat 04_counts/gene_counts.txt.summary`.
    *   Focus on the "Assigned" vs. "Unassigned" counts. A high proportion of assigned reads (>70-80%) is good.
2.  The main output `gene_counts.txt` is a tab-delimited matrix where rows are genes and columns are samples. The key column is the raw read count.

**Explicit Decision (based on pipeline completion):**
The summary file will be checked for a high assignment rate. The final output is the file `04_counts/gene_counts.txt`. This **gene-level count matrix is the end goal of this pipeline stage** and serves as the direct input for downstream statistical analysis tools like DESeq2 or edgeR, which will handle normalization and differential expression testing.

**Alternative:** **HTSeq-count** (`htseq-count`) is a popular alternative with similar functionality. Subread's `featureCounts` is generally faster, especially for multi-threaded jobs.

---

## Pipeline Conclusion

This workflow has processed raw RNA-Seq reads through quality control, adapter trimming, splice-aware alignment, and gene-level quantification. The final product is a count matrix (`gene_counts.txt`) where each cell represents the number of reads mapping to a specific gene in a specific sample.

**Next Steps:** This count matrix is the required input for statistical analysis in R/Bioconductor packages (e.g., DESeq2, edgeR, limma-voom). These tools perform normalization to account for differences in library size and composition, and test for statistically significant differences in gene expression between experimental conditions.

**Reproducibility Note:** To reproduce this analysis exactly, record the versions of all software used (e.g., `fastqc --version`, `hisat2 --version`) and the key parameters chosen at each step, particularly the strandedness (`-s`) flag in `featureCounts`.
