# VAT-Variant-Annotation-Tool
<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]



<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/github_username/repo_name">
    <img src="images/logo.png" alt="Logo" width="80" height="80">
  </a>

  <h3 align="center">VAT-Variant-Annotation-Tool</h3>

  <p align="center">
    VAT is a variant annotation tool that parses VCF files and fetches variant information from Ensembl Variant Effect Predictor (VEP)
    <br />
    <a href="https://github.com/github_username/repo_name/issues">Report Bug</a>
    ·
    <a href="https://github.com/github_username/repo_name/issues">Request Feature</a>
  </p>
</p>

<!-- TABLE OF CONTENTS -->
## Table of Contents

* [Introduction](#about-the-project)
* [Getting Started](#getting-started)
  * [Prerequisites](#prerequisites)
  * [Installation](#installation)
* [Usage](#usage)
	* [Example](#example)
	* [Output format] (#output)
* [Contact](#contact)



<!-- Introduction-->
## Introduction

VAT is a variant annotation tool that parses VCF files, annotates the following information and output annotation in a tsv file. 
1. Depth of sequence coverage at the site of variation.
2. Number of reads supporting the variant.
3. Percentage of reads supporting the variant versus those supporting reference reads.
4. Gene of the variant, type of variation (substitution, insertion, deleteion, CNV, etc.) and their effect (missense, silent, intergenic, etc.). The API
documentation is available here: https://rest.ensembl.org/#VEP
5. The minor allele frequency of the variant if available.
6. Additional annotations: rsid, minor_allele_freq, minor_allele, clin_sig, transcript_id, gene_id, impact, gene_symbol, biotype, polyphen_prediction, and sift_prediction. NOTE: there might be multiple transcript_id, gene_id, gene_symbol, polyphen_prediction, sift_prediction for each variant. 


<!-- GETTING STARTED -->
## Getting Started

To install VAT, python (>= 3.10), vcfpy, requests, json5 and jsonschema are required. 
### Prerequisites

```sh
conda create --name vat python=3.10
conda activate vat
pip install -r requirements.txt
```

### Installation

1. Clone the repo
```sh
git clone https://github.com/Jingwen7/VAT-Variant-Annotation-Tool.git
```


<!-- USAGE EXAMPLES -->
## Usage
```sh
usage: vat.py [-h] -i [INPUT_FILE_PATH] -o [OUTPUT_FILE_PATH]
```
### Example
```sh
python3 vat.py -i /project/mchaisso_100/cmb-16/jingwenr/tempus/test_vcf_data.txt -o /project/mchaisso_100/cmb-16/jingwenr/tempus/out.tsv
```
### Output format
VAT stores the annotation of each variant in a tsv file.


The following shows an example of the tsv file
```
CHROM	POS	REF	ALT	TYPE	total_read_depth	ref_read_depth	alt_read_depth	Ratio_supporting_reads_alt_vs_ref	rsid	most_severe_consequence	minor_allele_freq	minor_allele	clin_sig	transcript_id	gene_id	gene_symbol	impact	biotype	polyphen_prediction	sift_prediction
1	11090916	C	A	123	1	122	122.0	COSV68896102	missense_variant	NA	NA	NA	ENST00000400897,ENST00000607145	ENSG00000009724,ENSG00000271895	MASP2,RP4-635E18.8	MODERATE,MODIFIER	protein_coding,antisense	possibly_damaging,NA	deleterious,NA
1	1650832	A	G	314	163	151	0.9263803680981595	rs72909030	missense_variant	NA	NA	NA	ENST00000356200,ENST00000356937,ENST00000357760,ENST00000358779,ENST00000378633,ENST00000378635,ENST00000378638,ENST00000401096,ENST00000404249,ENST00000460465,ENST00000479362,ENST00000487462,ENST00000498810,ENST00000509982,ENST00000598846	ENSG00000008128,ENSG00000008128,ENSG00000008128,ENSG00000008128,ENSG00000008128,ENSG00000008128,ENSG00000008128,ENSG00000008128,ENSG00000008128,ENSG00000008128,ENSG00000008128,ENSG00000008128,ENSG00000008128,ENSG00000008128,ENSG00000268575	CDK11A,CDK11A,CDK11A,CDK11A,CDK11A,CDK11A,CDK11A,CDK11A,CDK11A,CDK11A,CDK11A,CDK11A,CDK11A,CDK11A,RP1-283E3.8	MODERATE,MODIFIER,MODERATE,MODERATE,MODERATE,MODERATE,MODERATE,MODIFIER,MODERATE,MODERATE,MODERATE,MODIFIER,MODIFIER,MODERATE,MODIFIER	protein_coding,retained_intron,protein_coding,protein_coding,protein_coding,protein_coding,protein_coding,protein_coding,protein_coding,nonsense_mediated_decay,protein_coding,retained_intron,retained_intron,nonsense_mediated_decay,processed_transcript	benign,NA,benign,benign,benign,benign,benign,NA,benign,benign,benign,NA,NA,benign,NA	tolerated_low_confidence,NA,tolerated_low_confidence,tolerated_low_confidence,tolerated_low_confidence,tolerated_low_confidence,tolerated_low_confidence,NA,tolerated_low_confidence,tolerated,tolerated_low_confidence,NA,NA,tolerated,NA
```
