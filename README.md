# VAT-Variant-Annotation-Tool

VAT-Variant-Annotation-Tool
<p align="left">
VAT is a variant annotation tool that parses VCF files and fetches variant information from Ensembl Variant Effect Predictor (VEP).
<br />
<a href="https://github.com/github_username/repo_name/issues">Report Bug</a>
·
<a href="https://github.com/github_username/repo_name/issues">Request Feature</a>
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

VAT is a variant annotation tool that parses VCF files, annotates the following information and outputs annotation to a tsv file. 
1. Depth of sequence coverage at the site of variation.
2. Number of reads supporting the variant.
3. Percentage of reads supporting the variant versus those supporting reference reads.
4. Gene of the variant, type of variation (substitution, insertion, deleteion, CNV, etc.) and their effect (missense, silent, intergenic, etc.). The API
documentation is available here: https://rest.ensembl.org/#VEP
5. The minor allele frequency of the variant if available.
6. Additional annotations: genotype, rsid, minor_allele_freq, minor_allele, clin_sig, transcript_id, gene_id, impact, gene_symbol, biotype, polyphen_prediction, and sift_prediction. NOTE: there might be multiple transcript_id, gene_id, gene_symbol, polyphen_prediction, sift_prediction for each variant. 


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
cd VAT-Variant-Annotation-Tool/
```


<!-- USAGE EXAMPLES -->
## Usage
```sh
usage: vat.py [-h] -i [INPUT_FILE_PATH] -o [OUTPUT_FILE_PATH]
```
### Command line
```sh
python3 vat.py -i data/test_vcf_data.vcf -o data/out.tsv
```
### Input format
A typical input vcf file is provided: https://github.com/Jingwen7/VAT-Variant-Annotation-Tool/blob/main/data/test_vcf_data.vcf

### Output format
VAT outputs the annotation of variants to a tsv file. Each row stores the annotation of a variant. 


The following shows an example of the tsv file
```
CHROM	POS	REF	ALT	GT	TYPE	total_read_depth	ref_read_depth	alt_read_depth	Ratio_supporting_reads_alt_vs_ref	rsid	most_severe_consequence	minor_allele_freq	minor_allele	clin_sig	transcript_id	gene_id	gene_symbol	impact	biotype	polyphen_prediction	sift_prediction
1	11090916	C	A	1/1	123	1	122	122.0	COSV68896102	missense_variant	NA	NA	NA	ENST00000400897,ENST00000607145	ENSG00000009724,ENSG00000271895	MASP2,RP4-635E18.8	MODERATE,MODIFIER	protein_coding,antisense	possibly_damaging,NA	deleterious,NA
1	11087524	G	A	1/1	156	2	154	77.0	rs1782455	synonymous_variant	0.3125	G	['benign']	ENST00000240185,ENST00000315091,ENST00000400897,ENST00000439080,ENST00000473869,ENST00000477447,ENST00000480464,ENST00000496840,ENST00000607145	ENSG00000120948,ENSG00000120948,ENSG00000009724,ENSG00000120948,ENSG00000120948,ENSG00000120948,ENSG00000120948,ENSG00000120948,ENSG00000271895	TARDBP,TARDBP,MASP2,TARDBP,TARDBP,TARDBP,TARDBP,TARDBP,RP4-635E18.8	MODIFIER,MODIFIER,LOW,MODIFIER,MODIFIER,MODIFIER,MODIFIER,MODIFIER,MODIFIER	protein_coding,protein_coding,protein_coding,protein_coding,nonsense_mediated_decay,nonsense_mediated_decay,processed_transcript,nonsense_mediated_decay,antisense	NA,NA,NA,NA,NA,NA,NA,NA,NA	NA,NA,NA,NA,NA,NA,NA,NA,NA
```

The following are the decription of each column in output tsv file. 
| Column Name                       | Description                                                                                  |
|-----------------------------------|----------------------------------------------------------------------------------------------|
| CHROM                             | -                                                                                            |
| POS                               | -                                                                                            |
| REF                               | -                                                                                            |
| ALT                               | -                                                                                            |
| TYPE                              | Type of the variant (ins, del, snp)                                                          |
| GT                                | Genotype of the variant                                                                      |
| total_read_depth                  | Total number of reads at this site                                                           |
| ref_read_depth                    | Total number of reads containing REF allele                                                  |
| alt_read_depth                    | Total number of reads contains ALT alleles                                                   |
| Ratio_supporting_reads_alt_vs_ref | Ratio of number of reads supporting ALT versus those supporting REF allele.                  |
| rsid                              | Unique label of snp                                                                          |
| most_severe_consequence           | -                                                                                            |
| minor_allele_freq                 | -                                                                                            |
| minor_allele                      | -                                                                                            |
| clin_sig                          | Interpretations of variant's significance to disease (https://www.ncbi.nlm.nih.gov/clinvar/) |
| transcript_id                     | -                                                                                            |
| gene_id                           | -                                                                                            |
| gene_symbol                       | -                                                                                            |
| impact                            | Simple assessment of the putative impact of the variant based on SnpEff                      |
| biotype                           | Gene biotype                                                                                 |
| polyphen_prediction               | Potential pathogenicity of a variant                                                         |
| sift_prediction                   | Prediction of whether an amino acid substitution is likely to affect protein function        |
