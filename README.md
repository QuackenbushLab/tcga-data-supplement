# Supplememtary code for the paper "Reproducible processing of TCGA regulatory networks"


In this repo we host the code to generate the data and figures for the paper
["Reproducible processing of TCGA regulatory networks"](https://www.biorxiv.org/content/10.1101/2024.11.05.622163v1). 

All the data is generated with the [tcga-data-nf](https://github.com/QuackenbushLab/tcga-data-nf) workflow. This folder holds sample files and analyses that can be
run thanks to the pipeline. 

## Content of the repo


```
.
├── LICENSE
├── README.md
├── config # sample configuration files
├── data
│   ├── conf
│   │   └── coad-subtype/ # configuration files for the COAD subtype application in the paper
│   └── external
│       ├── coad-subtype # subtype assignment for each TCGA-COAD sample
│       └── reactome_slim # reactome SLIM pathways used in the paper
├── envs # conda environments
├── notebooks
│   ├── colon_subtype_dragon.ipynb # DRAGON results in the paper
│   ├── colon_subtype_panda.ipynb # PANDA results in the paper
│   └── src # reusable functions
└── results # folder where all results are generated
```

## Full pipeline for TCGA-COAD subtype DRAGONs and PANDAs

First, we ran the full `tcga-data-nf` workflow with the configuration in 
`coad_subtype.config` and the metadata in  `full_coad_subtypes.json`.

```
$ nextflow run tcga-data-nf -profile conda --pipeline full -c coad_subtype.config
```

Results are stored into the `results/batch-coad-subtype-20240510/` folder which has the following structure:

```
├── tcga_coad_cms1
│   ├── analysis
│   │   ├── dragon
│   │   └── panda
│   ├── data_download
│   │   ├── clinical
│   │   ├── cnv
│   │   ├── methylation
│   │   ├── mutations
│   │   └── recount3
│   └── data_prepared
│       ├── methylation
│       └── recount3
├── tcga_coad_cms2
│   ...
├── tcga_coad_cms3
│   ...
└── tcga_coad_cms4
    ...
```

For each subtype, you'll find the downloaded data (`data_download`), the prepared data (`data_prepared`) and the
networks (`analysis`).

## Notebooks for TCGA-COAD subtype DRAGONs and PANDAs
aw
The notebooks reproduce the results in the paper. In order to run the code in them, you need to have the pre-processed
DRAGON and PANDA networks. 

You can either download the `batch-coad-subtype-20240510` folder, or run the workflow again to generate all the data. 

## Static version of the data

The data relative to this repo can be found on the 
Harvard Dataverse: Replication Data for: tcga-data-nf

```
@data{DVN/MCSSYJ_2024,
author = {Fanfani, Viola},
publisher = {Harvard Dataverse},
title = {{Replication Data for: tcga-data-nf}},
UNF = {UNF:6:TYixGNR1fJyPs/vReFVaPQ==},
year = {2024},
version = {V1},
doi = {10.7910/DVN/MCSSYJ},
url = {https://doi.org/10.7910/DVN/MCSSYJ}
}
```

## Examples: Data and configuration files for 12 common cancers

Data on AWS:
[tcga-data-nf-procumputed](https://us-east-2.console.aws.amazon.com/s3/buckets/tcga-data-nf-precomputed?region=us-east-2&bucketType=general&tab=objects).

In order to visualize and download this data, you need to have an active AWS account (a free tier one should suffice). 
For any additional help, please contact vfanfani@hsph.harvard.edu


We'll keep an updated list of exemplary configuration files inside the `config` folder.

**For the most updated structure of the configuration files always refer to the tests inside the tcga-data-nf repository**

### Full analysis

For examples of configuration files for a full analysis you can refer to those we used for the colon cancer application:  
1. Pipeline configurations: `data/conf/coad-subtype/coad_subtype.config`
2. Data configurations: `data/conf/coad-subtype/full_coad_subtypes.json`

### TCGA Downloads

We paste here the configuration files we used to download data from TCGA. These are also available alongside the data on
AWS.

**First round downloads**:

::warning:: These configuration files follow an older structure of the metadata, but they still include all relevant
information to understand what has been downloaded

1. Clinical data:  `config/download_clinical_tcgabiolinks_firstround.config`
2. Gene Expression: `config/download_expression_recount3_firstround.config`
3. Mutations: `config/download_mutation_tcgabiolinks_firstround.config`
4. Methylation: `config/download_methylation_firstround.config`

Files are at: 

**New Methylation**:

GDC data went through some ID changes/downgrading to legacy, so we re-downloaded and prepared all methylation data:
Configuration file: `conf/download_methylation.json`


### Prepare

#### Prepare Gene expression

We have pre-processed gene expression data for the following tumor types: BRCA, COAD, DLBC, KIRC, LAML, LIHC, PRAD, PAAD,
SKCM, STAD, LUAD, LUSC.

Configuration file (tcga-data-nf (0.0.10)): conf/expression_prepare.conf


Output files follow the naming:

> recount3_tcga_coad_purity06_normlogtpm_mintpm1_fracsamples000001_tissuetumor_batchtcgagdcplatform_adjtcgagdcplatform.txt
 where we write in the filename the parameters used to generate it. 

For instance, the file above is in logptm, has genes with at least 1 tpm in at least 0.000001 samples (we are basically
filtering out only 'all-zero' genes), and it has been corrected for gdc-platform.


#### Prepare Methylation

We have pre-processed methylation data for the following tumor types: BRCA, COAD, DLBC, KIRC, LAML, LIHC, PRAD, PAAD,
SKCM, STAD, LUAD, LUSC.

Configuration file (tcga-data-nf (0.0.13)): `conf/ methylation_prepare.conf`


### Analyze

We generated PANDA and LIONESS networks for 10 solid cancers: BRCA, COAD, KIRC, LIHC, LUAD, LUSC, PAAD, PRAD, SKCM,
STAD. 

We have used the prepared data with: 
- purity: 03
- normalization: logcpm
- gene filters: mintpm1, fracsamples01
- tissues: tissueall

## Authors

- Viola Fanfani, vfanfani@hsph.harvard.edu
