# Supplememtary code for the paper "Reproducible processing of TCGA regulatory networks"


In this repo we host the code to generate the data and figures for the paper
"Reproducible processing of TCGA regulatory networks". 


```
.
├── LICENSE
├── README.md
├── config
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


## Raw data for 10 common cancers

[How to download data guide](https://github.com/QuackenbushLab/tcga-data-supplement/data/manifests/manifests.md)



