# List of all available data from AWS


## Raw data

Raw data is stored inside the
`https://tcga-data-nf-precomputed.s3.us-east-2.amazonaws.com/raw-data/firstround-20221102/` bucket. 

Data is accessible to anyone and we have a complete list of the files inside the `aws-raw-data-20241105.txt` file. 

**How do I download the methylation data for PAAD?**

1. Go inside `aws-raw-data-20241105.txt` and locate the file of interest, in this case: `tcga_paad/methylation/tcga_paad_methylation_manifest.txt`
2. Go on your browser
3. Paste the address
   `https://tcga-data-nf-precomputed.s3.us-east-2.amazonaws.com/raw-data/firstround-20221102/tcga_paad/methylation/tcga_paad_methylation_manifest.txt`
   where you can see is made of <bucket> + <file_path>.
4. Download the file. 

AWS also offers programmatic access and a command line interface. 



