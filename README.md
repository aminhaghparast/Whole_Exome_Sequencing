# Whole-Exome-Sequencing





## Quickstart 

Install Nextflow by using the following command: 

    curl -s https://get.nextflow.io | bash 
    

Clone the repository by using the following command : 

    git clone https://github.com/aminhaghparast/Whole-Exome-Sequencing-.git

Install any of Docker or Singularity for full pipeline reproducibility.

test the pipeline on a minimal dataset with a single command : 

    nextflow run main.nf  -with-docker  aminhaghparast/wes_image:latest