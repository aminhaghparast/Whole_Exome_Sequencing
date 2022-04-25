process {
    withLabel: "process" { container 'aminhaghparast/wes_image:latest' }
    withLabel: "bwa" { container 'kathrinklee/bwa:latest' }
    withLabel: "gatk" { container 'broadinstitute/gatk:latest' }
    withLabel: "annotation" { container 'aminhaghparast/annovar_annotation:version1' }
    withLabel: "resource" { container 'oliversi/hg19' }
}

profiles {
    docker {
        docker.enabled = true
    }
}



