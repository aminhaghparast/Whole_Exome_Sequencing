

process REFERENCE_GENOME {
	publishDir "${params.outdir}/reference"
	label 'resource'
	
	output:
	file "ucsc.hg19.fasta"
	file "ucsc.hg19.dict"
	file "ucsc.hg19.fasta.fai"
	
    script:
	"""
	gunzip -dc /data/ucsc.hg19.fasta.gz > ucsc.hg19.fasta
	gunzip -dc /data/ucsc.hg19.dict.gz > ucsc.hg19.dict
	gunzip -dc /data/ucsc.hg19.fasta.fai.gz > ucsc.hg19.fasta.fai
	"""
}

process dbSNP {
	publishDir "${params.outdir}/reference"
	label 'resource'

	output:
	file "dbsnp_138.hg19.vcf" 
	file "dbsnp_138.hg19.vcf.idx"
    
    script:
	"""
	gunzip -dc /data/dbsnp_138.hg19.vcf.gz > dbsnp_138.hg19.vcf
	gunzip -dc /data/dbsnp_138.hg19.vcf.idx.gz > dbsnp_138.hg19.vcf.idx
	"""
}

process golden_indel {
	publishDir "${params.outdir}/reference"
	label 'resource'

	output:
	file "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
	file "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx"
    
    script:
	"""
	gunzip -dc /data/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz > Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
	gunzip -dc /data/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz > Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx
	"""
}

process hapmap {
	publishDir "${params.outdir}/reference"
	label 'resource'
    
	output:
	file "hapmap_3.3.hg19.sites.vcf"
	file "hapmap_3.3.hg19.sites.vcf.idx"

    script:
	"""
	gunzip -dc /data/hapmap_3.3.hg19.sites.vcf.gz > hapmap_3.3.hg19.sites.vcf
	gunzip -dc /data/hapmap_3.3.hg19.sites.vcf.idx.gz > hapmap_3.3.hg19.sites.vcf.idx
	"""
}

process omni {
	publishDir "${params.outdir}/reference"
	label 'resource'

	output:
	file "1000G_omni2.5.hg19.sites.vcf"
	file "1000G_omni2.5.hg19.sites.vcf.idx"

    script:
	"""
	gunzip -dc /data/1000G_omni2.5.hg19.sites.vcf.gz > 1000G_omni2.5.hg19.sites.vcf
	gunzip -dc /data/1000G_omni2.5.hg19.sites.vcf.idx.gz > 1000G_omni2.5.hg19.sites.vcf.idx
	"""
}

process phase1_SNPs {
	publishDir "${params.outdir}/reference"
	label 'resource'

	output:
	file "1000G_phase1.snps.high_confidence.hg19.sites.vcf"
	file "1000G_phase1.snps.high_confidence.hg19.sites.vcf.idx" 

    script:
	"""
	gunzip -dc /data/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz > 1000G_phase1.snps.high_confidence.hg19.sites.vcf
	gunzip -dc /data/1000G_phase1.snps.high_confidence.hg19.sites.vcf.idx.gz > 1000G_phase1.snps.high_confidence.hg19.sites.vcf.idx
	"""
}

process BWA_INDEX {
	publishDir "${params.outdir}/reference"
	label 'resource'

	output:
    tuple val("ucsc.hg19.fasta"), file ("ucsc.hg19.fasta.amb") , file ("ucsc.hg19.fasta.ann") , file ("ucsc.hg19.fasta.bwt") , file ("ucsc.hg19.fasta.pac")  , file ("ucsc.hg19.fasta.sa") 

    script:
	"""
	gunzip -dc /data/ucsc.hg19.fasta.amb.gz > ucsc.hg19.fasta.amb
	gunzip -dc /data/ucsc.hg19.fasta.ann.gz > ucsc.hg19.fasta.ann
	gunzip -dc /data/ucsc.hg19.fasta.bwt.gz > ucsc.hg19.fasta.bwt
	gunzip -dc /data/ucsc.hg19.fasta.pac.gz > ucsc.hg19.fasta.pac
	gunzip -dc /data/ucsc.hg19.fasta.sa.gz > ucsc.hg19.fasta.sa
	"""
}
