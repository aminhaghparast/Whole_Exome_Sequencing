# This Dockerfile has been used for creating aminhaghparast/annovar_annotation image in dockerhub ( 131 GB)
FROM ubuntu:20.04
RUN apt install wget
RUN apt-get install zip -y
RUN wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
RUN tar -xzvf annovar.latest.tar.gz && rm annovar.latest.tar.gz
RUN sudo mv annovar /opt/
RUN /opt/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene  /opt/annovar/humandb
RUN /opt/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gnomad211_genome /opt/annovar/humandb
RUN /opt/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03  /opt/annovar/humandb
RUN /opt/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2015aug  /opt/annovar/humandb
RUN /opt/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar gme  /opt/annovar/humandb 
RUN /opt/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar snp138  /opt/annovar/humandb
RUN /opt/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp150  /opt/annovar/humandb
RUN /opt/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20210501  /opt/annovar/humandb
RUN /opt/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar popfreq_all_20150413  /opt/annovar/humandb
RUN /opt/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp42a  /opt/annovar/humandb
RUN /opt/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar intervar_20180118 /opt/annovar/humandb
RUN /opt/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_all  /opt/annovar/humandb