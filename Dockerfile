# This Dockerfile has been used for creating aminhaghparast/wes_image image in dockerhub ( 3.94 GB)
FROM nfcore/base:1.14
RUN apt-get update
ENV DEBIAN_FRONTEND noninteractive
RUN apt install sudo -y
RUN apt install curl -y
RUN apt-get install zip -y
RUN mkdir data
RUN mkdir data/databases
COPY /bed_files/ /data/bed_files
COPY /environment/ /data/environment
COPY /main.nf .
RUN conda env create --quiet -f data/environment/environment.yml && conda clean -a
ENV PATH /opt/conda/envs/WES/bin:$PATH            #Add conda installation dir to PATH 
RUN mkdir Results
RUN conda env export --name WES > /Results/WES.yml



    


