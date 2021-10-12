#PREAMBLE
FROM alpine

WORKDIR /home/genomics
COPY . /home/genomics
RUN cd /home/genomics

RUN apk update \
	&& apk upgrade \
	&& apk add curl

ENV PATH /usr/local/bin:$PATH

#Install base things
RUN apk add --no-cache --upgrade bash curl libressl-dev curl-dev libxml2-dev gcc g++ git coreutils ncurses linux-headers libgit2-dev

#Install base-R things
RUN apk add --no-cache --upgrade fontconfig-dev harfbuzz-dev fribidi-dev tiff tiff-dev libxt-dev cairo-dev \
	&& apk add --no-cache --upgrade R R-dev R-doc \
	&& R -e "install.packages('devtools',dependencies=TRUE, repos='http://cran.rstudio.com/')" \
	&& R -e "install.packages('BiocManager', dependencies=TRUE, repos='http://cran.rstudio.com/'); BiocManager::install()"

#Install Aplications
RUN apk add build-base boost musl-dev make cmake zlib zlib-dev ncurses-dev boost-dev

#samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
  && tar jxf samtools-1.9.tar.bz2 \
  && rm samtools-1.9.tar.bz2  \
  && cd samtools-1.9 \
  && ./configure --prefix /usr/local \
  && make \
  && make install \
  && cd ..

RUN wget https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2 \
  && tar jxf bcftools-1.12.tar.bz2 \
  && rm bcftools-1.12.tar.bz2 \
  && cd bcftools-1.12 \
  && ./configure --prefix /usr/local \
  && make \
  && make install \
  && cd ..


#Install python from python alpine 
#FROM python:alpine

