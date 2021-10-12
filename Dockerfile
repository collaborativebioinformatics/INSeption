#PREAMBLE
FROM alpine

WORKDIR /home/genomics
COPY . /home/genomics
RUN cd /home/genomics

RUN apk update \
	&& apk upgrade \
	&& apk add curl

ENV PATH /usr/local/bin:/home/genomics/scripts:$PATH
ENV R_LIB=/usr/local/lib/R
ENV R_LIBS_USER=/usr/local/lib/R

RUN mkdir -p ~/R_libs/


# MAIN
#Install base things
RUN apk add --no-cache --upgrade bash curl libressl-dev curl-dev libxml2-dev gcc g++ git coreutils ncurses linux-headers libgit2-dev libbz2

#Install base-R things
# Because R plotting requires the X11 to run, we have to find a way to build a headless install, we use xvfb for that. 
RUN apk add --no-cache --upgrade fontconfig-dev harfbuzz-dev fribidi-dev tiff tiff-dev libxt-dev cairo-dev \
  && apk add --no-cache --upgrade R R-dev R-doc \
  && apk add --no-cache --upgrade xorg-server-dev xorg-server libpng libpng-dev xvfb xvfb-run \ 
  && apk add --no-cache --upgrade libxfont-dev libxfont libfontenc font-xfree86-type1 \
  && apk add --no-cache --upgrade font-adobe-source-code-pro font-adobe-utopia-type1 font-adobe-100dpi font-adobe-75dpi font-adobe-utopia-100dpi font-adobe-utopia-75dpi \
  && apk add --no-cache --upgrade msttcorefonts-installer fontconfig \
  && update-ms-fonts \
  && R -e "install.packages('devtools',dependencies=TRUE, repos='http://cran.rstudio.com/')" \
  && R -e "install.packages('BiocManager', dependencies=TRUE, repos='http://cran.rstudio.com/'); BiocManager::install()" \
  && R -e "install.packages('ggplot2',dependencies=TRUE, repos='http://cran.rstudio.com/')" \
  && R -e "install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')"

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

#bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.12/bcftools-1.12.tar.bz2 \
  && tar jxf bcftools-1.12.tar.bz2 \
  && rm bcftools-1.12.tar.bz2 \
  && cd bcftools-1.12 \
  && ./configure --prefix /usr/local \
  && make \
  && make install \
  && cd ..

# Now we work on the Python side of the house 
FROM python:alpine

RUN pip install matplotlib pandas 

# Spades
ENV PATH /home/genomics/SPAdes-3.12.0-Linux/bin/:$PATH

RUN wget http://cab.spbu.ru/files/release3.12.0/SPAdes-3.12.0-Linux.tar.gz \
  && tar -xzf SPAdes-3.12.0-Linux.tar.gz \
  && cd SPAdes-3.12.0-Linux/bin/

# flye
RUN git clone https://github.com/fenderglass/Flye \
  && cd Flye \
  && make
