FROM mskcc/roslin-variant-cmo-utils:1.9.15
LABEL maintainer="Nikhil Kumar (kumarn1@mskcc.org)" \
      version.image="1.0.0" \
      version.ngsfilters="1.4" \
      version.r="3.5.1" \
      version.alpine="3.8" \
      source.ngsfilters="https://github.com/mskcc/ngs-filters/releases/tag/v1.4" \
      source.r="https://pkgs.alpinelinux.org/package/edge/community/x86/R"

ENV NGS_FILTERS_VERSION 1.4
COPY runscript.sh /runscript.sh
COPY run_test.sh /run_test.sh
RUN sed -i 's/8/9/g' /etc/apk/repositories
RUN apk update
RUN apk upgrade
RUN sed -i 's/9/10/g' /etc/apk/repositories
RUN apk update
RUN apk upgrade
RUN sed -i 's/10/11/g' /etc/apk/repositories

RUN apk update
RUN apk del man py-pip py-setuptools python python-dev
RUN apk upgrade
RUN sed -i 's/11/12/g' /etc/apk/repositories

RUN apk update
RUN apk upgrade
RUN sed -i 's/12/13/g' /etc/apk/repositories

RUN apk update
RUN apk upgrade

RUN apk add python2 && ln -s /usr/bin/python2 python
RUN apk add py-pip
RUN pip install --no-cache --upgrade pip setuptools

RUN apk add --update \
      # for wget and bash
            && apk add ca-certificates openssl bash git \
      # for compilation (R packages)
            && apk add build-base musl-dev \
      # install R
            && apk add R R-dev 
      # download, unzip, install ngs-filters
RUN cd /tmp && wget https://github.com/mskcc/ngs-filters/archive/v${NGS_FILTERS_VERSION}.zip \
            && unzip v${NGS_FILTERS_VERSION}.zip \
      # install ngs-filters dependencies
            && cd /tmp/ngs-filters-${NGS_FILTERS_VERSION} \
            && Rscript --vanilla install-packages.R \
            && pip install -r requirements.txt \
            && apk add --update util-linux \
      # copy to /usr/bin/...
            && cp -r /tmp/ngs-filters-${NGS_FILTERS_VERSION} /usr/bin/ngs-filters/ \
      # clean up
            && rm -rf /var/cache/apk/* /tmp/* 
RUN chmod +x runscript.sh && \
	mv runscript.sh /usr/bin/ && \
	chmod +x run_test.sh && \
            exec /run_test.sh	

# disable per-user site-packages
ENV PYTHONNOUSERSITE set

RUN sed -i "s/getenv(None,/getenv('CMO_RESOURCE_CONFIG',/g" /usr/lib/python2.7/site-packages/cmo-1.9.15-py2.7.egg/cmo/util.py 
ENV CMO_RESOURCE_CONFIG=/genomic_resources.json
COPY genomic_resources.json genomic_resources.json


RUN sed -i 's|/opt/common/CentOS_6-dev/getbasecountsmultisample/v1.2.2/GetBaseCountsMultiSample|/usr/bin/GetBaseCountsMultiSample|g' /usr/bin/ngs-filters/maf_fillout.py
