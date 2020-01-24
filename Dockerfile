FROM centos:centos6.10 as builder
RUN yum update -y
RUN yum install -y gcc gcc-gfortran
WORKDIR /usr/src/myapp
COPY * ./
RUN chmod a+x launcher.sh
RUN chmod a+x datiGRADS.R
FROM builder
RUN make
CMD ["./launcher.sh"]
