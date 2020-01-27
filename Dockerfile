FROM centos:centos6.10 as builder
RUN yum update -y
RUN yum install -y gcc-4.4.7 gcc-gfortran-4.4.7
RUN yum install -y r-3.5.2
WORKDIR /usr/src/myapp
COPY * ./
RUN chmod a+x launcher.sh
RUN chmod a+x datiGRADS.R
FROM builder
RUN make
CMD ["./launcher.sh"]
