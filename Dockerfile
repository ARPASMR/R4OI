FROM centos:centos6.10 as builder
RUN yum update -y
RUN yum install -y mysql-devel
RUN yum install -y gcc-4.4.7 gcc-gfortran-4.4.7
RUN yum install -y epel-release
RUN yum install -y R-3.5.2 
RUN R -e "install.packages('DBI', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('RMySQL', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('readBrukerFlexData', repos = 'http://cran.us.r-project.org')"
RUN yum install -y grads-2.0
RUN yum install -y centos-release-scl
RUN yum install -y python27
RUN yum install -y setup-python27
#RUN yum install -y phyton3
#RUN yum install -y setup-python3
RUN pip install rasdapy
#RUN pip3 install rasdapy
WORKDIR /usr/src/myapp
COPY * ./
RUN chmod a+x launcher.sh
RUN chmod a+x datiGRADS.R
FROM builder
#RUN make
CMD ["./launcher.sh"]
