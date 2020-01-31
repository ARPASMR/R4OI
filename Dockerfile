FROM centos:centos7 as builder
RUN yum update -y
RUN yum install -y mysql-devel
RUN yum install -y gcc gcc-gfortran
RUN yum install -y epel-release
RUN yum install -y R
RUN R -e "install.packages('DBI', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('RMySQL', repos = 'http://cran.us.r-project.org')"
RUN R -e "install.packages('readBrukerFlexData', repos = 'http://cran.us.r-project.org')"
RUN yum install -y grads
#RUN yum install -y centos-release-scl
#RUN yum install -y python27
RUN yum install -y python-pip
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
