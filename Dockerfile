FROM centos:centos6.10
RUN yum update -y
RUN yum install -y gcc
WORKDIR /usr/src/myapp
COPY * .
RUN gfortran -static --verbose -L/usr/lib/x86_64-redhat-linux5E/lib64 t2m19.f90 subs19.f90 subsf89.f SpatialStuff.f90 \ 
              coords.f90 ll_gb.f90 ll_utm.f90 -o t2m19
RUN gfortran -static --verbose -L/usr/lib/x86_64-redhat-linux5E/lib64 rhtd19.f90  SpatialStuff.f90 subsf89.f subs19.f90 \
              coords.f90 ll_utm.f90 ll_gb.f90 -o rhtd19
RUN gfortran -static --verbose -L/usr/lib/x86_64-redhat-linux5E/lib64 plzln19.f90 SpatialStuff.f90 subs19.f90 subsf89.f \
              coords.f90 ll_utm.f90 ll_gb.f90 -o plzln19              
