FROM debian:stretch-slim
RUN apt-get update
RUN apt-get install -y gcc
WORKDIR /usr/src/myapp
COPY *.f90 .
COPY *.f .
COPY *.sh .
COPY oro* .
