#!/bin/bash
source s3put.sh 
numsec=3600 # 1 h 
sleep $numsec
while [ 1 ]
do
  if [ $SECONDS -ge $numsec ]
  then
    ora=$(date --date "3 hours ago" +%Y%m%d%H)
    ./datiGRADS.R $ora
    ./t2m19 $ora 1 datiGRADS.dat 3 2 ./ ./
    ./rhtd19 $ora 1 datiGRADS.dat 3 3 2 ./ ${ora:0:8}t2m_s.dat ${ora:0:8}t2m_g.dat ./
    ./plzln19 $ora 1 datiGRADS.dat 3 1 2 ./ ${ora:0:8}t2m_s.dat ./
    # upload in minio
    for dati in $ora_*.txt ; do 
      putS3 $dati "/" "analisi"
      if [[ "$?" != "0" ]] then
        echo "Errore nel caricare su MINIO il file " $dati
      else 
        "Caricato su MINIO il file " $dati
        rm -v $dati
      fi
    done
    # inserire update DBMETEO
    SECONDS=0
    sleep $numsec
  fi
done
