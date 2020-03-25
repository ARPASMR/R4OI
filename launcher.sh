#!/bin/bash
source s3put.sh 
numsec=600 # 1 h 
sleep $numsec
while [ 1 ]
do
  if [ $SECONDS -ge $numsec ]
  then
    ora=$(date --date "3 hours ago" +%Y%m%d%H)
    echo "Chiedo i dati per "$ora
    ./datiGRADS.R $ora
    ./t2m19 $ora 1 datiGRADS.dat 3 2 ./ ./
    ./rhtd19 $ora 1 datiGRADS.dat 3 3 2 ./ ${ora:0:8}t2m_s.dat ${ora:0:8}t2m_g.dat ./
    ./plzln19 $ora 1 datiGRADS.dat 3 1 2 ./ ${ora:0:8}t2m_s.dat ./
    echo "Finito interpolazione per ora: "$ora
    echo "Upload su minio..."
    # upload in minio
    for dati in $ora_*.txt ; do 
      echo "putS3 "." $dati "/" "analisi""
      putS3 "." $dati "" "analisi"
      if [[ "$?" != "0" ]] 
      then
        echo "Errore nel caricare su MINIO il file " $dati
      else 
        echo "Caricato su MINIO il file " $dati
        rm -v $dati
      fi
      rm -v temperatura_*.csv
      rm -v pluviometrizln_*.csv
      rm -v umidita_*.csv
    done
    # inserire update DBMETEO
    SECONDS=0
    sleep $numsec
  fi
done
