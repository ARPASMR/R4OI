#!/bin/bash
source s3put.sh 
numsec=1800 # 1/2 h 
sleep $numsec
while [ 1 ]
do
  if [ $SECONDS -ge $numsec ]
  then
    ora=$(date --date "2 hours ago" +%Y%m%d%H)
    echo "Chiedo i dati per "$ora
    ./datiGRADS.R $ora
    ./t2m19 $ora 1 datiGRADS.dat 3 2 ./ ./  > t2m.log
    ./rhtd19 $ora 1 datiGRADS.dat 3 3 2 ./ ${ora:0:8}t2m_s.dat ${ora:0:8}t2m_g.dat ./  > rhtd.log
    ./plzln19 $ora 1 datiGRADS.dat 3 1 2 ./ ${ora:0:8}t2m_s.dat ./ > plzln.log
    echo "Finito interpolazione per ora: "$ora
    echo "########"
    echo "File presenti:  "
    ls -1 *
    echo "########"
    echo "Upload su minio..."
    # upload in minio
    for dati in *_32632_$ora.txt ; do
      echo "file per upload:  "$dati 
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
      rm ${ora:0:8}t2m_s.dat
      rm ${ora:0:8}t2m_g.dat
      rm ${ora:0:8}tdrh_s.dat
      rm ${ora:0:8}tdrh_g.dat
      rm ${ora:0:8}plzln_s.dat
      rm ${ora:0:8}plzln_g.dat
    done
    # inserire update DBMETEO
    SECONDS=0
    sleep $numsec
  fi
done
