#!/bin/bash
numsec=1800 # 1/2 h 
sleep $numsec
while [ 1 ]
do
  if [ $SECONDS -ge $numsec ]
  then
    ora=$(date --date "2 hours ago" +%Y%m%d%H)
    echo "Chiedo i dati per "$ora
    ./datiGRADS.R $ora
    if [[ "$?" -gt "0" ]]
    then
      logger --id --stderr --server $SYSLOG_MASTER -P $SYSLOG_PORT -T -p user.err -t OSSERVAZIONI_ELABORAZIONE "TAGMETEO R4OI errore nella richiesta DBMETEO o formattazione dati per $ora"
      echo "Errore nella richiesta DBMETEO o formattazione dati per $ora"
    fi
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
      s3cmd --config=config_minio.txt put $dati s3://analisi
      if [[ "$?" -gt "1" ]] 
      then
        echo "Errore nel caricare su MINIO il file " $dati
        logger --id --stderr --server $SYSLOG_MASTER -P $SYSLOG_PORT -T -p user.warning -t OSSERVAZIONI_ELABORAZIONE "TAGMETEO R4OI errore upload su minio del file $dati"
      else 
        echo "Caricato su MINIO il file " $dati
        rm -v $dati
      fi
    done
    rm -v temperatura_*.csv
    rm -v pluviometrizln_*.csv
    rm -v umidita_*.csv
    rm -v ${ora:0:8}t2m_s.dat
    rm -v ${ora:0:8}t2m_g.dat
    rm -v ${ora:0:8}tdrh_s.dat
    rm -v ${ora:0:8}tdrh_g.dat
    rm -v ${ora:0:8}plzln_s.dat
    rm -v ${ora:0:8}plzln_g.dat
    rm -v RMS.txt 
    # inserire update DBMETEO
    SECONDS=0
    sleep $numsec
  fi
done
