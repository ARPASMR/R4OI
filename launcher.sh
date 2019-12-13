#!/bin/bash
numsec=3600 # 1 h 
sleep $numsec
while [ 1 ]
do
  if [ $SECONDS -ge $numsec ]
  then
    ora=$(date --date "1 hours ago" +%Y%m%d%H)
    ./datiGRADS.R $ora
    ./t2m19 $ora 1 datiGRADS.dat 3 2 ./ ./
    # inserire interpolazione rhtd
    # inserire interpolazione pluv
    # inserire upload rasdaman
    # inserire update DBMETEO
    SECONDS=0
    sleep $numsec
  fi
done
