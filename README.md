# R4OI

Extracts hourly data and station info from mysql DB and formats data as grads station file (1), runs OI F90 code (modified 2019) (2), updates rasdaman DB with output gridded analysis (3) and mysql DB with station point analysis (4).

1. R script to query mysql DB and wite station data files: datiGRADS.R
    - anagrafica.dat, anagrafica.ctl
    - datiGRADS.dat, datiGRADS.ctl (temperature,precipitazioni,umidità relativa)
    - pluvrisk.txt

2. Optimal interpolation 
    - 2.1 temperature: t2m19.f90
    - 2.2 relative humidity: rhtd19.f90
    - 2.3 precipitation: 
    
3. update rasdaman:   boh.py

4. update MySQL DB:    boh.R
