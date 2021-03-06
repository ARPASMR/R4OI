#! /usr/bin/env Rscript
#
# script per estrarre l'anagrafica in R e scrivere file stazione di Grads
# 
# MS, ottobre 2019
#

undef<- -999.9
main_dir<-"./"
stnmap<-"/usr/bin/stnmap"

file_log<-paste(main_dir,"prova.log",sep="")

file_dat_anagrafica<-paste(main_dir,"anagrafica.dat",sep="")
anag<-file(description=file_dat_anagrafica,open="wb")
file_ctl_anagrafica<-paste(main_dir,"anagrafica.ctl",sep="")
file_map_anagrafica<-paste(main_dir,"anagrafica.map",sep="")

file_dat_precipitazione<-paste(main_dir,"datiGRADS.dat",sep="")
preci<-file(description=file_dat_precipitazione,open="wb")
file_ctl_precipitazione<-paste(main_dir,"datiGRADS.ctl",sep="")
file_map_precipitazione<-paste(main_dir,"datiGRADS.map",sep="")
file_risk<-paste(main_dir,"pluvrisk.txt",sep="")

library(DBI)
library(RMySQL)
library(readBrukerFlexData) # serve una funzione della libreria per convertire double -> single


#________________________________________________________
#
cat ( "InIZio-=---=-- =---=-----=---==----iNiZio-----=\n" , file = file_log,append=FALSE)
#_______________________________________________________


# leggo la data richiesta
args <- commandArgs(trailing = TRUE)
if (length(args)!=1 | nchar(args[1],type="char")!=10 ) { 
	cat( "Manca la data o il formato è sbagliato! ", date()," \n\n" , file = file_log,append=FALSE)
	cat( "Usage: ./datiGRADS.R aaaammddhh " ," \n\n" , file = file_log,append=TRUE)
	quit(status=1)
} else { data_ora <- args[1] }

# controllo data
anno<-substr(data_ora,1,4)
mese<-substr(data_ora,5,6)
mesi_grads<-c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
mese_grads<-mesi_grads[as.integer(mese)]
giorno<-substr(data_ora,7,8)
ora<-substr(data_ora,9,10)

# scrivo la data per il DB meteo
data_DB <-paste("'",anno,"-",mese,"-",giorno," ",ora,"'",sep="")

cat ( "---------------------------------- \n" , file = file_log,append=T)
cat ( "ESTRAZIONE DATI DAL DB ",data_DB ," \n" , file = file_log,append=T)
cat ( "---------------------------------- \n" , file = file_log,append=T)
#___________________________________________________
#    COLLEGAMENTO AL DB
#___________________________________________________


cat("collegamento al DB...",file=file_log,append=T)
drv<-dbDriver("MySQL")
conn<-try(dbConnect(drv, user=as.character(Sys.getenv("MYSQL_USR")), password=as.character(Sys.getenv("MYSQL_PWD")), dbname=as.character(Sys.getenv("MYSQL_DBNAME")), host=as.character(Sys.getenv("MYSQL_HOST"))))

if (inherits(conn,"try-error")) {
  print( "ERRORE nell'apertura della connessione al DB \n")
  print( "chiusura connessione malriuscita ed uscita dal programma \n")
  dbDisconnect(conn)
  rm(conn)
  dbUnloadDriver(drv)
  quit(status=1)
}



#cat("collegamento al DB\n",file=file_log,append=T)
#MySQL(max.con=16,fetch.default.rec=500,force.reload=FALSE)

#definisco driver
#drv<-dbDriver("MySQL")

#apro connessione con il db descritto nei parametri del gruppo "tabella_rif"
#nel file "/home/meteo/.my.cnf
#conn<-dbConnect(drv,group="Visualizzazione_Sinergico")

#___________________________________________________
#    query anagrafica
#___________________________________________________

cat ( "---------------------------------- \n" , file = file_log,append=T)
cat ( "----------ANAGRAFICA-------------- \n" , file = file_log,append=T)
cat ( "---------------------------------- \n" , file = file_log,append=T)

# stazioni e tipologia

query<- paste("select a.IDstazione, b.IDsensore,NOMEtipologia from A_Stazioni as a, A_Sensori as b, A_Sensori2Destinazione as c  where a.IDstazione=b.IDstazione and b.IDsensore=c.IDsensore and Destinazione=11 and c.DataInizio < ",data_DB," and (c.DataFine is NULL OR c.DataFine > ",data_DB,") ",sep="")

q <- try(dbGetQuery(conn, query),silent=TRUE)
  if (inherits(q,"try-error")) {
    cat(q,"\n",file=file_log,append=T)
    quit(status=1)
  }

# numero stazioni con almeno una destinazione OI
query_num<- paste("select a.IDstazione, count(a.IDstazione) as sensori,Y(a.CoordUTM),X(a.CoordUTM),Quota,UrbanWeight,IDrete from A_Stazioni as a, A_Sensori as b, A_Sensori2Destinazione as c  where a.IDstazione=b.IDstazione and b.IDsensore=c.IDsensore and Destinazione=11 and c.DataInizio < ",data_DB," and (c.DataFine is NULL OR c.DataFine > ",data_DB,") group by a.IDstazione",sep="")

num <- try(dbGetQuery(conn, query_num),silent=TRUE)
  if (inherits(num,"try-error")) {
    cat(num,"\n",file=file_log,append=T)
    quit(status=1)
  }

num_staz <- length(num[,1])

cat(paste("Numero di stazioni con almeno una destinazione OI: ",num_staz,"\n \n", sep=""),file = file_log,append=T)

# scrivo il file dati binario ed il ctl
# NOTE:-----------------------------------
# R memorizza tutti i reali come doppia precisione (8 byte), ci vuole una specifica libreria per cnvertirli in singola
# R memorizza gli interi come 4 byte, ma sono numeric: bisogna specificare as.integer prima di scrivere!
# ----------------------------------------


for(i in 1:num_staz) {
# HEADER-----------------------------------
	writeBin(paste("Lo",sprintf("%05d", num[i,1]),sep=""),anag,size=8)
	writeBin(readBrukerFlexData:::.double2singlePrecision(num[i,3]),anag,size=4)
        writeBin(readBrukerFlexData:::.double2singlePrecision(num[i,4]),anag,size=4)
# t - Time in grid-relative units  
        writeBin(readBrukerFlexData:::.double2singlePrecision(0.0),anag,size=4)
# nlev - Number of data groups following the header: 1 surface group + vert levls
#        writeBin(readBrukerFlexData:::.double2singlePrecision(1),anag,size=4) questa non funziona!!
        writeBin(as.integer(1),anag,size=4) # questa si!!
# flag - If set to 1, then there are surface variables following the header.
        writeBin(as.integer(1),anag,size=4)
# Variabili---------------------------------
# Quota
	writeBin(readBrukerFlexData:::.double2singlePrecision(num[i,5]),anag,size=4)
# Indice urbanita
        if ( is.na(num[i,6]) ) { 
		writeBin(readBrukerFlexData:::.double2singlePrecision(undef),anag,size=4) 
	} else { 
		writeBin(readBrukerFlexData:::.double2singlePrecision(num[i,6]),anag,size=4) 
	}
# ID rete : 1 - RRQA, 2 - CMG, 4 - INM, 5 - altro fuori LO, 6 - altro in Lo 
	writeBin(readBrukerFlexData:::.double2singlePrecision(num[i,7]),anag,size=4)
# ID Termometro
	if ( length(q[which(q[,3]=="T" & q[,1] == num[i,1]),2])==0 ) { 
		writeBin(readBrukerFlexData:::.double2singlePrecision(undef),anag,size=4)
	} else if (length(q[which(q[,3]=="T" & q[,1] == num[i,1]),2])>1 )  {
		cat(paste(" Trovato sensore doppio nella stazione: ID ",num[i,1],"\n",sep=""), file=file_log, append=T)
		quit(status=1)
	} else { 
		writeBin(readBrukerFlexData:::.double2singlePrecision(q[which(q[,3]=="T" & q[,1] == num[i,1]),2]),anag,size=4)
	}
#	ID Igrometro	
	if ( length(q[which(q[,3]=="UR" & q[,1] == num[i,1]),2])==0 ) { 
		writeBin(readBrukerFlexData:::.double2singlePrecision(undef),anag,size=4)
	} else if (length(q[which(q[,3]=="UR" & q[,1] == num[i,1]),2])>1 )  {
		cat(paste(" Trovato sensore doppio nella stazione: ID ",num[i,1],"\n",sep=""), file=file_log, append=T)
		quit(status=1)
	} else { 
		writeBin(readBrukerFlexData:::.double2singlePrecision(q[which(q[,3]=="UR" & q[,1] == num[i,1]),2]),anag,size=4)
	}
#       ID pluviometro
	if ( length(q[which(q[,3]=="PP" & q[,1] == num[i,1]),2])==0 ) { 
		writeBin(readBrukerFlexData:::.double2singlePrecision(undef),anag,size=4)
	} else if (length(q[which(q[,3]=="PP" & q[,1] == num[i,1]),2])>1 )  {
		cat(paste(" Trovato sensore doppio nella stazione: ID ",num[i,1],"\n",sep=""), file=file_log, append=T)
		quit(status=1)
	} else { 
		writeBin(readBrukerFlexData:::.double2singlePrecision(q[which(q[,3]=="PP" & q[,1] == num[i,1]),2]),anag,size=4)
	}
}
# time group terminator - nell'anagrafica ho un solo istante, quindi finisco il file...
writeBin(paste("Lo",sprintf("%05d", num[i,1]),sep=""),anag,size=8)
writeBin(readBrukerFlexData:::.double2singlePrecision(num[num_staz,3]),anag,size=4)
writeBin(readBrukerFlexData:::.double2singlePrecision(num[num_staz,4]),anag,size=4)
writeBin(readBrukerFlexData:::.double2singlePrecision(0.0),anag,size=4)
writeBin(as.integer(0),anag,size=4) 
writeBin(as.integer(1),anag,size=4)
cat( "File anagrafica.dat scritto, chiudo il file" ," \n\n" , file = file_log,append=TRUE)
close(anag)

# scrivo il *.ctl 
cat(paste("DSET   ",file_dat_anagrafica,"
DTYPE  station 
STNMAP ",file_map_anagrafica,"
UNDEF  -999.9
TITLE  Anagrafica Stazioni 
TDEF   1 linear ",ora,"Z",giorno,mese_grads,anno," 1hr
VARS 6
elev  0  99  Quota stazione sul livello del mare
whi   0  99  Indice di urbanita [0-1]
rete  0  99  ID rete: 1 - RRQA, 2 - CMG, 4 - INM
temp  0  99  ID sensore T 
ur    0  99  ID sensore UR 
prec  0  99  ID sensore PP
ENDVARS",sep=""), file=file_ctl_anagrafica)

# genero il file di mappa
command<-paste(stnmap," -i ",file_ctl_anagrafica,sep="")
#command<-paste(stnmap," -i ",file_ctl_anagrafica, " >> ", file_log,sep="")
output<-system(command)
if (output!=0) {
    cat("Errore nel generare il file stn map per Grads ","\n",file=file_log,append=T)
    quit(status=1)
}
cat( "File anagrafica.ctl e file anagrafica.map scritti correttamente " ," \n\n" , file = file_log,append=TRUE)


#___________________________________________________
#
#    query temperature
#__________________________________________________

cat ( "---------------------------------- \n" , file = file_log,append=T)
cat ( "----------TEMPERATURE------------- \n" , file = file_log,append=T)
cat ( "---------------------------------- \n" , file = file_log,append=T)
cat( paste("Chiedo al DB METEO le misure dei termometri per l'ora  ",data_DB," \n\n",sep="") , file = file_log,append=TRUE)

id_termometri <-q[which(q[,3]=="T"),2]
lista_termometri <- paste(id_termometri,collapse=",")
query_temp<-paste("select staz.IDstazione, Misura from M_Termometri_",anno," as mis, A_Stazioni as staz, A_Sensori as sens where staz.IDstazione=sens.IDstazione and sens.IDsensore=mis.IDsensore and flag_manuale!='E' and flag_automatica!='F' and Data_e_ora=",data_DB," and sens.IDsensore in (",lista_termometri,")  and mis.IDsensore not in (select IDsensore from A_ListaNera where (DataInizio< ",data_DB," and DataFine is NULL) or (DataInizio< ",data_DB," and DataFine>",data_DB,"));",sep="")
temp <- try(dbGetQuery(conn, query_temp),silent=TRUE)
  if (inherits(temp,"try-error")) {
    cat(temp,"\n",file=file_log,append=T)
    quit(status=1)
  }


#___________________________________________________
#
#    query umidita relativa
#__________________________________________________

cat ( "---------------------------------- \n" , file = file_log,append=T)
cat ( "----------UMIDITA REL------------- \n" , file = file_log,append=T)
cat ( "---------------------------------- \n" , file = file_log,append=T)
cat( paste("Chiedo al DB METEO le misure degli igrometri per l'ora  ",data_DB," \n\n",sep="") , file = file_log,append=TRUE)

id_igrometri <-q[which(q[,3]=="UR"),2]
lista_igrometri <- paste(id_igrometri,collapse=",")
query_urel<-paste("select staz.IDstazione, Misura from M_Igrometri_",anno," as mis, A_Stazioni as staz, A_Sensori as sens where staz.IDstazione=sens.IDstazione and sens.IDsensore=mis.IDsensore and flag_manuale!='E' and flag_automatica!='F' and Data_e_ora=",data_DB," and sens.IDsensore in (",lista_igrometri,") and mis.IDsensore not in (select IDsensore from A_ListaNera where (DataInizio< ",data_DB," and DataFine is NULL) or (DataInizio< ",data_DB," and DataFine>",data_DB,")); ",sep="")
urel <- try(dbGetQuery(conn, query_urel),silent=TRUE)
  if (inherits(urel,"try-error")) {
    cat(urel,"\n",file=file_log,append=T)
    quit(status=1)
  }


#___________________________________________________
#
#    query precipitazioni
#__________________________________________________

cat ( "---------------------------------- \n" , file = file_log,append=T)
cat ( "----------PRECIPITAZIONI---------- \n" , file = file_log,append=T)
cat ( "---------------------------------- \n" , file = file_log,append=T)
cat( paste("Chiedo al DB METEO le misure dei pluviometri per l'ora  ",data_DB," \n\n",sep="") , file = file_log,append=TRUE)

id_pluviometri <-q[which(q[,3]=="PP"),2]
lista_pluviometri <- paste(id_pluviometri,collapse=",")
query_pluv<-paste("select staz.IDstazione, Misura from M_Pluviometri_",anno," as mis, A_Stazioni as staz, A_Sensori as sens where staz.IDstazione=sens.IDstazione and sens.IDsensore=mis.IDsensore and flag_manuale!='E' and flag_automatica!='F' and Data_e_ora=",data_DB," and sens.IDsensore in (",lista_pluviometri,") and mis.IDsensore not in (select IDsensore from A_ListaNera where (DataInizio< ",data_DB," and DataFine is NULL) or (DataInizio< ",data_DB," and DataFine>",data_DB,")); " ,sep="")
pluv <- try(dbGetQuery(conn, query_pluv),silent=TRUE)
  if (inherits(pluv,"try-error")) {
    cat(temp,"\n",file=file_log,append=T)
    quit(status=1)
  }

#_____________________________________________________
#
#scrivo il file grads su tutte le stazioni con almeno una destinazione OI
#
#_____________________________________________________

cat ( "---------------------------------- \n" , file = file_log,append=T)
cat ( "----Scrittura file GRADS---------- \n" , file = file_log,append=T)
cat ( "---------------------------------- \n" , file = file_log,append=T)


for(i in 1:num_staz) {
# HEADER-----------------------------------
	writeBin(paste("Lo",sprintf("%05d", num[i,1]),sep=""),preci,size=8)
	writeBin(readBrukerFlexData:::.double2singlePrecision(num[i,3]),preci,size=4)
        writeBin(readBrukerFlexData:::.double2singlePrecision(num[i,4]),preci,size=4)
# t - Time in grid-relative units  
        writeBin(readBrukerFlexData:::.double2singlePrecision(0.0),preci,size=4)
# nlev - Number of data groups following the header: 1 surface group + vert levls
#        writeBin(readBrukerFlexData:::.double2singlePrecision(1),anag,size=4) questa non funziona!!
        writeBin(as.integer(1),preci,size=4) # questa si!!
# flag - If set to 1, then there are surface variables following the header.
        writeBin(as.integer(1),preci,size=4)
# Variabili---------------------------------
#       PP
        if ( length(pluv[which(pluv[,1] == num[i,1]),1])==0 ) { 
		writeBin( readBrukerFlexData:::.double2singlePrecision(undef),preci,size=4)
        } else {
                writeBin( readBrukerFlexData:::.double2singlePrecision(pluv[which(pluv[,1]== num[i,1]),2]),preci,size=4)
        }
#       T
        if ( length(temp[which(temp[,1] == num[i,1]),1])==0 ) {
                writeBin( readBrukerFlexData:::.double2singlePrecision(undef),preci,size=4)
        } else {
                writeBin( readBrukerFlexData:::.double2singlePrecision(temp[which(temp[,1]== num[i,1]),2]),preci,size=4)
        }
#       UR
        if ( length(urel[which(urel[,1] == num[i,1]),1])==0 ) {
                writeBin( readBrukerFlexData:::.double2singlePrecision(undef),preci,size=4)
        } else {
                writeBin( readBrukerFlexData:::.double2singlePrecision(urel[which(urel[,1]== num[i,1]),2]),preci,size=4)
        }
}
# time group terminator - nell'anagrafica ho un solo istante, quindi finisco il file...
writeBin(paste("Lo",sprintf("%05d", num[i,1]),sep=""),preci,size=8)
writeBin(readBrukerFlexData:::.double2singlePrecision(num[num_staz,3]),preci,size=4)
writeBin(readBrukerFlexData:::.double2singlePrecision(num[num_staz,4]),preci,size=4)
writeBin(readBrukerFlexData:::.double2singlePrecision(0.0),preci,size=4)
writeBin(as.integer(0),preci,size=4) 
writeBin(as.integer(1),preci,size=4)
cat( "File datiGRADS.dat scritto, chiudo il file" ," \n\n" , file = file_log,append=TRUE)
close(preci)

#______________________________________________
#
# scrivo il *.ctl
#______________________________________________

cat ( "---------------------------------- \n" , file = file_log,append=T)
cat ( "--Scrittura del file *.ctl-------- \n" , file = file_log,append=T)
cat ( "---------------------------------- \n" , file = file_log,append=T)

cat(paste("DSET   ",file_dat_precipitazione,"
DTYPE  station
STNMAP ",file_map_precipitazione,"
UNDEF  -999.9
TITLE  Dati precipitazione
TDEF   1 linear ",ora,"Z",giorno,mese_grads,anno," 1hr
VARS 3 
prec  0  99  precipitazione
temp  0  99  temperatura
urel  0  99  umidita relativa
ENDVARS",sep=""), file=file_ctl_precipitazione)

#______________________________________________
#
# genero il file di mappa
#_____________________________________________



command<-paste(stnmap," -i ",file_ctl_precipitazione,sep="")
#command<-paste(stnmap," -i ",file_ctl_precipitazione, " >> ", file_log,sep="")
output<-system(command)
if (output!=0) {
    cat("Errore nel generare il file stn map per Grads ","\n",file=file_log,append=T)
    quit(status=1)
}
cat( "File datiGRADS.ctl e file datiGRADS.map scritti correttamente " ," \n\n" , file = file_log,append=TRUE)


#______________________________________________________
#
# genero un file 'pluvrisk.txt' con elenco dei pluviometri riscaldati
#______________________________________________________

query_risk <- paste( "select a.IDstazione from A_Sensori as a, A_Sensori_specifiche as b where a.IDsensore=b.IDsensore and a.IDsensore in (",lista_pluviometri,") and b.RiscVent='Yes' and (DataDisistallazione is NULL or DataDisistallazione='0000-00-00') ; " ,sep="")
#print(query_risk)
risk <- try(dbGetQuery(conn, query_risk),silent=TRUE)
  if (inherits(risk,"try-error")) {
    cat(temp,"\n",file=file_log,append=T)
    quit(status=1)
  }

write(risk$IDstazione,file_risk,ncolumns=1)


#___________________________________________________
#    DISCONNESSIONE DAL DB
#___________________________________________________

# chiudo db
cat ( "chiudo DB \n" , file = file_log , append = TRUE )
RetCode<-try(dbDisconnect(conn),silent=TRUE)
if (inherits(RetCode,"try-error")) {
  cat(RetCode,"\n",file=file_log,append=T)
  quit(status=1)
}

rm(conn)
dbUnloadDriver(drv)


cat ( "PROGRAMMA ESEGUITO CON SUCCESSO alle ", date()," \n" , file = file_log , append = TRUE )
cat ( "FIne --=---=-- getcsv_recenti.R =---=-----=---==-------fInE-=\n" , file = file_log,append=T)

quit(status=0)


