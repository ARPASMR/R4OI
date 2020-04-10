program plzln19
!-------------------------------------------------------------------------------
! CODICE ORIGINALE:
! analisi di precipitazione da pluviometri - ARPA Lombardia
!
! $ plzln11 yyyymmddhh ninst <filedati> nvar nposPLUV nposTEMP <ana_oro_dir> \
!           <file_analisi_T_stazioni> <output_dir>
!
! orography grid Lombardia GB 1500m 177x174
!
! Francesco Uboldi, 2011
!-------------------------------------------------------------------------------
! CSV contiene anche il residuo sqrt((yo-ya)(yo-yaV)) , ma senza fare SCT (senza loop QC)
! ao         :   UNDEF           0        >0
! pluviometro:   NO            DRY        WET
! idiwet     :   idiw         idiw       idiw   ! i valori IDI, WET e DRY sono
! ididry     :   idid         idid       idid   ! calcolati per tutte le stazioni
! aa         :   (idiw/d)        0       awet   ! se ao=UNDEF aa dipende dagli IDI, WET e DRY
!
! CV
! idiwetV    :  =idiw        =idiw      idiwV   ! per pluviometro DRY: idiwV=idiw
! ididryV    :  =idid        ididV      =idid   ! per pluviometro WET: ididV=idid
! aaV        : =aa(idiVw/d) (idiVw/d) (idiVw/d) ! aaV in ogni caso dipende dai CV-IDI, WET e DRY
!
! aa(idiwV,ididV) = 0 se idid*0.6 >= idiw
!                 = aa-WET altrimenti
!
! aav(idiwV,ididV) =  0 se CV-IDI-DRY*0.6 >= CV-IDI-WET
!                  = aaV-WET altrimenti
!
! quindi nel CSV ho anche idiV-wet e idiV-dry:
! Id; Year/Mo/Dy Hr;      yo;       ya;      yaV;      res;  idi-wet;  idi-dry; idiV-wet; idiV-dry;
!
! dei 4 idi, 2 sono sempre uguali:
!   o idi-wet=idiV-wet (se la stazione e' DRY)
!   o idi-dry=idiV-dry (se la stazione e' WET)
!-------------------------------------------------------------------------------
!CODICE MODIFICATO
!
! Modificate routine lettura e scrittura per lavorare su griglia UTM 251x249
! con coord. stazione in UTM 
! Non produce file cumulate giornaliere
!
! Marta, dicembre 2019
!-------------------------------------------------------------------------------
! calcolo IDI da pluviometri con misura 0.D0--> XIDIdry
! poi uso solo pluviometri com misura positiva WET
! converto rain (mm/h) in refl
! interpolo refl
! riconverto refl in rain
! assegno rain solo se XIDI > XIDIdry, altrimenti zero
! no decorrelazioni verticali
!-------------------------------------------------------------------------------
 use SpatialStuff
 implicit none
!-------------------------------------------------------------------------------
 integer, parameter :: KTOT=1000, NTOT=24 & ! max 1000 stazioni e 1 giorno
&                    , INX=251, INY=249, IM=INX*INY, IEFF=IM  ! non conosco IEFF, uso il massimo valore

 real(8), parameter :: UNDEF=-9999.d0

! at least KKMIN obs each time
 integer, parameter :: KKMIN= 20

!-------------------------------------------------------------------------------
! grid
 real(4), dimension(IM) :: g4x, g4y, g4hgt,g4hui,g4msk,g4lon,g4lat,gxa,gidiw,gidid!, gxacum24h

 real(8), dimension(IEFF) :: xa, xlon, xlat, xhgt, dx!, xacum24h

 real(8) :: xum, arg, rah, ggh
 integer :: i, II, ig, IIG

!-----------------------------------------------------------------------
! i-th line of matrix G:
 real(8), dimension(IEFF,KTOT) :: GG
 integer :: norec, ntrec, nwrec
! logical itisthere

!-----------------------------------------------------------------------
! station archive
 character(8), dimension(KTOT) :: atid
 real(8), dimension(KTOT) :: ax, ay, alon, alat, ahgt &
&                          , ao, aa, aav, az, ady &
&                          , aidid, aidiw, aididv, aidiwv &
&                          , adiald, adialw !, anvdiasrnv

 integer, dimension(KTOT) :: jarkwet, jarkdry
 integer :: KKA

!-----------------------------------------------------------------------
! DRY, WET
 integer :: KKWET, KKDRY, IW
 real(8), dimension(KTOT) :: zidid, didid
 real(8), dimension(KTOT) :: yd, zidiw, didiw
 real(8) :: ydmax

 real(8), dimension(IM) :: xidiw, xidid

!-----------------------------------------------------------------------
! functions:
 real(8) :: rainrefl, reflrain

!-----------------------------------------------------------------------
! observations R4 variables (GrADS format I/O)
 real(4), dimension(KTOT) :: so, sa, sav, s4x, s4y, shgt
 real(4), dimension(KTOT) :: tsta
 character(8), dimension (KTOT) :: stid, tstid
 integer :: KK, KKT, nvars, lpospluv, lpostemp

!-----------------------------------------------------------------------
! DQC SCT residual:
 real(8) :: absd
 real(8), dimension(KTOT) :: asctres

 character(200) :: arow

!-----------------------------------------------------------------------
! oi07 specific variables:
 real(8) :: ssh
 integer :: k2, j2

! previous index array:
 integer, dimension(KTOT) :: jarkwp, jarkdp
 integer :: KKDP, KKWP

! station-station b.f. error covariance matrix
 real(8), dimension(KTOT,KTOT) :: ZS, SRNVD, SRNVDP, SRNVW, SRNVWP
 real(8), dimension(KTOT) :: diald, dialw
 real(8) :: dialdmx, dialwmx

!-----------------------------------------------------------------------
! analysis parameters and auxiliary variables
 real(8), parameter :: eps= 0.07d0, sigd= 10.d0

!-----------------------------------------------------------------------
 real(8) :: zum, res, qres, pres, resw, qresw, presw, cvscwet, cvscdry

!-----------------------------------------------------------------------
! files:
 character(200) :: STAGFI, OROGRIFI, PRSKFI, STANFI, GRIDFI, DATDIR, OUTDIR &
&   , FULFISTAG, FULFIORO, FULFIOB, FULFIST, FULFIGR, CSVFI, FULFICSV &
&   , FULFICUM, CUMFI, FULFITST, FULFIPRSK, FULFIASCIIa, FULFIASCIId &
&   , FULFIASCIIw

!-----------------------------------------------------------------------
 character(10) :: adate
 character(12) :: datehour
 character(12), dimension(NTOT) :: datehr
 integer :: NN
 integer :: lhour, lmin, lyear, lmonth, lday, ldayfeb
 integer :: lodir, ladir

!-----------------------------------------------------------------------
 character(6) :: ann
 character(2) :: anvrs, apospluv, apostemp

 logical :: griglia

!-----------------------------------------------------------------------
! functions:
 real(8) :: distan
 integer :: lstrim

!-----------------------------------------------------------------------
 integer :: j, k, m, n

!-----------------------------------------------------------------------
 ! character(1) :: ans

!===================================================================================================
 adate= 'yyyymmddhh'
 datehour= 'yearmodyhrmn'
 FULFIOB= ''
 OUTDIR= ''
 DATDIR= './'
 NN= 24
 griglia= .true.

!-----------------------------------------------------------------------
 call getarg(1,adate)
 if (adate(1:1) == ' ') then
  write (6,*) 'plzln19 yyyymmddhh ninst <filedati> nvar nposPLUV nposTEMP' &
& //' <ana_oro_dir> <file_analisi_T_stazioni> <output_dir>'
  write (6,*) ' '
  write (6,*) ' OPPURE, per stimare solo su stazioni, "-" prima di ninst:'
  write (6,*) 'plzln19 yyyymmddhh -ninst <filedati> nvar nposPLUV nposTEMP' &
& //' <ana_oro_dir> <file_analisi_T_stazioni> <output_dir>'
  stop 1
 end if

 call getarg(2,ann)
 read (ann,*) NN
 if (NN < 0) then
  griglia= .false.
  NN= -NN
  write (ann,'(i6)') NN
  ann= adjustl(ann)
  write (6,*) ' '
  write (6,*) 'stima solo su stazioni!'
  write (6,*) ' '
  call flush(6)
 end if 

 call getarg(3,FULFIOB)
 call getarg(4,anvrs)
 call getarg(5,apospluv)
 call getarg(6,apostemp)
 call getarg(7,DATDIR)
 call getarg(8,FULFITST)
 call getarg(9,OUTDIR)

 read (anvrs,*) nvars
 read (apospluv,*) lpospluv
 read (apostemp,*) lpostemp

!-----------------------------------------------------------------------
 lodir= lstrim(OUTDIR)
 if (OUTDIR(lodir:lodir) /= '/') then
  lodir= lodir +1
  OUTDIR(lodir:lodir)='/'
 end if
 write (6,*) 'output directory: ', OUTDIR(1:lodir)

 if (DATDIR(1:1) == ' ') DATDIR= './'
 ladir= lstrim(DATDIR)
 if (DATDIR(ladir:ladir) /= '/') then
  ladir= ladir +1
  DATDIR(ladir:ladir)='/'
 end if
 write (6,*) 'anag-orog-matrix directory: ', DATDIR(1:ladir)
 write (6,*) ' '

!-----------------------------------------------------------------------
! costruzione nomi file:
 OROGRIFI= 'oromaskutm_mauri.dat'
 STAGFI=   'anagrafica.dat'
 PRSKFI=   'pluvrisk.txt'

 STANFI= adate(1:8) //'plzln_s.dat'
 GRIDFI= adate(1:8) //'plzln_g.dat'
! CUMFI= adate(1:8) // 'CUMplzln_g.dat'

 CSVFI= 'pluviometrizln_' //adate //'_' //trim(ann) //'.csv'

!-----------------------------------------------------------------------
! DIRECTORY:
! DATDIR:
 FULFIPRSK= DATDIR(1:ladir) //PRSKFI
 FULFISTAG= DATDIR(1:ladir) //STAGFI
 FULFIORO = DATDIR(1:ladir) //OROGRIFI

! OUTDIR:
 FULFIGR= OUTDIR(1:lodir) //GRIDFI
 FULFIST= OUTDIR(1:lodir) //STANFI
! FULFICUM= OUTDIR(1:lodir) //CUMFI
 FULFICSV= OUTDIR(1:lodir) //CSVFI

!-----------------------------------------------------------------------
 if (griglia) write (6,*) 'orography and grid file    : ', trim(FULFIORO)
 write (6,*) 'station archive file       : ', trim(FULFISTAG)
 write (6,*) 'observation file           : ', trim(FULFIOB)
 write (6,*) 'T station analysis file    : ', trim(FULFITST)
 write (6,*) 'station analysis file      : ', trim(FULFIST)
 if (griglia) write (6,*) 'grid analysis file - binary  : ', trim(FULFIGR)
 if (griglia) write (6,*) 'grid analysis file - ascii   : ', trim(FULFIASCIIa)
! if (griglia) write (6,*) 'grid cum-24h analysis file : ', trim(FULFICUM)
 write (6,*) 'CSV station file           : ', trim(FULFICSV)
 write (6,*) ' '

!-----------------------------------------------------------------------
 write (6,'(a,3f10.2)') 'eps, sigd ', eps, sigd
 write (6,*) ' ok?'
 ! read (5,'(a1)') ans
 ! write (6,*) ' '
 call flush(6)

!-----------------------------------------------------------------------
! data iniziale e durata evento
 read (adate(1:4),'(i4)') lyear
 read (adate(5:6),'(i2)') lmonth
 read (adate(7:8),'(i2)') lday
 read (adate(9:10),'(i2)') lhour
 lmin= 0
 datehour= adate //'00'
 write (6,*) 'data ora iniziale: ' //datehour(1:4) //'/'  //datehour(5:6) &
& //'/'  //datehour(7:8) //' '  //datehour(9:10) //'h' //datehour(11:12) 
 write (6,*) 'N. istanti: ', NN

!-----------------------------------------------------------------------
 ldayfeb= 28
 if (mod(lyear,4) == 0) ldayfeb= 29
 if (mod(lyear,100) == 0) ldayfeb= 28
 if (mod(lyear,400) == 0) ldayfeb= 29

!-----------------------------------------------------------------------
! prima preparo la stringa con le date
 datehr= ''
 n= 1
 datehr(n)= datehour
 write (6,'(a,i3,1x,a12)') 'n datehr: ', n, datehr(n)

 do n=2,NN
  read (datehr(n-1)( 1: 4),'(i4)') lyear
  read (datehr(n-1)( 5: 6),'(i2)') lmonth
  read (datehr(n-1)( 7: 8),'(i2)') lday
  read (datehr(n-1)( 9:10),'(i2)') lhour
  read (datehr(n-1)(11:12),'(i2)') lmin
  lhour= lhour +1

! gestione del cambio data 
  if (lhour == 24) then
   lhour= 0
   if (lday == 31 .and.lmonth == 12) then
    lday= 1
    lmonth= 1
    lyear= lyear +1
   else if ((lmonth == 2. .and. lday == ldayfeb) .or. lday == 31 .or. &
&           (lday == 30 .and. &
&           (lmonth == 11 .or. lmonth == 4 .or. lmonth == 6 .or. lmonth == 9)) &
&          ) then
    lday= 1
    lmonth= lmonth +1
   else
    lday= lday +1
   end if
  end if

  write (datehr(n)( 1: 4),'(i4.4)') lyear
  write (datehr(n)( 5: 6),'(i2.2)') lmonth
  write (datehr(n)( 7: 8),'(i2.2)') lday
  write (datehr(n)( 9:10),'(i2.2)') lhour
  write (datehr(n)(11:12),'(i2.2)') lmin

  write (6,'(a,i3,1x,a12)') 'n datehr: ', n, datehr(n)
 end do
 write (6,*) 'ho la stringa con le date'
 write (6,*) ' '
 call flush(6)
 ! write (6,*) ' ok?'
 ! read (5,'(a1)') ans

!-----------------------------------------------------------------------
! read station archive:
!-----------------------------------------------------------------------
 write (6,*) 'now read station archive'
 call anagreadUTM(FULFISTAG,atid,ax,ay,alon,alat,ahgt, KKA)
 write (6,*) 'archive read'

! !deb
!  write (6,'(A4,1X,A8,1X,2A9,2A9,A6)') &
! &   'k','ID   ',' GBx  ',' GBy  ',' lon  ',' lat  ',' hgt '
!  do k=1,KKA
!   write (6,'(i4,1x,a8,1x,2f9.0,2f9.4,f6.0)') &
! &   k,atid(k),ax(k),ay(k),alon(k),alat(k),ahgt(k)
!   if (mod(k,40) == 0) then
!    write (6,*) ' ok?'
!    read (5,'(a1)') ans
!   end if
!  end do
 write (6,*) 'KKA=',KKA
 write (6,*) ' '
 call flush(6)
 ! write (6,*) ' ok?'
 ! read (5,'(a1)') ans

!-----------------------------------------------------------------------
! compute ZS matrix
!-----------------------------------------------------------------------
 write (6,*) 'building matrix ZS..sigd eps', sigd, eps
 ZS= 0.d0
 do k=1,KKA
  do k2=1,KKA
   ! horizontal covariance:
   rah= distan(alon(k),alat(k), alon(k2),alat(k2))
   ssh= 0.d0
   arg= rah /sigd
   if (arg < 7.d0) ssh= exp(-arg*arg/2.d0)
   ! 2d covariance:
   ZS(k,k2)= ssh
  end do
  ZS(k,k)= ZS(k,k) +eps
 end do
 write (6,*) ' ho finito ZS KKA KKAxKKAx8', KKA, KKA*KKA*8
 write (6,*) ' '
 call flush(6)
 ! write (6,*) ' ok?'
 ! read (5,'(a1)') ans

 !------------------------------------------------------------------------------
 if (griglia) then

  !-----------------------------------------------------------------------
  ! read grid and orography
  ! G4X G4Y : original grid coordinates (lon-lat; kmetric, Gauss-Boaga...)
  ! G4LON G4LAT : lon lat of gridpoints
  ! G4HGT: height of gridpoints
  ! G4MSK: admin boundary mask
  !-----------------------------------------------------------------------
  write (6,*) 'now read orography and grid'
  call GetUTMOrography(FULFIORO, IIG, g4x,g4y, g4lon,g4lat,g4hgt,g4msk,g4hui)
  write (6,*) 'total number of gridpoints IM IIG: ',IM,IIG
  ! write (6,*) ' ok?'
  ! read (5,'(a1)') ans

  !-----------------------------------------------------------------------
  ! assign to R8 variables
!  xacum24h= UNDEF
  i= 0
  do ig=1,IIG
   ! punti interni attivi
   if (g4msk(ig) > 0.d0) then
    i= i +1
    xlon(i)= dble(g4lon(ig))
    xlat(i)= dble(g4lat(ig))
    xhgt(i)= dble(g4hgt(ig))
!    xacum24h(i)= 0.d0
   end if
  end do
  II= i
  write (6,*) 'lon lat oro ok'
  write (6,*) 'number of active gridpoints II IEFF: ',II,IEFF
  write (6,*) ' '
  call flush(6)
  ! write (6,*) ' ok?'
  ! read (5,'(a1)') ans

  !-----------------------------------------------------------------------
  ! compute G matrix
  !-----------------------------------------------------------------------
  write (6,*) 'building G matrix...', sigd
  ! GG for active gridpoints only (inside administrative boundary)
  do i=1,II
   do k=1,KKA
    rah= distan(xlon(i),xlat(i), alon(k),alat(k))
    ggh= 0.d0
    arg= rah /sigd
    if (arg < 7.d0) ggh= exp(-arg*arg/2.d0)
    ! 2d covariance:
    GG(i,k)= ggh
   end do
  end do
  write (6,*) ' ho finito G II KKA IIxKKAx8=', II, KKA, II*KKA*8
  write (6,*) ' '
  call flush(6)
  ! write (6,*) ' ok?'
  ! read (5,'(a1)') ans

 end if ! griglia

!-----------------------------------------------------------------------
! inizializza KKP, SRNVP, JARKP, ovvero
! numero di stazioni, matrice inversa, e vettore di indici del passo precedente:
 KKDP= 0
 KKWP= 0
 SRNVDP= UNDEF
 SRNVWP= UNDEF
 jarkdp= 0
 jarkwp= 0

!-------------------------------------------------------------------------------
! this is the time loop:
 write (6,*) 'now start time loop'
 write (6,*) ' '
 call flush(6)

 open (10, file= OUTDIR(1:lodir) //'RMS.txt')
 write (10, '(A5,1X,A12,1X,5A12)') '   n', '  datehour  ', ' D-CVscore' &
&, ' W-rms(yo)', ' W-rms(ya)', ' Wrms(ya-yo)', ' W-CVscore'
 call flush(10)

 open (12, file= FULFICSV)
! anche il residuo sqrt((yo-ya)(yo-yaV)) , ma senza fare SCT ; CV-idi-wet e CV-idi-dry
!       123456789012345678901
!       12345;12345678901234;
 arow= '   Id; Year/Mo/Dy Hr;' &
&    //'      yo;       ya;      yaV;      res;  idi-wet;  idi-dry; idiV-wet; idiV-dry;'
!       12345678;123456789;123456789;123456789;123456789;123456789;123456789;123456789;
!       2345678901234567890123456789012345678901234567890123456789012345678901234567890
! totale 100 caratteri
 write (12,'(a100)') arow(1:100)
 call flush(12)

 ! read /write record number:
 norec= 0 ! plobsread
 ntrec= 0 ! sttread
 nwrec= 0

!-------------------------------------------------------------------------------
 do n=1,NN
!-------------------------------------------------------------------------------

  write (6,*) 'n=', n
  read (datehr(n)( 1: 4),'(i4)') lyear
  read (datehr(n)( 5: 6),'(i2)') lmonth
  read (datehr(n)( 7: 8),'(i2)') lday
  read (datehr(n)( 9:10),'(i2)') lhour
  read (datehr(n)(11:12),'(i2)') lmin
  write (6,*) 'date: ' //datehr(n)(1:4) //'/' //datehr(n)(5:6) &
& //'/' //datehr(n)(7:8) //' ' //datehr(n)(9:10) //'h' //datehr(n)(11:12)
  write (6,*) ' '
  call flush(6)
!-----------------------------------------------------------------------
! nomi file per output ascii
  FULFIASCIIa = OUTDIR(1:lodir) // 'prec_ana_32632_' // datehr(n)(1:10) // '.txt'
  FULFIASCIId = OUTDIR(1:lodir) // 'pdry_idi_32632_' // datehr(n)(1:10) // '.txt' 
  FULFIASCIIw = OUTDIR(1:lodir) // 'pwet_idi_32632_' // datehr(n)(1:10) // '.txt' 

  ! write (6,'(A)') ' ok?'
  ! read (5,'(a1)') ans
!-----------------------------------------------------------------------
! data-ora per file CSV:
  arow= ''
  arow(6:7)= '; '
  write (arow(8:21),'(a)') datehr(n)(1:4) //'/' //datehr(n)(5:6) &
& //'/' //datehr(n)(7:8) //' ' //datehr(n)(9:10) //';'

  !-----------------------------------------------------------------------
  ! read temperature analysis on station:
  call sttread(FULFITST, ntrec, KKT, tstid, tsta)
  write (6,*) 'ho letto TSTA'

  ! !deb
  ! do k=1,KKT
  !  write (6,'(A,A8,A2,F8.3)') ' id Ta: [',tstid(k),'] ',tsta(k)
  !  if (mod(k,40) == 0) then
  !   write (6,*) ' ok?'
  !   read (5,'(a1)') ans
  !  end if
  ! end do

  write (6,*) 'KKT KKA : ', KKT, KKA
  write (6,*) ' '
  call flush(6)
  ! write (6,*) ' ok?'
  ! read (5,'(a1)') ans

  !-----------------------------------------------------------------------
  ! INPUT observations (STID yo):
  write (6,*) 'now read observations...'
  KK= 0
  stid= ''
  so= UNDEF

  call plobsread(FULFIOB, norec, nvars, lpospluv, lpostemp, FULFIPRSK, KKT, tstid, tsta &
&                , KK, stid, so)
  write (6,*) 'back from plobsread'

  ! !deb
  !  do k=1,KK
  !   write (6,'(a,i5,a1,a8,a1,f12.2)') ' k id yo', k,' ',stid(k),' ',so(k)
  !   if (mod(k,20) == 0) then
  !    write (6,*) ' ok?'
  !    read (5,'(a1)') ans
  !   end if
  !  end do

  write (6,*) 'KK= ', KK
  write (6,*) ' '
  call flush(6)
  ! write (6,*) ' ok?'
  ! read (5,'(a1)') ans

  !-----------------------------------------------------------------------
  ! index of station k in archive
  ! data not present in archive are eliminated
  write (6,*) ' assegno jarkwet jarkdry e ao'
  KKWET= 0
  KKDRY= 0
  jarkwet= 0
  jarkdry= 0
  ao= UNDEF
  k= 1
  do while (k <= KK)
   j= 1
   do while (j <= KKA .and. stid(k) /= atid(j))
    j= j +1
   end do

   if (j > KKA) then

! blacklist:
!    if (j > KKA &
!&   .or. stid(k) == 'Lo00041 ' &
!     &   .or. (stid(k) == 'Lo00503 ' .and. n >= 18)
!&      ) then
!diagnostica:
!    if (stid(k) == 'Lo00649 ') write (6,*) 'blacklist  649***',atid(j)

    write (6,*) 'stazione fuori archivio: ',k,' ', stid(k)
    KK= KK -1
    do m=k,KK
     stid(m)= stid(m+1)
     so(m)= so(m+1)
    end do

   else

    ao(j)= so(k) ! observed value

    ! WET/DRY:
    if (ao(j) > 0.d0) then
     KKWET= KKWET +1
     jarkwet(KKWET)= j   ! index array
     ao(j)= reflrain(ao(j))
     ! write (6,'(i4,1x,a8,a1,a8,1x,i4,f9.2)') KKWET,stid(k),'=',atid(j),jarkwet(KKWET), ao(j)
    else
     KKDRY= KKDRY +1
     jarkdry(KKDRY)= j   ! index array
     ! write (6,'(i4,1x,a8,a1,a8,1x,i4,f9.2)') KKDRY,stid(k),'=',atid(j),jarkdry(KKDRY), ao(j)
    end if

    k= k +1
   end if

   ! if (mod(k,20) == 0) then
   !  write (6,*) ' ok?'
   !  read (5,'(a1)') ans
   ! end if
  end do

  write (6,*) 'assegnati jarkwet jarkdry, KKWET KKDRY: ', KKWET, KKDRY
  write (6,*) ' KK= KKDRY+KKWET', KK, (KKDRY+KKWET)
  write (6,*) ' '
  call flush(6)
  ! write (6,*) ' ok?'
  ! read (5,'(a1)') ans

!-------------------------------------------------------------------------------
  write (6,*) ' valori DRY:'
  write (6,'(2a4,1x,a8,1x,4a10,a7,a8)') 'k','j','id ','x','y','lon','lat','hgt','yo'
  do k=1,KKDRY
   j= jarkdry(k)
   write (6,'(2i4,1x,a8,1x,2f10.0,2f10.5,f7.0,f8.1)') &
&   k, j, atid(j), ax(j), ay(j), alon(j), alat(j), ahgt(j), ao(j)
  end do
  write (6,*) 'KKDRY=', KKDRY
  write (6,*) ' '
  call flush(6)
  ! write (6,*) ' ok?'
  ! read (5,'(a1)') ans

  write (6,*) ' valori WET:'
  write (6,'(2a4,1x,a8,1x,4a10,a7,a8)') 'k','j','id ','x','y','lon','lat','hgt','yo'
  do k=1,KKWET
   j= jarkwet(k)
   write (6,'(2i4,1x,a8,1x,2f10.0,2f10.5,f7.0,f8.1)') &
&   k, j, atid(j), ax(j), ay(j), alon(j), alat(j), ahgt(j), ao(j)
  end do
  write (6,*) 'KKWET=', KKWET
  write (6,*) ' '
  call flush(6)
  ! write (6,*) ' ok?'
  ! read (5,'(a1)') ans

  !-----------------------------------------------------------------------------
  if ((KKDRY +KKWET) < KKMIN) then
   KKDP= 0
   jarkdp= 0
   SRNVDP= UNDEF

   KKWP= 0
   jarkwp= 0
   SRNVWP= UNDEF

   do j=1,KKA ! saran pochi, ma ripristino i valori originali:
    if (ao(j) > 0) ao(j)= rainrefl(ao(j))
   end do

   aidid= 0.d0
   aididv= 0.d0
   aidiw= 0.d0
   aidiwv= 0.d0
   aa= UNDEF
   aav= UNDEF
   asctres= UNDEF

   ! CV-score file:
   write (10, '(i5,1x,a12,1x,5i12)') n, datehr(n) &
&                                  , int(UNDEF), int(UNDEF), int(UNDEF), int(UNDEF), int(UNDEF)
   call flush(10)

   write (6,*) ' !!! KK < KKMIN !!! KKDRY KKWET KKMIN : ',KKDRY, KKWET, KKMIN
   call flush(6)
   ! write (6,*) ' ok?'
   ! read (5,'(a1)') ans

  !-----------------------------------------------------------------------------
  else
   ! (KKDRY+KKWET) >= KKMIN:

   !-----------------------------------------------------------------------
   ! DRY
   if (KKDRY == 0) then
    KKDP= 0
    jarkdp= 0
    SRNVDP= UNDEF

    aidid= 0.d0
    aididv= 0.d0

   else 
    ! KKDRY > 0; (KKDRY+KKWET) >= KKMIN:

    ! update matrix (S+R)^-1
    write (6,*) 'DRY: updating inverse matrix (S+R)^-1 ...', KKDP, '-->', KKDRY
    ! reset SRNVD (not the previous SRNVDP):
    SRNVD= UNDEF
    call udsrnv(KTOT, KKDP, jarkdp, SRNVDP, KKDRY, jarkdry, ZS, SRNVD)
    write (6,*) 'DRY: back from updating inverse matrix'
    write (6,*) ' '
    call flush(6)
    ! write (6,*) ' ok?'
    ! read (5,'(a1)') ans

    ! elementi diagonali di W=HK
    write (6,*) 'DRY: calcolo diag(W), eps=', eps
    dialdmx= 0.d0
    diald= UNDEF
    didid= UNDEF
    do k=1,KKDRY
     ! Wkk:
     diald(k)= 1.d0 -eps *SRNVD(k,k)
     if (diald(k) > dialdmx) dialdmx= diald(k)
     !didi:
     didid(k)= 1.d0
    end do
    ! write (6,'(5f16.8)') (diald(m),m=1,KKDRY)
    write (6,*) 'DRY: dialdry max=', dialdmx
    ! write (6,*) ' ok?'
    ! read (5,'(a1)') ans
    write (6,*) ' '
    call flush(6)

    !zidid
    zidid= UNDEF
    write (6,*) 'DRY: now compute zidi=(S+R)^-1*[1 1...1] '
    do k=1,KKDRY
     zidid(k)= 0.d0
     do k2=1,KKDRY
      zidid(k)= zidid(k) +SRNVD(k,k2) *didid(k2)
     end do
    end do
    write (6,*) 'DRY: zidid ok '
    write (6,*) ' '
    call flush(6)
    ! write (6,*) ' ok?'
    ! read (5,'(a1)') ans

    ! compute IDI-DRY - stimo comunque su TUTTE le stazioni dell'archivio:
    write (6,*)'IDI-DRY on obs sites...'
    aidid= 0.d0
    do j=1,KKA
     zum= 0.d0
     do k2=1,KKDRY
      j2= jarkdry(k2)
      zum= zum+ ZS(j,j2)* zidid(k2)
      ! sulla diagonale di ZS c'e' 1+eps, devo togliere eps*yz
      if (j2 == j) then
       zum= zum -eps *zidid(k2) 
      end if
     end do
     aidid(j)= zum
    end do
    write (6,*) 'IDI-DRY ok'
    write (6,*) ' '
    call flush(6)
    ! write (6,*) ' ok?'
    ! read (5,'(a1)') ans

    ! DRY: Cross Validation
    write (6,*) 'DRY: computing CV-IDI...'
    aididv= aidid
    adiald= UNDEF
    do k=1,KKDRY
     j= jarkdry(k)
     ! tre versioni della stessa cosa:
     ! 1:  aididv(j)= aidid(j) -(diald(k)/(1.d0 -diald(k))) *(1.d0 -aidid(j))
     ! 2:  aididv(j)= (aidid(j) -diald(k)) /(1.d0 -diald(k))
     ! 3:
     aididv(j)=  1.d0 +(aidid(j) -1.d0) /(1.d0 -diald(k))
     adiald(j)= diald(k)
    end do
    write (6,*) 'DRY: CV-IDI computed'
    write (6,*) ' '
    call flush(6)
    ! write (6,*) ' ok?'
    ! read (5,'(a1)') ans

   end if ! KKDRY == 0 

   !-----------------------------------------------------------------------
   ! WET
   if (KKWET == 0) then
    KKWP= 0
    jarkwp= 0
    SRNVWP= UNDEF

    aa= 0.d0
    aav= 0.d0
    aidiw= 0.d0
    aidiwv= 0.d0

   else
    ! KKWET > 0; (KKDRY+KKWET) >= KKMIN:

    ! update matrix (S+R)^-1
     write (6,*) 'WET: updating inverse matrix (S+R)^-1 ...', KKWP, '-->', KKWET
    ! reset SRNVW (not the previous SRNVWP):
    SRNVW= UNDEF
    call udsrnv(KTOT, KKWP, jarkwp, SRNVWP, KKWET, jarkwet, ZS, SRNVW)
    write (6,*) 'WET: back from updating inverse matrix'
    write (6,*) ' '
    call flush(6)
    ! write (6,*) ' ok?'
    ! read (5,'(a1)') ans

    ! elementi diagonali di W=HK
    write (6,*) 'WET: calcolo diag(W), eps=', eps
    dialwmx= 0.d0
    dialw= UNDEF
    do k=1,KKWET
     ! Wkk:
     dialw(k)= 1.d0 -eps *SRNVW(k,k)
     if (dialw(k) > dialwmx) dialwmx= dialw(k)
    end do
    ! write (6,'(5f16.8)') (dialw(m),m=1,KKWET)
    write (6,*) 'WET: dialwet max=', dialwmx
    ! write (6,*) ' ok?'
    ! read (5,'(a1)') ans
    write (6,*) ' '
    call flush(6)

    ! innovation
    write (6,*) 'WET: compute innovation'
    ydmax= 0.d0
    res= 0.d0
    yd= UNDEF
    didiw= UNDEF
    do k=1,KKWET
     yd(k)= ao(jarkwet(k))
     if (abs(yd(k)) > ydmax) ydmax= abs(yd(k))
     res= res+ yd(k)*yd(k)
     didiw(k)= 1.d0
    end do
    if (KKWET > 0) then
     res= sqrt(res/KKWET)
    else
     res= UNDEF
    end if
    write (6,*) 'WET: rms & max (yo-yb)=', res, ydmax
    write (6,*) ' '
    call flush(6)
    ! write (6,*) ' ok?'
    ! read (5,'(a1)') ans
  
    ! compute ZIDI
    zidiw= UNDEF
    write (6,*) 'WET: now compute zidi=(S+R)^-1*[1 1...1] '
    do k=1,KKWET
     zidiw(k)= 0.d0
     do k2=1,KKWET
      zidiw(k)= zidiw(k) +SRNVW(k,k2) *didiw(k2)
     end do
    end do
    write (6,*) 'WET: zidiw ok '

    ! compute YZ
    write (6,*) 'WET: now compute yz=(S+R)^-1*d '
    az= UNDEF
    do k=1,KKWET
     j= jarkwet(k)
     az(j)= 0.d0
     do k2=1,KKWET
      az(j)= az(j) +SRNVW(k,k2) *yd(k2)
     end do
    end do
    write (6,*) 'WET: yz ok ',KKWET
    write (6,*) ' '
    call flush(6)
    ! write (6,*) ' ok?'
    ! read (5,'(a1)') ans

    ! compute analysis on obs sites
    ! stimo comunque su TUTTE le stazioni dell'archivio:
    write (6,*)'analysis on obs sites...'
    ady= UNDEF
    aidiw= 0.d0
    do j=1,KKA
     ady(j)= 0.d0
     zum= 0.d0
     do k2=1,KKWET
      j2= jarkwet(k2)
      ady(j)= ady(j)+ ZS(j,j2)* az(j2)
      zum= zum+ ZS(j,j2)* zidiw(k2)
      ! sulla diagonale di ZS c'e' 1+eps, devo togliere eps*yz
      if (j2 == j) then
       ady(j)=  ady(j) -eps *az(j2)
       zum= zum -eps *zidiw(k2) 
      end if
     end do
     aidiw(j)= zum

     if (ao(j) == UNDEF) then
      ! se non ho ao, decido con gli IDI:
      if (0.6*aidid(j) >= aidiw(j)) then
       aa(j)= 0.d0
      else
       aa(j)= ady(j)
      end if
     else
      if (ao(j) > 0) aa(j)= ady(j)
      if (ao(j) == 0) aa(j)= 0.d0
     end if

    end do ! j=1,KKA
    
    ! residui:
    qres= 0.d0
    pres= 0.d0
    do k=1,KKWET
     j= jarkwet(k)
     qres= qres+ (aa(j))**2
     pres= pres+ (aa(j)-ao(j))**2
    end do
    qres= sqrt(qres/KKWET)
    pres= sqrt(pres/KKWET)
    write (6,*) 'WET: rms(yo   ) ',  res
    write (6,*) 'WET: rms(ya   ) ', qres
    write (6,*) 'WET: rms(ya-yo) ', pres
    write (6,*) ' '
    call flush(6)
    ! write (6,*) ' ok?'
    ! read (5,'(a1)') ans

    ! Cross Validation
    write (6,*) 'WET: computing CV analysis...'
    aav= aa
    aidiwv= aidiw
    adialw= UNDEF
    do k=1,KKWET
     j= jarkwet(k)
     ! tre versioni della stessa cosa:
     ! 1:
     ! aav(j) =  aa(j) -(dialw(k)/(1.d0 -dialw(k))) *(ao(j) -aa(j))
     ! aidiwv(j)= aidiw(j) -(dialw(k)/(1.d0 -dialw(k))) *(1.d0 -aidiw(j))
     ! 2:
     ! aav(j) =  (aa(j) -dialw(k) *ao(j)) /(1.d0 -dialw(k))
     ! aidiwv(j)= (aidiw(j) -dialw(k)) /(1.d0 -dialw(k))
     ! 3:
     aav(j) = ao(j) +(aa(j) -ao(j)) /(1.d0 -dialw(k))
     aidiwv(j)=  1.d0 +(aidiw(j) -1.d0) /(1.d0 -dialw(k))

     ! CV, infine decido con gli IDI:
     if (0.6*aididv(j) >= aidiwv(j)) then
      aav(j)= 0.d0
     end if

     adialw(j)= dialw(k)
    end do

    ! stazioni "DRY" aa=ao=0, ma come CV diventano come UNDEF e aav dipende dagli IDI:
    do k=1,KKDRY
     j= jarkdry(k)
     if (0.6*aididv(j) < aidiwv(j)) then
      aav(j)= ady(j) ! ce l'ho ancora in memoria
      write (6,*) 'DRY-CV -> WET (idiV) - CVIDI-d CVIDI-w aav', aididv(j), aidiwv(j), aav(j)
     end if
    end do
  
    write (6,*) 'WET: CV analysis computed'
    write (6,*) ' '
    call flush(6)
    ! write (6,*) ' ok?'
    ! read (5,'(a1)') ans

   end if ! KKWET == 0

   !-----------------------------------------------------------------------
   ! residuo SCT
   write (6,*) ' calcolo residuo SCT'
   asctres= UNDEF

   do k=1,KKDRY
    asctres(jarkdry(k))= 0.d0
   end do

   do k=1,KKWET
    j=jarkwet(k)
    absd= abs( (ao(j) -aa(j)) *(ao(j) -aav(j)) )
    asctres(j)= sqrt(absd)
   end do

   write (6,*) ' '
   call flush(6)
   ! write (6,*) ' ok?'
   ! read (5,'(a1)') ans

   !-----------------------------------------------------------------------
   ! back to rain:
   do j=1,KKA
    if (ao(j) > 0) ao(j)= rainrefl(ao(j))
    if (aa(j) > 0) aa(j)= rainrefl(aa(j))
    if (aav(j) > 0) aav(j)= rainrefl(aav(j))
    if (asctres(j) > 0) asctres(j)= rainrefl(asctres(j))
   end do

   !-----------------------------------------------------------------------
   ! CV-score RAIN
   !-----------------------------------------------------------------------
   cvscwet= UNDEF
   resw= UNDEF
   qresw= UNDEF
   presw= UNDEF
   if (KKWET > 0) then
    cvscwet= 0.d0
    resw= 0.d0
    qresw= 0.d0
    presw= 0.d0
    do k=1,KKWET
     j= jarkwet(k)
     cvscwet= cvscwet +(ao(j) -aav(j))**2
     resw= resw +ao(j)**2
     qresw= qresw +aa(j)**2
     presw= presw +(ao(j) -aa(j))**2
    end do
    cvscwet= sqrt(cvscwet /KKWET)
    resw= sqrt(resw /KKWET)
    qresw= sqrt(qresw /KKWET)
    presw= sqrt(presw /KKWET)
   end if

   cvscdry= UNDEF
   if (KKDRY > 0) then
    cvscdry= 0.d0
    do k=1,KKDRY
     j= jarkdry(k)
     cvscdry= cvscdry +(ao(j) -aav(j))**2
    end do
    cvscdry= sqrt(cvscdry /KKDRY)
   endif

! write (10, '(A5,1X,A12,1X,8A12)') '   n', '  datehour  ' &
!&, ' D rms(yo)', ' D rms(ya)', ' Drms(ya-yo)', ' D CV-score' &
!&, ' W rms(yo)', ' W rms(ya)', ' Wrms(ya-yo)', ' W CV-score'

   ! write CVscwet CVskdry:
   if (KKDRY > 0 .and. KKWET >0) then
    write (10, '(i5,1x,a12,1x,5f12.4)') n, datehr(n), cvscdry, resw, qresw, presw, cvscwet
   else if (KKWET > 0) then
    write (10, '(i5,1x,a12,1x,i12,4f12.4)') n, datehr(n), int(UNDEF), resw, qresw, presw, cvscwet
   else if (KKDRY > 0) then
    write (10, '(i5,1x,a12,1x,f12.4,4i12)') n, datehr(n), cvscdry &
&                                         , int(UNDEF), int(UNDEF), int(UNDEF), int(UNDEF)
   end if
   call flush(10)

   write (6, '(a,i5,1x,a12,1x,5f12.4)') 'CVscore: ', n, datehr(n), cvscdry &
&                                                               , resw, qresw, presw, cvscwet
   write (6,*) ' '
   call flush(6)
   ! write (6,*) ' ok?'
   ! read (5,'(a1)') ans

  end if ! (KKDRY+KKWET) < KKMIN
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! output CSV station file
  !          123456789012345678901
  !          12345;12345678901234;
  !   arow= '   Id; Year/Mo/Dy Hr;' &
  ! &     //'      yo;       ya;      yaV;      res;  idi-wet;  idi-dry; idiV-wet; idiV-dry;'
  !          12345678;123456789;123456789;123456789;123456789;123456789;123456789;123456789;
  !          2345678901234567890123456789012345678901234567890123456789012345678901234567890
  !           1         2
  !              3         4         5         6         7         8         9        10
  do j=1,KKA
   arow(1:6)= '     ;'
   arow(1:5)= atid(j)(3:7)

   arow(22:30)= '        ;'
   write (arow(22:29),'(f8.1)') ao(j)
   if (ao(j) == UNDEF)  write (arow(22:29),'(i8)') int(UNDEF)

   arow(31:40)= '         ;'
   write (arow(31:39),'(f9.2)') aa(j)
   if (aa(j) == UNDEF)  write (arow(31:39),'(i9)') int(UNDEF)

   arow(41:50)= '         ;'
   write (arow(41:49),'(f9.2)') aav(j)
   if (aav(j) == UNDEF)  write (arow(41:49),'(i9)') int(UNDEF)

   arow(51:60)= '         ;'
   write (arow(51:59),'(f9.2)') asctres(j)
   if (asctres(j) == UNDEF)  write (arow(51:59),'(i9)') int(UNDEF)

   arow(61:70)= '         ;'
   write (arow(61:69),'(f9.4)') aidiw(j)
   if (aidiw(j) == UNDEF)  write (arow(61:69),'(i9)') int(UNDEF)

   arow(71:80)= '         ;'
   write (arow(71:79),'(f9.4)') aidid(j)
   if (aidid(j) == UNDEF)  write (arow(71:79),'(i9)') int(UNDEF)

   arow(81:90)= '         ;'
   write (arow(81:89),'(f9.4)') aidiwv(j)
   if (aidiwv(j) == UNDEF)  write (arow(81:89),'(i9)') int(UNDEF)

   arow(91:100)= '         ;'
   write (arow(91:99),'(f9.4)') aididv(j)
   if (aididv(j) == UNDEF)  write (arow(91:99),'(i9)') int(UNDEF)

   write (12,'(a100)') arow(1:100)
   call flush(12)
  end do ! j
  write (6,*) 'analysis etc. written on CSV station file'
  write (6,*) ' '
  call flush(6)
  ! write (6,*) ' ok?'
  ! read (5,'(a1)') ans

  !-----------------------------------------------------------------------------
  ! GRID
  if (griglia) then

   !-----------------------------------------------------------------------------
   if ((KKDRY +KKWET) < KKMIN) then
    dx= UNDEF
    xidid= 0.d0
    xidiw= 0.d0
    xa= UNDEF

   !-----------------------------------------------------------------------------
   else
    ! (KKDRY+KKWET) >= KKMIN:

    !-----------------------------------------------------------------------
    ! analysis estimate on grid
    dx= UNDEF
    xidid= 0.d0
    xidiw= 0.d0
    xa= 0.d0

    write (6,*) ' computing analysis on grid x=G*z...', datehr(n)
    write (6,*) ' calcolo G*z...', KKA, KKDRY, KKWET, II
    res= 0.d0
    IW= 0
    do i=1,II ! GG for active gridpoints only

     ! DRY AREA IDI
     xum= 0.d0
     do k=1,KKDRY
      j= jarkdry(k)
      xum= xum+ GG(i,j)*zidid(k)
     end do
     xidid(i)= xum

     ! WET AREA
     dx(i)= 0.d0
     xum= 0.d0
     do k=1,KKWET
      j= jarkwet(k)
      dx(i)= dx(i) +GG(i,j) *az(j)
      xum= xum+ GG(i,j)*zidiw(k)
     end do
     xidiw(i)= xum

     if (0.6*xidid(i) >= xidiw(i)) then
      xa(i)= 0.d0
      dx(i)= UNDEF
     else
      xa(i)= 0.d0
      if (dx(i) > 0.d0) xa(i)= rainrefl(dx(i))
      res= res+ xa(i)**2
      IW= IW +1
     end if

     ! cum 24h:
!     xacum24h(i)= xacum24h(i) +xa(i)

    end do

    res= sqrt(res/IW)
    write (6,*) ' analysis on grid ok, rms(xaw) ', res
    write (6,'(a,5f14.6)') 'xa ', xa(1), xa(II/3), xa(II/2), xa(2*II/3), xa(II)
!    write (6,'(a,5f14.6)') 'xacum24h ' &
!&   , xacum24h(1), xacum24h(II/3), xacum24h(II/2), xacum24h(2*II/3), xacum24h(II)
    write (6,*) ' '
    call flush(6)
    ! write (6,*) ' ok?'
    ! read (5,'(a1)') ans

   end if ! (KKDRY+KKWET) < KKMIN
   !-----------------------------------------------------------------------

  end if ! griglia
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! output obs analysis estimate, le salvo tutte con gli UNDEF:
  stid= ''
  s4x= UNDEF
  s4y= UNDEF
  shgt= UNDEF
  so= UNDEF
  sa= UNDEF
  sav= UNDEF
  do j=1,KKA
   stid(j)= atid(j)
   s4x(j)= ax(j)
   s4y(j)= ay(j)
   shgt(j)= ahgt(j)
   so(j)= ao(j)
   sa(j)= aa(j)
   sav(j)= aav(j)
  end do

  write (6,*) 'writing obs analysis...', KKA, FULFIST
  call plstzwri(FULFIST, nwrec, KKA, stid, s4x, s4y, shgt, so, sa, sav)

  write (6,*) ' analysis etc. written on GrADS station file '
  write (6,*) ' '
  call flush(6)
  ! write (6,*) ' ok?'
  ! read (5,'(a1)') ans

  !-----------------------------------------------------------------------------
  if (griglia) then
   ! output analysis on grid
   i= 0
   do ig=1,IIG
    gxa(ig)= sngl(UNDEF)
    gidiw(ig)= sngl(UNDEF)
    gidid(ig)= sngl(UNDEF)
    if (g4msk(ig) > 0.d0) then
     i= i +1
     gxa(ig)= sngl(xa(i))
     gidiw(ig)= sngl(xidiw(i))
     gidid(ig)= sngl(xidid(i))
    end if
   end do
   write (6,*) 'writing grid analysis...', i, FULFIGR
   call plgriwri(n, FULFIGR, gxa, gidiw, gidid)
   write (6,*) ' grid analysis etc. written on binary file '
   call asciiUTM(FULFIASCIIa,gxa)
   call asciiUTM(FULFIASCIId,gidid)
   call asciiUTM(FULFIASCIIw,gidiw)
   write (6,*) ' grid analysis etc. written on ascii file '
   write (6,*) ' '
   call flush(6)
   ! write (6,*) ' ok?'
   ! read (5,'(a1)') ans

  end if ! griglia

  !-----------------------------------------------------------------------
  write (6,*) ' '
  write (6,*) 'end datehour: ' //datehr(n)(1:4) //' ' //datehr(n)(5:6) &
&  //' ' //datehr(n)(7:8) //' ' //datehr(n)(9:10) //'h' //datehr(n)(11:12)
  write (6,*) '----------------------------------------------------------------'
  write (6,*) ' '
  call flush(6)
  call flush(10)
  call flush(12)
  ! write (6,*) ' ok?'
  ! read (5,'(a1)') ans

!-----------------------------------------------------------------------
 end do ! end do n
!-----------------------------------------------------------------------

 !-----------------------------------------------------------------------
 ! close CVscore file:
 close (10)
 ! close file CSV:
 close (12)

!-----------------------------------------------------------------------
! output cum 24h
! if (griglia) then
!
!  ! output analysis on grid
!  i= 0
!  do ig=1,IIG
!   gxacum24h(ig)= sngl(UNDEF)
!   if (g4msk(ig) > 0.d0) then
!    i= i +1
!    gxacum24h(ig)= sngl(xacum24h(i))
!   end if
!  end do
!  write (6,*) 'writing grid analysis of 24-h accumulated precipitation...', i, FULFICUM
!  open (9, file= FULFICUM, form= 'UNFORMATTED', access= 'DIRECT', recl=IM*4)
!  write (9,rec=1) gxacum24h
!  close (9)
!  write (6,*) 'cum 24h written'
! end if ! griglia

 !-------------------------------------------------------------------------------
 stop
end program plzln19
!*******************************************************************************
! SUBROUTINES AND FUNCTIONS (here)
! sttread
! plobsread
! plstzwri
! plgriwri
! asciiUTM
! FUNCTIONS:
! rainrefl
! reflrain
!***************************************************************************************************
subroutine sttread(FULFIST, nrec, KK, stid, sa)
!-------------------------------------------------------------------------------
! read temperature analysis and CV analysis from GrADS "station" file
! from .ctl:
! VARS    6
! elev    0  99  station elvation ASL
! To      0  99  T observation
! Tb      0  99  T background
! Ta      0  99  T analysis
! TaV     0  99  CV - T analysis
! uhi     0  99  uhi
!
! header 28 bytes= NHEA*4 bytes   NHEA=7
! 1 ( 1):  4 stid(1:4) character(4) stid character(8)
! 2 ( 2):  8 stid(5:8) character(4)
! 2 ( 3): 12 lat/y    real(4)
! 3 ( 4): 16 lon/x    real(4)
! 4 ( 5): 20 tim      real(4)
! 5 ( 6): 24 nlev     integer(4)
! 6 ( 7): 28 nflag    integer(4)
! 
! record terminator: header with nlev=0
!-------------------------------------------------------------------------------
 implicit none
 integer, parameter :: KTOT=1000
 real(4), parameter :: UNDEF=-9999.d0

!-----------------------------------------------------------------------
 character(200) :: FULFIST
 ! integer, save :: nrec
 integer :: nrec
 integer :: KK
 character(8), dimension (KTOT) :: stid
 real(4), dimension (KTOT) :: sa ! , sav
 !FRANCESCO PROVE 20091102 SENZA UHI: integer, parameter :: nvars= 5
 integer, parameter :: nvars= 6

!-----------------------------------------------------------------------
 real(4), allocatable :: puff(:)

 real(4) :: tim, s4lox, s4lay
 integer :: nlev, nflag
 character(8) :: aid

!-----------------------------------------------------------------------
 integer j, k
 integer iid

!-----------------------------------------------------------------------
 integer :: lstr
! functions:
 integer :: lstrim

 !character(1) :: ans

!===============================================================================
 stid= 'Lo00000 '
 aid= ''

 !write (6,*) 'sttread: file :',trim(FULFIST)
 !write (6,*) 'sttread: nrec :',nrec
 !write (6,*) ' ok?'
 !read (5,'(a1)') ans

!-----------------------------------------------------------------------
 allocate(puff(nvars))

 open (8, file= FULFIST, form= 'UNFORMATTED', access= 'DIRECT', recl=4)

!-----------------------------------------------------------------------
! nrec is the record number of the last record previously read

! first obs:
 k= 1

!-----------------------------------------------------------------------
! read first header:
 nrec= nrec +1
 read (8,rec=nrec) aid(1:4)
 nrec= nrec +1
 read (8,rec=nrec) aid(5:8)
 nrec= nrec +1
 read (8,rec=nrec) s4lay
 nrec= nrec +1
 read (8,rec=nrec) s4lox
 nrec= nrec +1
 read (8,rec=nrec) tim
 nrec= nrec +1
 read (8,rec=nrec) nlev
 nrec= nrec +1
 read (8,rec=nrec) nflag

 !write (6,*) 'ID: ',aid
 !write (6,*) 's4lox: ',s4lox
 !write (6,*) 's4lay: ',s4lay
 !write (6,*) 'tim: ',tim
 !write (6,*) 'nlev: ',nlev
 !write (6,*) 'nflag: ',nflag
 !write (6,*) ' ok?'
 !read (5,'(a1)') ans

!-----------------------------------------------------------------------
! check for record terminator:
 do while (nlev /= 0)

 ! read variables:
  do j=1,nvars
   nrec= nrec +1
   read (8,rec=nrec) puff(j)
  end do

!-----------------------------------------------------------------------
! assign variables:
  ! write (6,'(a,a8,a,2f12.3)') 'ID: [' , aid ,'] obs=', puff(4), puff(7)

  ! STation ID:
  stid(k)= aid

  !lstr= lstrim(aid)
  ! write (6,'(a1,a8,a1,i3)') '[',aid,']', lstr
  !read (aid(3:lstr),*) iid                        ! 3: file GrADS analisi aid='Lo.....'
  !write (stid(k)(3:7),'(i5.5)') iid
  ! write (6,*) 'iid=',iid
  ! write (6,'(a1,a8,a2,i4)') '[',stid(k),']=', iid

  ! observation:
  ! 1 elev
  ! 2 To
  ! 3 Tb
  sa(k)= puff(4)
  ! sav(k)= puff(5)
  ! 6 uhi

  ! write (6,'(a,4f10.2)') 'sttread: puff', puff(1), puff(2), puff(3), puff(4)
  ! write (6,*) 'sttread: puff ', puff
  ! write (6,'(a,i5,1x,a8,2f9.2)') 'sttread: k id sa sav', k, stid(k), sa(k), sav(k)
  ! write (6,*) ' ok?'
  ! read (5,'(a1)') ans

  ! next obs:
  k= k +1

  !-----------------------------------------------------------------------
  ! read next header:
  nrec= nrec +1
  read (8,rec=nrec) aid(1:4)
  nrec= nrec +1
  read (8,rec=nrec) aid(5:8)
  nrec= nrec +1
  read (8,rec=nrec) s4lay
  nrec= nrec +1
  read (8,rec=nrec) s4lox
  nrec= nrec +1
  read (8,rec=nrec) tim
  nrec= nrec +1
  read (8,rec=nrec) nlev
  nrec= nrec +1
  read (8,rec=nrec) nflag

 end do ! go and check for record terminator

 ! record terminator found; number of observation at present time:
 KK= k -1
 close (8)
 deallocate(puff)
 write (6,*) 'sttread: read KK=',KK

!-------------------------------------------------------------------------------
 return
end subroutine sttread
!***************************************************************************************************
subroutine plobsread(FULFIOB, nrec, nvars, lpospluv, lpostemp, FULFIPRSK, KKT, tstid, tsta &
&                  , KK, stid, yo)
!-------------------------------------------------------------------------------
! read observations from GrADS "station" file
! eliminate missing data
! each observation is composed by the header (7*4byte) +puff (nvars*4byte)
! lpos position of desired variable after header = position in .ctl file
! observed variable = puff(lpos)
!
! header 28 bytes= NHEA*4 bytes   NHEA=7
! 1 ( 1):  4 stid(1:4) character(4) stid character(8)
! 2 ( 2):  8 stid(5:8) character(4)
! 2 ( 3): 12 lat/y    real(4)
! 3 ( 4): 16 lon/x    real(4)
! 4 ( 5): 20 tim      real(4)
! 5 ( 6): 24 nlev     integer(4)
! 6 ( 7): 28 nflag    integer(4)
! 
! record terminator: header with nlev=0
!-------------------------------------------------------------------------------
 implicit none
 integer, parameter :: KTOT=1000

! UNDEF for obs data is different from UNDEF in the code!
 real(4), parameter :: UNDEF=-999.9

!-----------------------------------------------------------------------
 character(200) :: FULFIOB, FULFIPRSK
 !integer, save :: nrec
 integer :: nrec
 integer :: KK, KKT
 character(8), dimension (KTOT) :: stid, tstid
 real(4), dimension (KTOT) :: yo, tsta, ytemp
 integer :: nvars, lpospluv, lpostemp

!-----------------------------------------------------------------------
 real(4), allocatable :: puff(:)

 real(4) :: tim, s4lox, s4lay
 integer :: nlev, nflag
 character(8) :: aid

!-----------------------------------------------------------------------
 real(4) :: yomin, yomax
 integer j, k
 integer iid

!-----------------------------------------------------------------------
! pluviometri non riscaldati:
 real(4), parameter :: tempnonrisk= 0.0d0
 character(8) :: stidrisk(KTOT)
 integer :: nriskid, ios, KSCA
 integer :: KRISK

!-----------------------------------------------------------------------
! analisi temperatura:
 real(8), parameter :: TAUNDEF=-9999.0

!-----------------------------------------------------------------------
 integer :: lstr
! functions:
 integer :: lstrim

 ! character(1) :: ans

!===============================================================================
! precipitation mm in 1 h
 yomin= 0.d0
 yomax= 200.d0

!-----------------------------------------------------------------------
 stid= 'Lo00000 '
 aid= ''

 ! write (6,*) 'plobsread: file :',trim(FULFIOB)
 ! write (6,*) 'plobsread: nrec :',nrec
 ! write (6,*) 'plobsread: nvars:',nvars
 ! write (6,*) 'plobsread: lpos pluv, temp:',lpospluv, lpostemp
 ! write (6,*) ' ok?'
 ! read (5,'(a1)') ans

!-----------------------------------------------------------------------
 allocate(puff(nvars))

 open (8, file= FULFIOB, form= 'UNFORMATTED', access= 'DIRECT', recl=4)

!-----------------------------------------------------------------------
! nrec is the record number of the last record previously read

! first obs:
 k= 1

!-----------------------------------------------------------------------
! read first header:
 nrec= nrec +1
 read (8,rec=nrec) aid(1:4)
 nrec= nrec +1
 read (8,rec=nrec) aid(5:8)
 nrec= nrec +1
 read (8,rec=nrec) s4lay
 nrec= nrec +1
 read (8,rec=nrec) s4lox
 nrec= nrec +1
 read (8,rec=nrec) tim
 nrec= nrec +1
 read (8,rec=nrec) nlev
 nrec= nrec +1
 read (8,rec=nrec) nflag

 !write (6,*) 'ID: ',aid
 !write (6,*) 's4lox: ',s4lox
 !write (6,*) 's4lay: ',s4lay
 !write (6,*) 'tim: ',tim
 !write (6,*) 'nlev: ',nlev
 !write (6,*) 'nflag: ',nflag
 !write (6,*) ' ok?'
 !read (5,'(a1)') ans

!-----------------------------------------------------------------------
! check for record terminator:
 do while (nlev /= 0)

 ! read variables:
  do j=1,nvars
   nrec= nrec +1
   read (8,rec=nrec) puff(j)
  end do

!-----------------------------------------------------------------------
! assign variables:

  ! STation ID:
  stid(k)=aid
  !stid(k)= 'Lo00000 '
  !lstr= lstrim(aid)
  !write (6,'(a1,a4,a1,a4)') ' ',bid(1),':', bid(2)
  !write (6,'(a1,a8,a1,i3)') '[',aid,']', lstr
  !read (aid(1:lstr),*) iid
  !write (stid(k)(3:7),'(i5.5)') iid
  !write (6,*) 'iid=',iid
  !write (6,'(a1,a8,a2,i4)') '[',stid(k),']=', iid

  ! observation:
  yo(k)= puff(lpospluv)
  ytemp(k)= puff(lpostemp)

  !write (6,'(a,4f10.2)') 'plobsread: puff', puff(1), puff(2), puff(3), puff(4)
  !write (6,*) 'plobsread: puff ', puff
!  if (yo(k) /= UNDEF) write (6,'(a,i5,1x,a8,2f8.1)') 'plobsread: k id pluv temp'&
!&                   , k, stid(k), yo(k), ytemp(k)
  !write (6,*) ' ok?'
  !read (5,'(a1)') ans

  ! next obs:
  k= k +1

  !-----------------------------------------------------------------------
  ! read next header:
  nrec= nrec +1
  read (8,rec=nrec) aid(1:4)
  nrec= nrec +1
  read (8,rec=nrec) aid(5:8)
  nrec= nrec +1
  read (8,rec=nrec) s4lay
  nrec= nrec +1
  read (8,rec=nrec) s4lox
  nrec= nrec +1
  read (8,rec=nrec) tim
  nrec= nrec +1
  read (8,rec=nrec) nlev
  nrec= nrec +1
  read (8,rec=nrec) nflag

 end do ! go and check for record terminator

 ! record terminator found; number of observation at present time:
 KK= k -1
 close (8)
 deallocate(puff)
 write (6,*) 'plobsread: read KK=',KK
 ! write (6,*) ' ok?'
 ! read (5,'(a1)') ans

!-----------------------------------------------------------------------
! check for non-heated raingauges file= 'pluvrisk.txt'
 open (3, file= FULFIPRSK)
 k= 0
 do
  read (3, *, iostat=ios) nriskid
  if (ios /= 0) exit
  k= k +1
  stidrisk(k)= 'Lo00000 '
  write (stidrisk(k)(3:7),'(i5.5)') nriskid
 end do
 close (3)
 KRISK= k

!-----------------------------------------------------------------------
! uso l'analisi di temperatura tsta
! faccio scansione dei pluviometri, se ytemp=UNDEF cerco l'analisi
! di temperatura tsta su quella stazione (ID) e metto ytemp= tsta

 do k=1,KK
  if (ytemp(k) < (UNDEF+0.5)) then ! cerco analisi solo se non ho temp osservata
   j= 1
   do while (j <= KKT .and. stid(k) /= tstid(j))
    j= j +1
   end do
   ! se j <= KKT (ovvero l'ho trovata) e tsta /= UNDEF, la assegno a ytemp(k):
   if (j <= KKT .and. tsta(j) > TAUNDEF) then ! TAUNDEF e' diverso da UNDEF
    write (6,*) 'trovata Ta per pluviometro senza To, id To Ta: ', stid(k), ytemp(k), tsta(j)
    ytemp(k)= tsta(j)
   end if
   ! se j > KKT non l'ho trovata e non faccio niente
  end if
 end do

!-----------------------------------------------------------------------
! controllo pluviometri non riscaldati (il controllo di plausibilita' e' fatto sotto)

 KSCA= 0
 ! write (6,*) ' '
 write (6,*) 'pluviometri riscaldati: ',KRISK
 write (6,*) 'controllo temperatura pluviometri non riscaldati :', (KK-KRISK)
 write (6,*) ' scarto se rain < 0.1  e  T <', tempnonrisk

 do k=1,KK

  ! cerco il pluviometro nella lista dei riscaldati:
  j= 1
  do while (j <= KRISK .and. stid(k) /= stidrisk(j))
   j= j +1
  end do

  if (j < KRISK) cycle ! pluviometro riscaldato: ok, vedi il prossimo

  write (6,*) 'pluviometro non riscaldato, id rain temp: ', stid(k), yo(k), ytemp(k)

  ! se arrivo qui, questo pluviometro non e' riscaldato, quindi:
  ! se (precipitazione<0.1) e fa troppo freddo, oppure non c'e' indicazione di temperatura
  !    (ovvero ne' termometro To ne' analisi Ta - UNDEF e' comunque minore di tempnonrisk),
  ! allora scarto questo dato

  if (yo(k) >= 0.0d0 .and. yo(k) < 0.1d0 .and. ytemp(k) < tempnonrisk) then

   write (6,*) 'plobsread: SCARTO PLUVIOMETRO NON RISCALDATO id rain temp: ' &
&  , stid(k), yo(k), ytemp(k)

   yo(k)= UNDEF  ! cosi' lo scarto
   KSCA= KSCA +1
  end if

 end do

 write (6,*) 'scartati ', KSCA, '   su', KK
 ! write (6,*) 'ok?'
 ! read (5,'(a1)') ans

!-----------------------------------------------------------------------
! check for missing data and climatological range
! check on station elevation is done in ANAGREAD

 write (6,*) 'plobsread: missing data and climatological range check min-max: ', yomin, yomax
 k= 0
 do while (k < KK)
  k= k +1

  if (yo(k) < yomin .or. yo(k) > yomax) then  ! anche UNDEF < yomin

   if (yo(k) /= UNDEF) write (6,*) 'plobsread: climat check failed id obs: ', stid(k), yo(k)

   KK= KK -1
   do j=k,KK
    stid(j)= stid(j+1)
    yo(j)= yo(j+1)
    ytemp(j)= ytemp(j+1) ! allineo anche ytemp anche se non serve piu'
   end do
   k= k -1

  end if

 end do
 write (6,*) 'plobsread: data climatologically checked and'
 write (6,*) '           missing data eliminated, KK=',KK
 ! write (6,*) ' ok?'
 ! read (5,'(a1)') ans

!-----------------------------------------------------------------------
! UNDEF elsewhere
 do j=KK+1,KTOT
  stid(j)= ' '
  yo(j)= UNDEF
  ytemp(j)= UNDEF
 end do

!-----------------------------------------------------------------------
 do k=1,KK
   write (6,'(a,i5,1x,a8,1x,2f9.1)') 'plobsread: k id pluv temp', k, stid(k), yo(k), ytemp(k)
 end do
 write (6,*) 'plobsread: KK=', KK
 ! write (6,*) ' ok?'
 ! read (5,'(a1)') ans

!-------------------------------------------------------------------------------
 return
end subroutine plobsread
!*******************************************************************************
subroutine plstzwri(FULFIST, nrec, KK, stid, s4x, s4y, shgt, so, sa, sav)
!-------------------------------------------------------------------------------
! write obs, analysis etc. on GrADS "station" file. "pl" : no background
!-------------------------------------------------------------------------------
 implicit none
 integer, parameter :: KTOT=1000
 real(4), parameter :: UNDEF=-9999.0

 character(200) :: FULFIST
 ! integer, save :: nrec
 integer :: nrec
 integer :: KK
 character(8), dimension (KTOT) :: stid
 real(4), dimension (KTOT) :: s4x, s4y, shgt, so, sa, sav

 real(4) :: tim
 integer :: nlev, nflag

 integer :: k

! character(1) :: ans

!===============================================================================
 open (8, file= FULFIST, form= 'UNFORMATTED', access= 'DIRECT', recl=4)

!-------------------------------------------------------------------------------
 tim= 0.
 nflag= 1
 nlev= 1

 !-----------------------------------------
 ! nrec is the record number of last previously written record
 do k=1,KK
  ! header:
  nrec= nrec +1
  write (8,rec=nrec) stid(k)(1:4)
  nrec= nrec +1
  write (8,rec=nrec) stid(k)(5:8)
  nrec= nrec +1
  write (8,rec=nrec) s4y(k)
  nrec= nrec +1
  write (8,rec=nrec) s4x(k)
  nrec= nrec +1
  write (8,rec=nrec) tim
  nrec= nrec +1
  write (8,rec=nrec) nlev
  nrec= nrec +1
  write (8,rec=nrec) nflag
  ! variables:
  nrec= nrec +1
  write (8,rec=nrec) shgt(k)
  ! observations:
  nrec= nrec +1
  write (8,rec=nrec) so(k)
  ! analysis:
  nrec= nrec +1
  write (8,rec=nrec) sa(k)
  ! CV analysis:
  nrec= nrec +1
  write (8,rec=nrec) sav(k)
 end do

 !-----------------------------------------
 ! write record terminator:
 nlev= 0
 nrec= nrec +1
 write (8,rec=nrec) stid(1)(1:4)
 nrec= nrec +1
 write (8,rec=nrec) stid(1)(5:8)
 nrec= nrec +1
 write (8,rec=nrec) s4y(1)
 nrec= nrec +1
 write (8,rec=nrec) s4x(1)
 nrec= nrec +1
 write (8,rec=nrec) tim
 nrec= nrec +1
 write (8,rec=nrec) nlev
 nrec= nrec +1
 write (8,rec=nrec) nflag

 !-----------------------------------------
 close (8)
 write (6,*) 'plstzwri: KK and last nrec=', KK, nrec

 !write (6,*) 'plstzwri: return '
 !write (6,*) 'ok?'
 !read (5,'(a1)') ans

!-------------------------------------------------------------------------------
 return
end subroutine plstzwri
!*******************************************************************************
subroutine plgriwri(n, FULFIGR, xa, xidiw, xidid)
!-------------------------------------------------------------------------------
 implicit none
 integer, parameter :: INX=251, INY=249, IM=INX*INY

 integer :: n
 character(200) :: FULFIGR
 real(4), dimension(IM) :: xa, xidiw, xidid

 integer :: nrec
! character(1) :: ans

!===============================================================================
 nrec= (n -1) *3

 WRITE (6,*) 'plgriwri: n IM nrec ', n, IM, nrec

 open (9, file= FULFIGR, form= 'UNFORMATTED', access= 'DIRECT', recl=IM*4)

! xa, xidiw, xidid
 nrec= nrec +1
 write (9,rec=nrec) xa
 nrec= nrec +1
 write (9,rec=nrec) xidiw
 nrec= nrec +1
 write (9,rec=nrec) xidid

 close (9)
 write (6,*) 'plgriwri: last nrec ', nrec
 !write (6,*) 'ok?'
 !read (5,'(a1)') ans

!-------------------------------------------------------------------------------
 return
end subroutine plgriwri
!*******************************************************************************
subroutine asciiUTM(FULFIASCII, x)
!-------------------------------------------------------------------------------
 implicit none
! here info from main code (could also be passed in call)
 integer, parameter :: INX=251, INY=249, IM=INX*INY
 real,parameter :: dXstart=  436000.0 , dYstart= 4918000.0 ,dDx= 1000.0

!input
 character(200) :: FULFIASCII
 real(4), dimension(IM) :: x

! della subroutine
 character(1) :: ans
 integer :: i,j,iErr
 CHARACTER(300000) :: sBuffer
 real(4) :: rUNDEF = -9999.0

!===============================================================================
 iErr=0
 open(15, file= FULFIASCII,status='REPLACE', IOSTAT=iErr)
 IF (iErr/=0) THEN
     PRINT *,'subroutine asciiUTM: error in opening file ',FULFIASCII
     STOP 1
 ENDIF
! scrivo l'header
 WRITE(15,"('ncols',1X,I10)",IOSTAT=iErr) iNX
 WRITE(15,"('nrows',1X,I10)",IOSTAT=iErr) iNY
 WRITE(15,"('xllcorner',1X,F20.5)",IOSTAT=iErr) dXstart
 WRITE(15,"('yllcorner',1X,F20.5)",IOSTAT=iErr) dYstart
 WRITE(15,"('cellsize',1X,F20.5)",IOSTAT=iErr) dDX
 WRITE(15,"('NODATA_value',1X,F20.5)",IOSTAT=iErr) rUNDEF
 DO i=1,iNY
    WRITE(sBuffer,*,IOSTAT=iErr) (x((INY-i)*INX+j),j=1,iNX)
    WRITE(15,"(A)",IOSTAT=iErr) TRIM(ADJUSTL(sBuffer))
 ENDDO
close (15)

end subroutine asciiUTM
!*******************************************************************************
real(8) function rainrefl(refl)
!-------------------------------------------------------------------------------
! refl= (log(prec/0.1))*3/0.46+10
! prec=0.1*exp(0.46*(refl-10)/3)
! refl [dbz] prec [mm/h]
!-------------------------------------------------------------------------------

 implicit none
 real(8) :: refl

 rainrefl= 0.1d0 *exp(0.1535d0 *(refl -10.d0))

 return
end function rainrefl
!*******************************************************************************
real(8) function reflrain(prec)
!-------------------------------------------------------------------------------
! refl=(LOG(prec/0.1))*3/0.46+10
! prec=0.1*exp(0.46*(refl-10)/3)
! refl [dbz] prec [mm/h]
!-------------------------------------------------------------------------------

 implicit none
 real(8) :: prec

 reflrain= (log(prec /0.1d0)) /0.1535d0 +10.d0

 return
end function reflrain
!*******************************************************************************
