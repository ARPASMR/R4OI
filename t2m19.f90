program t2m19
!-------------------------------------------------------------------------------
! analisi di temperatura - ARPA Lombardia
!
! $ t2m19 yyyymmddhh ninst <filedati> nvar npos <ana_oro_dir> <output_dir>
!
! orography grid Lombardia GB 1500m 177x174
!
! CODICE ORIGINALE:
! Francesco Uboldi, 2011
!-----------------------------------------------------------------------------
! MS, nov 2019
!  - mod anaguhiread e obsread per leggere dati e anagrafica scritti da R (UTM)
!  - mod parametri per lavorare su griglia UTM:
!    251x249  1.kmX1.km (come PRISMA) oromaskutm_mauri.dat 
!    (e grazie mille a Maurizio Favaron!)
!-------------------------------------------------------------------------------
 use SpatialStuff
 implicit none
!-------------------------------------------------------------------------------
 integer, parameter :: KTOT=1000, NTOT=240 & ! max 1 giorno e 1000 stazioni 
&                    , INX=251, INY=249, IM=INX*INY, IEFF=IM ! griglia UTM, non so IEFF
!&                    , INX=177, INY=174, IM=INX*INY, IEFF=11107 ! IEFF...li ho contati

 real(8), parameter :: UNDEF=-9999.d0

! at least KKMIN obs each time
 integer, parameter :: KKMIN= 20

!-------------------------------------------------------------------------------
! grid
 real(4), dimension(IM) :: g4x, g4y, g4hgt, g4msk, g4lon, g4lat, gxb, gxa, gidi

 real(8), dimension(IEFF) :: xb, xa, xidi, xlon, xlat, xhgt, dx

 real(8) :: xum, arg, distan, rah, rav, ggh, ggv
 integer :: i, II, ig, IIG

!-----------------------------------------------------------------------
! i-th line of matrix G:
 real(8), dimension(IEFF,KTOT) :: GG
 integer :: norec, nwrec
! logical itisthere

!-----------------------------------------------------------------------
! station archive
 character(8), dimension(KTOT) :: atid
 real(8), dimension(KTOT) :: ax, ay, alon, alat, ahgt &
&                          , ao, ab, aa, aav, az, aidi, aidiv &
&                          , anvdiasrnv, adial

 integer, dimension(KTOT) :: jark
 integer :: KKA

!-----------------------------------------------------------------------
! observations R4 variables (GrADS format I/O)
 real(4), dimension(KTOT) :: so, sb, sa, sav, s4x, s4y, shgt, slon, slat
 character(8), dimension (KTOT) :: stid
 integer :: KK, nvars, lpos

!-----------------------------------------------------------------------
! observations R8 variables
 real(8), dimension(KTOT) :: yd, zidi, didi
 real(8) :: ydmax

!-----------------------------------------------------------------------
! background field
 real(8) :: avlon, avlat, zinv, yinv, alpa, beta, gama, alpb, betb, gamb
 real(8), parameter :: dzbf= 100.d0
 real(8) :: bfabov, bfbelo, zabov, zbelo

!-----------------------------------------------------------------------
! DQC BF check, SCT
! bftol[hPa],  scttol[hPa^2]
 real(8), parameter :: bftol= 12.d0 &       ! BF tol [°C]
&                    , scttol= 17.d0*1.29  ! CV tol [°C**2]= T2*sigo2=(4.68°C)**2
 real(8) :: absd, absdmx
 real(8), dimension(KTOT) :: abfres, asctres
 character(5), dimension(KTOT) :: aflag
 integer nsca, kdmx
 character(2) :: abftol, ascttol
 character(200) :: arow

!-----------------------------------------------------------------------
! oi07 specific variables:
 real(8) :: ssh, ssv, dy
 integer :: k2, j2

! previous index array:
 integer, dimension(KTOT) :: jarkp
 integer :: KKP

! station-station b.f. error covariance matrix
 real(8), dimension(KTOT,KTOT) :: ZS, SRNV, SRNVP
 real(8), dimension(KTOT) :: dial, diasrnv
 real(8) :: dialmx, diasmn

!-----------------------------------------------------------------------
! analysis parameters and auxiliary variables
 real(8), parameter :: eps= 0.5d0, sigd= 20.d0, sigz= 500.d0

 real(8) :: zum, res, qres, pres

!-----------------------------------------------------------------------
 real(8) :: cvscore

!-----------------------------------------------------------------------
!UHI Urban Heat Island 
 real(8), dimension(KTOT) :: auhi
 real(8), dimension(IEFF) :: xuhi
 real(4), dimension(KTOT) :: suhi
 real(4), dimension(IM) :: g4uhi
 real(8), parameter :: attmin= 0.4d0
 real(8) :: att
 ! integer :: ix, jy
 ! area MI per CVscore
 character(200) :: FULFIAMI
 integer, dimension(KTOT) :: jarkami
 integer :: KAMI, koami, iaux
 character(8) :: amid
 real*8 :: cvskami

!-----------------------------------------------------------------------
! files:
 character(200) :: STAGFI, OROGRIFI, STANFI, GRIDFI, DATDIR, OUTDIR &
&   , FULFISTAG, FULFIORO, FULFIOB, FULFIST, FULFIGR, CSVFI, FULFICSV&
&   , FULFIASCIIa, FULFIASCIIb, FULFIASCIIi

!-----------------------------------------------------------------------
 character(10) :: adate
 character(12) :: datehour
 character(12), dimension(NTOT) :: datehr
 integer :: NN
 integer :: lhour, lmin, lyear, lmonth, lday, ldayfeb
 integer :: lodir, ladir

!-----------------------------------------------------------------------
 character(6) :: ann
 character(2) :: anvrs, apos

 logical :: griglia

!-----------------------------------------------------------------------
! functions:
 integer :: lstrim

!-----------------------------------------------------------------------
 integer :: j, k, m, n

!-----------------------------------------------------------------------
 character(1) :: ans

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
  write (6,*) 't2m19 yyyymmddhh ninst <filedati> nvar npos' &
& //' <ana_oro_dir> <output_dir>'
  write (6,*) ' '
  write (6,*) ' OPPURE, per stimare solo su stazioni, "-" prima di ninst:'
  write (6,*) 't2m19 yyyymmddhh -ninst <filedati> nvar npos' &
& //' <ana_oro_dir> <output_dir>'
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
 call getarg(5,apos)
 call getarg(6,DATDIR)
 call getarg(7,OUTDIR)

 read (anvrs,*) nvars
 read (apos,*) lpos

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

 STANFI= adate(1:8) //'t2m_s.dat'
 GRIDFI= adate(1:8) //'t2m_g.dat'

 write (abftol,'(I2.2)') int(bftol+0.5)
 write (ascttol,'(I2.2)') int(sqrt(scttol+0.5))
 CSVFI= 'temperatura_' //adate //'_' //trim(ann) //'_bf' //abftol  //'_sct' //ascttol //'.csv'

!-----------------------------------------------------------------------
! DIRECTORY:
! DATDIR:
 FULFISTAG= DATDIR(1:ladir) //STAGFI
 FULFIORO = DATDIR(1:ladir) //OROGRIFI
 FULFIAMI = DATDIR(1:ladir) //'stazioni_areamilano.txt'

! OUTDIR:
 FULFIGR= OUTDIR(1:lodir) //GRIDFI
 FULFIST= OUTDIR(1:lodir) //STANFI
 FULFICSV = OUTDIR(1:lodir) //CSVFI

!-----------------------------------------------------------------------
 if (griglia) write (6,*) 'orography and grid file: ', trim(FULFIORO)
 write (6,*) 'station archive file   : ', trim(FULFISTAG)
 write (6,*) 'MI area station file   : ', trim(FULFIAMI)
 write (6,*) 'observation file       : ', trim(FULFIOB)
 write (6,*) 'station analysis file  : ', trim(FULFIST)
 if (griglia) write (6,*) 'grid analysis file     : ', trim(FULFIGR)
 write (6,*) 'CSV station file       : ', trim(FULFICSV)
 write (6,*) ' '

!-----------------------------------------------------------------------
 write (6,'(a,3f10.2)') 'eps, sigd, sigz', eps, sigd, sigz
 !write (6,*) ' ok?'
 !read (5,'(a1)') ans
 write (6,*) ' '
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
 !write (6,*) ' ok?'
 !read (5,'(a1)') ans

!-----------------------------------------------------------------------
! read station archive:
!-----------------------------------------------------------------------
 write (6,*) 'now read station archive'
 call anaguhird(FULFISTAG,atid,ax,ay,alon,alat,ahgt,auhi,KKA)
 write (6,*) 'archive read, KKA=',KKA
 

!deb
! write (6,'(A4,1X,A8,1X,2A9,2A9,A6)') &
!&   'k','ID   ',' UTMx  ',' UTMy  ',' lon  ',' lat  ',' hgt '
! do k=1,KKA
!  write (6,'(i4,1x,a8,1x,2f9.0,2f9.4,f6.0,f5.0)') &
!&   k,atid(k),ax(k),ay(k),alon(k),alat(k),ahgt(k),(100*auhi(k))
!  !if (mod(k,20) == 0) then
!  ! write (6,*) ' ok?'
!  ! read (5,'(a1)') ans
!  !end if
! end do
 write (6,*) ' '
 call flush(6)
 !write (6,*) ' ok?'
 !read (5,'(a1)') ans

!-----------------------------------------------------------------------
! leggi stazioni area milanese (grande) per CVscore
 write (6,*) ' ora leggo le stazioni dell''area di Milano'
 open (3, file= FULFIAMI)
 k= 0
 do while (.True.)
  read (3,'(a8)', iostat=iaux) amid
  if (iaux < 0) exit
  j= 1
  do while (j <= KKA .and. atid(j) /= amID)
   j= j +1
  end do
  if (j > KKA) then
   write (6,*) 'non trovo in archivio la: ',amid
  else
   k= k +1
   jarkami(k)= j
   ! write (6,*) 'area MI ',k, ' id= ',atid(j),' j=',j
  end if
 end do
 close (3)
 KAMI= k
 write (6,*) ' kami=',KAMI
 !write (6,*) ' ok?'
 !read (5,'(a1)') ans

!-----------------------------------------------------------------------
! compute ZS matrix
!-----------------------------------------------------------------------
 write (6,*) 'building matrix ZS..sigd sigz eps', sigd, sigz, eps
 ZS= 0.d0
 do k=1,KKA
  do k2=1,KKA
   ! horizontal covariance:
   rah= distan(alon(k),alat(k), alon(k2),alat(k2))
   ssh= 0.d0
   arg= rah /sigd
   if (arg < 7.d0) ssh= exp(-arg*arg/2.d0)
   ! vertical covariance:
   rav= abs(ahgt(k)- ahgt(k2))
   ssv= 0.d0
   arg= rav /sigz
   if (arg < 7.d0) ssv= exp(-arg*arg/2.d0)
   ! factorized 3d covariance:
   ZS(k,k2)= ssv*ssh
   ! UHI - decrease correlations between stations in/outside urban areas:
   att= 1.d0 -(1.d0-attmin) *abs(auhi(k)-auhi(k2))
   ZS(k,k2)= ZS(k,k2) *att
  end do
  ZS(k,k)= ZS(k,k) +eps
 end do
 write (6,*) ' ho finito ZS KKA KKAxKKAx8', KKA, KKA*KKA*8
 write (6,*) ' '
 call flush(6)
 !write (6,*) ' ok?'
 !read (5,'(a1)') ans

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
  call GetUTMOrography(FULFIORO, IIG, g4x, g4y, g4lon, g4lat, g4hgt, g4msk,g4uhi)
  write (6,*) 'total number of gridpoints IM IIG: ',IM,IIG
  !write (6,*) ' ok?'
  !read (5,'(a1)') ans

  !-----------------------------------------------------------------------
  ! assign to R8 variables
  i= 0
  do ig=1,IIG
   ! punti interni attivi
   if (g4msk(ig) > 0.d0) then
    i= i +1
    xlon(i)= dble(g4lon(ig))
    xlat(i)= dble(g4lat(ig))
    xhgt(i)= dble(g4hgt(ig))
    xuhi(i)= dble(g4uhi(ig))
   end if
  end do
  II= i
  write (6,*) 'lon lat oro ok'
  write (6,*) 'number of active gridpoints II IEFF: ',II,IEFF
  write (6,*) ' '
  call flush(6)
  !write (6,*) ' ok?'
  !read (5,'(a1)') ans

  !-----------------------------------------------------------------------
  ! compute G matrix
  !-----------------------------------------------------------------------
  write (6,*) 'building G matrix...', sigd, sigz
  ! GG for active gridpoints only (inside administrative boundary)
  do i=1,II
   do k=1,KKA
    rah= distan(xlon(i),xlat(i), alon(k),alat(k))
    ggh= 0.d0
    arg= rah /sigd
    if (arg < 7.d0) ggh= exp(-arg*arg/2.d0)
    ! vertical cov:
    rav= abs(xhgt(i)- ahgt(k))
    ggv= 0.d0
    arg= rav /sigz
    if (arg < 7.d0) ggv= exp(-arg*arg/2.d0)
    ! factorized 3d cov:
    GG(i,k)= ggv* ggh
    ! UHI - decrease correlations between stations in/outside urban areas:
    att= 1.d0 -(1.d0-attmin) *abs(xuhi(i)-auhi(k))
    GG(i,k)= GG(i,k) *att
   end do
  end do
  write (6,*) ' ho finito G II KKA IIxKKAx8=', II, KKA, II*KKA*8
  write (6,*) ' '
  call flush(6)
  !write (6,*) ' ok?'
  !read (5,'(a1)') ans

 end if ! griglia

!-----------------------------------------------------------------------
! inizializza KKP, SRNVP, JARKP, ovvero
! numero di stazioni, matrice inversa, e vettore di indici del passo precedente:
 KKP= 0
 SRNVP= UNDEF
 jarkp= 0

!-------------------------------------------------------------------------------
! this is the time loop:
 write (6,*) 'now start time loop'
 write (6,*) ' '
 call flush(6)

 open (10, file= OUTDIR(1:lodir) //'RMS.txt')
 write (10, '(A5,1X,A12,1X,5A12)') '   n', '  datehour  ' &
&, ' rms(yo-yb)', ' rms(ya-yb)', ' rms(ya-yo)', '  CV-score', ' CVsc-areaMI'
 call flush(10)

 open (12, file= FULFICSV)
 arow= &
& '   Id; Year/Mo/Dy Hr;      yo;       yb; flag;      res;       ya;      yaV;     yIDI;    yIDIV;'
!  12345;12345678901234;12345678;123456789;12345;123456789;123456789;123456789;123456789;123456789;
!  123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456
! totale 96 caratteri

 write (12,'(a96)') arow(1:96)
 call flush(12)

 ! read /write record number:
 norec= 0 ! obsread
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
! nome file per output ascii
  FULFIASCIIa = OUTDIR(1:lodir) // 't2m_ana_32632_' // datehr(n)(1:10) // '.txt'
  FULFIASCIIb = OUTDIR(1:lodir) // 't2m_bkg_32632_' // datehr(n)(1:10) // '.txt'
  FULFIASCIIi = OUTDIR(1:lodir) // 't2m_idi_32632_' // datehr(n)(1:10) // '.txt'
  write (6,'(A)') 'Nomi file ascii: ',FULFIASCIIa,FULFIASCIIb,FULFIASCIIi
  !write (6,'(A)') ' ok?'
  !read (5,'(a1)') ans

!-----------------------------------------------------------------------
! data-ora per file CSV:
  arow= ''
  arow(6:7)= '; '
  write (arow(8:21),'(a)') datehr(n)(1:4) //'/' //datehr(n)(5:6) &
& //'/' //datehr(n)(7:8) //' ' //datehr(n)(9:10) //';'

  !-----------------------------------------------------------------------
  ! INPUT observations (STID yo):
  write (6,*) 'now read observations...'
  stid= ''
  so= UNDEF

  call obsread(FULFIOB, norec, KK, stid, so, nvars, lpos)
  write (6,*) 'back from obsread', KK
  write (6,*) ' '
  call flush(6)

  !deb
  do k=1,KK
   write (6,'(a,i5,a1,a8,a1,f12.2)') ' k id yo', k,' ',stid(k),' ',so(k)
!   if (mod(k,20) == 0) then
!    write (6,*) ' ok?'
!    read (5,'(a1)') ans
!   end if
  end do
!  write (6,*) ' ok?'
!  read (5,'(a1)') ans

  !-----------------------------------------------------------------------
  ! index of station k in archive
  ! data not present in archive are eliminated
  write (6,*) ' assegno jark e ao'
  jark= 0
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
! per caso 02 nov 2009:
!&   .or. stid(k) == 'Lo00041 ' &
! per caso foehn:
!&   .or. stid(k) == 'Lo00649 ' &
!&   .or. (stid(k) == 'Lo00012 ' .and.n >= 12 .and.n <= 17) &
! per caso nebbia:
!&   .or. stid(k) == 'Lo00649 ' &
!     &   .or. (stid(k) == 'Lo00503 ' .and. n >= 18)
!&      ) then
!diagnostica:
!    if (stid(k) == 'Lo00012 ') write (6,*) 'blacklist  012***',atid(j),n
!    if (stid(k) == 'Lo00503 ') write (6,*) 'blacklist  503***',atid(j),n
!    if (stid(k) == 'Lo00649 ') write (6,*) 'blacklist  649***',atid(j)
!    if (stid(k) == 'Lo00041 ') write (6,*) 'blacklist   41***',atid(j), n

    write (6,*) 'stazione fuori archivio: ',k,' ', stid(k)
    KK= KK -1
    do m=k,KK
     stid(m)= stid(m+1)
     so(m)= so(m+1)
    end do
   else

    !write (6,'(i4,1x,a8,a1,a8,1x,i4)') k,stid(k),'=',atid(jark(k)),jark(k)

    jark(k)= j   ! index array
    ao(j)= so(k) ! observed value

    k= k +1
   end if
   !if (mod(k,20) == 0) then
   ! write (6,*) ' ok?'
   ! read (5,'(a1)') ans
   !end if
  end do

  write (6,*) 'assegnato jark'
  write (6,*) ' '
  call flush(6)
  !write (6,*) ' ok?'
  !read (5,'(a1)') ans

!-------------------------------------------------------------------------------
! quality control loop starts here
!-------------------------------------------------------------------------------
  nsca= 0
  abfres= UNDEF
  asctres= UNDEF
  aflag= ' null'
  do k=1,KK
   aflag(jark(k))= ' pass'
  end do

  do while (KK >= KKMIN)
    write (6,*) ' start new QC BFcheck/SCT loop, KK=', KK
    write (6,*) ' '
    call flush(6)

    write (6,*) ' valori prima di calcolare il BF:'
    write (6,'(2a4,1x,a8,1x,4a10,a7,a8)') &
&    'k','j','id','x','y','lon','lat','hgt','yo'
    !write (6,*) ' ok?'
    !read (5,'(a1)') ans

   ! assign variables for background field:
   so= UNDEF
   slon= UNDEF
   slat= UNDEF
   shgt= UNDEF
   do k=1,KK
    j= jark(k)
    ! auxiliary arrays for the background field:
    so(k)= ao(j)
    slon(k)= alon(j)
    slat(k)= alat(j)
    shgt(k)= ahgt(j)

    write (6,'(2i4,1x,a8,1x,2f10.0,2f10.5,f7.0,f8.1)') &
&    k, j, atid(j), ax(j), ay(j), alon(j), alat(j), ahgt(j), ao(j)
    !if (mod(k,20) == 0) then
    ! write (6,*) ' ok?'
    ! read (5,'(a1)') ans
    !end if
   end do
   write (6,*) ' '
   call flush(6)
   !write (6,*) ' ok?'
   !read (5,'(a1)') ans

   !-----------------------------------------------------------------------
   ! background field
   !-----------------------------------------------------------------------

   ! inversione XZ
   call xzinv(KK, so, slon, shgt, avlon, dzbf, zinv, yinv, alpa, gama, alpb, gamb, res)
   avlat=0.d0
   beta= 0.d0
   betb= 0.d0
   if (gama > 0.) then
    write (6,*) 'SITUAZIONE COMPLESSA!'
   end if
   if (gama == gamb) then
    write (6,*) 'senza inversione res=', res
    write (6,*) 'avZ avT', zinv, yinv
    write (6,'(A,3F16.7)') 'Tx Ty Tz', alpa, beta, gama
   else if (gama*gamb < 0.) then
    write (6,*) 'inversione res=', res
    write (6,*) 'Zinv Tinv', zinv, yinv
    write (6,'(A,3F16.7)') 'above Tx Ty Tz', alpa, beta, gama
    write (6,'(A,3F16.7)') 'below Tx Ty Tz', alpb, betb, gamb
   else
    write (6,*) 'cambio pendenza res=', res
    write (6,*) 'Zinv Tinv', zinv, yinv
    write (6,'(A,3F16.7)') 'above Tx Ty Tz', alpa, beta, gama
    write (6,'(A,3F16.7)') 'below Tx Ty Tz', alpb, betb, gamb
   end if

   !-----------------------------------------------------------------------
   ! estimate bf on station points
   ab= UNDEF

   zabov= zinv +dzbf
   zbelo= zinv -dzbf

   do j=1,KKA
    if (ahgt(j) > zabov) then
     ab(j)= alpa *(alon(j) -avlon) &
&          +beta *(alat(j) -avlat) &
&          +gama *(ahgt(j) -zinv) &
&          +yinv
    else if (ahgt(j) < zbelo) then
     ab(j)= alpb *(alon(j) -avlon) &
&          +betb *(alat(j) -avlat) &
&          +gamb *(ahgt(j) -zinv) &
&          +yinv
    else
     bfabov= alpa *(alon(j) -avlon) &
&           +beta *(alat(j) -avlat) &
&           +gama *(zabov -zinv) +yinv
     bfbelo= alpb *(alon(j) -avlon) &
&           +betb *(alat(j) -avlat) &
&           +gamb *(zbelo -zinv) +yinv
     ab(j)= (bfabov *(ahgt(j) -zbelo) +bfbelo *(zabov -ahgt(j))) &
&           /(zabov -zbelo)
    end if
   end do
   write (6,*) ' '
   call flush(6)
   !write (6,*) ' ok?'
   !read (5,'(a1)') ans

   !-----------------------------------------------------------------------
   ! BF check
   write (6,*) 'BF check, tol=', bftol, '°C'
   absdmx= 0.d0
   kdmx= 0
   do k=1,KK
    j= jark(k)

    absd= abs(ao(j)-ab(j))
    abfres(j)= absd ! comunque

    if (absd > absdmx) then
     absdmx= absd
     kdmx= k
    end if
   end do

   if (absdmx > bftol) then
    ! SCARTO 
    nsca= nsca +1
    j= jark(kdmx)
    aflag(j)= '   BF'
    write (6,'(a,2i5,1x,a8,1x,f6.0,3f7.1)') 'BF check failed nsca j id h yo yb |yo-yb|: ' &
&        , nsca, j, atid(j), ahgt(j), ao(j), ab(j), abfres(j)

    ! LO ELIMINO: 
    KK= KK -1
    write (6,*) 'bfcheck: KDMX', kdmx
    do k=kdmx,KK
     jark(k)= jark(k+1)
    end do
    jark(KK+1)= 0
    write (6,*) 'bfcheck: KK=',KK, 'now cycle'
    write (6,*) ' '
    call flush(6)
    !write (6,*) ' ok?'
    !read (5,'(a1)') ans

    cycle

   else
    j= jark(kdmx)
    write (6,'(a,i3,1x,a8,1x,f6.0,3f10.2)') 'BF check passed: nsca id h yo yb res ' &
&    , nsca, atid(j), ahgt(j), ao(j), ab(j), abfres(j)
    write (6,*) ' '
    call flush(6)
    !write (6,*) ' ok?'
    !read (5,'(a1)') ans
   end if

   write (6,*) 'BF check done, KK=', KK
   write (6,*) 'res=', res, sqrt(res/KK)
   write (6,*) ' '

   write (6,'(2a4,1x,a8,1x,4a10,a7,2a8)') &
&   'k','j','id ','x','y','lon','lat','hgt','yo','yb'
   do k=1,KK
    j= jark(k)
    write (6,'(2i4,1x,a8,1x,2f10.0,2f10.5,f7.0,f8.1,f8.2)') &
&    k, j, atid(j), ax(j), ay(j), alon(j), alat(j), ahgt(j), ao(j), ab(j)
   end do
   write (6,*) 'dopo check BF, KK=', KK
   write (6,*) ' '
   call flush(6)
   !write (6,*) ' ok?'
   !read (5,'(a1)') ans

   !-----------------------------------------------------------------------
   ! update matrix (S+R)^-1
   write (6,*) 'updating inverse matrix (S+R)^-1 ...', KKP, '-->', KK
   !write (6,*) ' ok?'
   !read (5,'(a1)') ans
   ! reset SRNV (not the previous SRNVP):
   SRNV= UNDEF
   call udsrnv(KTOT, KKP, jarkp, SRNVP, KK, jark, ZS, SRNV)
   write (6,*) 'back from updating inverse matrix'
   write (6,*) ' '
   call flush(6)
   !write (6,*) ' ok?'
   !read (5,'(a1)') ans

   !-----------------------------------------------------------------------
   ! elementi diagonali di W=HK
   write (6,*) 'calcolo diag(W), eps=', eps
   dialmx= 0.d0
   diasmn= 1.d10
   dial= UNDEF
   diasrnv= UNDEF
   do k=1,KK
    ! Wkk:
    dial(k)= 1.d0 -eps *SRNV(k,k)
    ! zk:
    diasrnv(k)= SRNV(k,k)
    if (dial(k) > dialmx) dialmx= dial(k)
    if (diasrnv(k) < diasmn) diasmn= diasrnv(k)
   end do
   ! write (6,'(5f16.8)') (dial(m),m=1,KK)
   write (6,*) 'dial max=', dialmx
   ! write (6,*) ' ok?'
   ! read (5,'(a1)') ans
   ! write (6,*) 'diag(SRNV)= '
   ! write (6,'(5f16.8)') (diasrnv(m),m=1,KK)
   write (6,*) 'diasrnv min',diasmn
   write (6,*) ' '
   call flush(6)
   ! write (6,*) ' ok?'
   ! read (5,'(a1)') ans

   !-----------------------------------------------------------------------
   ! innovation
   write (6,*) 'compute innovation'
   ydmax= 0.d0
   res= 0.d0
   yd= UNDEF
   didi= UNDEF
   do k=1,KK
    yd(k)= ao(jark(k))- ab(jark(k))

    if (abs(yd(k)) > ydmax) ydmax= abs(yd(k))
    res= res+ yd(k)*yd(k)
    didi(k)= 1.d0
   end do
   if (KK > 0) then
    res= sqrt(res/KK)
   else
    res= UNDEF
   end if
   write (6,*) 'rms & max (yo-yb)=', res, ydmax
   write (6,*) ' '
   call flush(6)
   !write (6,*) ' ok?'
   !read (5,'(a1)') ans

   !-----------------------------------------------------------------------
   ! compute YZ and ZIDI
   ! yz
   write (6,*) 'now compute yz=(S+R)^-1*d '
   az= UNDEF
   do k=1,KK
    j= jark(k)
    az(j)= 0.d0
    do k2=1,KK
     az(j)= az(j) +SRNV(k,k2) *yd(k2)
    end do
   end do
   write (6,*) 'yz ok ',KK

  !zidi
   zidi= UNDEF
   write (6,*) 'now compute zidi=(S+R)^-1*[1 1...1] '
   do k=1,KK
    zidi(k)= 0.d0
    do k2=1,KK
     zidi(k)= zidi(k) +SRNV(k,k2) *didi(k2)
    end do
   end do
   write (6,*) 'zidi ok '
   write (6,*) ' '
   call flush(6)
   !write (6,*) ' ok?'
   !read (5,'(a1)') ans

!-------------------------------------------------------------------------------
! compute analysis on obs sites
!-------------------------------------------------------------------------------
  ! stimo comunque su TUTTE le stazioni dell'archivio:
   write (6,*)'analysis on obs sites...'
   aa= UNDEF
   aidi= UNDEF
   do j=1,KKA
    dy= 0.d0
    zum= 0.d0
    do k2=1,KK
     j2= jark(k2)
     dy= dy+ ZS(j,j2)* az(j2)
     zum= zum+ ZS(j,j2)* zidi(k2)
     ! sulla diagonale di ZS c'e' 1+eps, devo togliere eps*yz
     if (j2 == j) then
      dy=  dy -eps *az(j2)
      zum= zum -eps *zidi(k2) 
     end if
    end do
    ! add bf:
    aa(j)= ab(j) +dy
    aidi(j)= zum
   end do
   ! residui:
   qres= 0.d0
   pres= 0.d0
   do k=1,KK
    j= jark(k)
    qres= qres+ (aa(j)-ab(j))**2
    pres= pres+ (aa(j)-ao(j))**2
   end do
   qres= sqrt(qres/KK)
   pres= sqrt(pres/KK)
   write (6,*) 'rms(yo-yb) ',  res
   write (6,*) 'rms(ya-yb) ', qres
   write (6,*) 'rms(ya-yo) ', pres
   write (6,*) ' '
   call flush(6)
   !write (6,*) ' ok?'
   !read (5,'(a1)') ans

!-------------------------------------------------------------------------------
! Cross Validation
   write (6,*) 'computing CV analysis...'
   aav= aa
   aidiv= aidi
   adial= UNDEF
   anvdiasrnv= UNDEF
   do k=1,KK
    j= jark(k)
    ! tre versioni della stessa cosa:
    ! 1:
    ! aav(j) =  aa(j) -(dial(k)/(1.d0 -dial(k))) *(ao(j) -aa(j))
    ! aidiv(j)= aidi(j) -(dial(k)/(1.d0 -dial(k))) *(1.d0 -aidi(j))
    ! 2:
    ! aav(j) =  (aa(j) -dial(k) *ao(j)) /(1.d0 -dial(k))
    ! aidiv(j)= (aidi(j) -dial(k)) /(1.d0 -dial(k))
    ! 3:
    aav(j) = ao(j) +(aa(j) -ao(j)) /(1.d0 -dial(k))
    aidiv(j)=  1.d0 +(aidi(j) -1.d0) /(1.d0 -dial(k))

    adial(j)= dial(k)
    anvdiasrnv(j)= 1.d0 /diasrnv(k)
   end do
   write (6,*) 'CV analysis computed'
   write (6,*) ' '
   call flush(6)
   !write (6,*) ' ok?'
   !read (5,'(a1)') ans

!-----------------------------------------------------------------------
! SCT (CV check)
!-----------------------------------------------------------------------
   write (6,*) 'SCT: tol=', scttol, '°C^2 ', sqrt(scttol), '°C'
   absdmx= 0.d0
   kdmx= 0
   do k=1,KK
    j=jark(k)

    absd= (ao(j) -aa(j)) *(ao(j) -aav(j))
    asctres(j)= sqrt(absd)

    if (absd > absdmx) then
     absdmx= absd
     kdmx= k
    end if
   end do

   if (absdmx > scttol) then
    ! SCARTO 
    nsca= nsca +1
    j= jark(kdmx)
    aflag(j)= '  SCT'
    write (6,'(a,2i5,1x,a8,1x,f6.0,4f7.1)') 'SCT failed: nsca j id hgt tdo tda tdav res: ' &
&     , nsca, j, atid(j), ahgt(j), ao(j), aa(j), aav(j), asctres(j)

    ! LO ELIMINO: 
    KK= KK -1
    do k=kdmx,KK
     jark(k)= jark(k+1)
    end do
    jark(KK+1)= 0
    write (6,*) 'SCT: KK=',KK
    write (6,*) ' '
    call flush(6)
    !write (6,*) ' ok?'
    !read (5,'(a1)') ans

   else
    j= jark(kdmx)
    write (6,'(a,i5,1x,a8,5f12.3)') 'SCT passed: nsca id h yo ya yav ' &
&     , nsca, atid(j), ahgt(j), ao(j), aa(j), aav(j), asctres(j)
    exit
   endif

  !-----------------------------------------------------------------------
  ! end do-while BFcheck and SCT :
  end do ! while (KK > KKMIN)
  !-----------------------------------------------------------------------
  write (6,*) ' '
  call flush(6)

  !-----------------------------------------------------------------------
  if (KK >= KKMIN) then
   write (6,*) 'BF check and SCT passed, KK=', KK
   write (6,'(a,i2)') 'numero di scarti: ',nsca
   write (6,'(a,a8,1x,f6.0,3f7.1)') 'prima stazione non scartata: ' &
&  , atid(j),ahgt(j),ao(j),aa(j),aav(j)

  else
   ab= UNDEF
   aa= UNDEF
   aav= UNDEF
   aidi= 0.d0
   aidiv= 0.d0
   res= UNDEF
   qres= UNDEF
   pres= UNDEF
   write (6,'(2(a,i4))') ' KK=', KK, ' < KKMIN= ',KKMIN
   write (6,*) ' analisi <-- UNDEF'
  end if

  write (6,*) ' '
  call flush(6)
  !write (6,*) ' ok?'
  !read (5,'(a1)') ans

  !-----------------------------------------------------------------------
  ! output CSV station file
! arow= &
!&'   Id; Year/Mo/Dy Hr;      yo;       yb; flag;      res;       ya;      yaV;     yIDI;    yIDIV;'
!  12345;12345678901234;12345678;123456789;12345;123456789;123456789;123456789;123456789;123456789;
!  123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456
!           1         2         3         4         5         6         7         8         9
  do j=1,KKA
   arow(1:6)= '     ;'
   arow(1:5)= atid(j)(3:7)

   arow(22:30)= '        ;'
   write (arow(22:29),'(f8.1)') ao(j)
   if (ao(j) == UNDEF)  write (arow(22:29),'(i8)') int(UNDEF)

   arow(31:40)= '         ;'
   write (arow(31:39),'(f9.2)') ab(j)
   if (ab(j) == UNDEF)  write (arow(31:39),'(i9)') int(UNDEF)

   arow(41:46)= '     ;'
   arow(41:45)= aflag(j)

   arow(47:56)= '         ;'
   if (aflag(j) == ' null') then
    write (arow(47:55),'(i9)') int(UNDEF)
   else if (aflag(j) == '   BF') then
    write (arow(47:55),'(f9.2)') abfres(j)
   else
    write (arow(47:55),'(f9.2)') asctres(j)
   end if

   arow(57:66)= '         ;'
   write (arow(57:65),'(f9.2)') aa(j)
   if (aa(j) == UNDEF)  write (arow(57:65),'(i9)') int(UNDEF)

   arow(67:76)= '         ;'
   write (arow(67:75),'(f9.2)') aav(j)
   if (aav(j) == UNDEF)  write (arow(67:75),'(i9)') int(UNDEF)

   arow(77:86)= '         ;'
   write (arow(77:85),'(f9.4)') aidi(j)
   if (aidi(j) == UNDEF)  write (arow(77:85),'(i9)') int(UNDEF)

   arow(87:96)= '         ;'
   write (arow(87:95),'(f9.4)') aidiv(j)
   if (aidiv(j) == UNDEF)  write (arow(87:95),'(i9)') int(UNDEF)

   write (12,'(a96)') arow(1:96)
   call flush(12)
  end do ! j
  write (6,*) 'analysis etc. written on CSV station file'
  write (6,*) ' '
  call flush(6)
  !write (6,*) ' ok?'
  !read (5,'(a1)') ans

  !-----------------------------------------------------------------------
  ! CV-score
  !-----------------------------------------------------------------------
  if (KK >= KKMIN) then
   cvscore= 0.d0
   do k=1,KK
    j= jark(k)
    cvscore= cvscore +(ao(j) -aav(j))**2
   end do
   cvscore= sqrt(cvscore /KK)
  else
   cvscore= UNDEF
  endif

  ! CV-score area MI:
  cvskami= 0.d0
  koami= 0
  do k=1,KAMI
   j= jarkami(k)
   if (aflag(j)/=' pass' .or. aav(j)==UNDEF) cycle
   koami= koami +1
   cvskami= cvskami +(ao(j) -aav(j))**2
  end do
  if (koami >= 1) then
   cvskami= sqrt(cvskami /koami)
  else
   cvskami= UNDEF
  end if

  write (6,*) 'cvskami KOAMI KAMI n',cvskami,koami, KAMI, n
  write (6,*) ' '
  call flush(6)

  ! write CVscore CVskami:
  if (KK >= KKMIN .and. cvskami/=UNDEF) then
   write (10, '(i5,1x,a12,1x,5f12.4)') n, datehr(n),  res, qres, pres, cvscore, cvskami
   write (6, '(a,i5,1x,a12,1x,5f12.4)') 'CVscore: ', n, datehr(n), res, qres, pres, cvscore, cvskami
  else if (KK >= KKMIN) then  ! aMI UNDEF:
   write (10, '(i5,1x,a12,1x,4f12.4,i12)') n, datehr(n),  res, qres, pres, cvscore, int(UNDEF)
   write (6, '(a,i5,1x,a12,1x,4f12.4,f12.0)') 'CVscore: ', n, datehr(n), res, qres, pres, cvscore &
&   , cvskami
  else ! tutti UNDEF:
   write (10, '(i5,1x,a12,1x,5i12)') n, datehr(n), int(UNDEF), int(UNDEF), int(UNDEF)&
&   , int(UNDEF), int(UNDEF)
   write (6, '(a,i5,1x,a12,1x,5f12.0)') 'CVscore: ', n, datehr(n), res, qres, pres, cvscore, cvskami
  end if
  call flush(10)

  write (6,*) ' '
  call flush(6)
  !write (6,*) ' ok?'
  !read (5,'(a1)') ans

  !-----------------------------------------------------------------------------
  ! GRID
  if (griglia .and. KK>=KKMIN) then

   !-----------------------------------------------------------------------
   ! BF on grid points
   write (6,*) 'estimating BF on grid points...'
   xb= UNDEF
   write (6,*) 'Zinv Tinv zabov', zinv, yinv, zabov
   write (6,'(A,3F16.7)') 'above Tx Ty Tz', alpa, beta, gama
   write (6,'(A,3F16.7)') 'below Tx Ty Tz', alpb, betb, gamb
   do i=1,II
    if (xhgt(i) > zabov) then
     xb(i)= alpa *(xlon(i) -avlon) &
&          +beta *(xlat(i) -avlat) &
&          +gama *(xhgt(i) -zinv) &
&          +yinv
    else if (xhgt(i) < zbelo) then
     xb(i)= alpb *(xlon(i) -avlon) &
&          +betb *(xlat(i) -avlat) &
&          +gamb *(xhgt(i) -zinv) &
&          +yinv
    else
     bfabov= alpa *(xlon(i) -avlon) &
&           +beta *(xlat(i) -avlat) &
&           +gama *(zabov -zinv) +yinv
     bfbelo= alpb *(xlon(i) -avlon) &
&           +betb *(xlat(i) -avlat) &
&           +gamb *(zbelo -zinv) +yinv
     xb(i)= (bfabov *(xhgt(i) -zbelo) +bfbelo *(zabov -xhgt(i))) &
&           /(zabov -zbelo)
    end if
   end do
   write (6,*) ' ok bf'
   write (6,*) ' '
   call flush(6)
   !write (6,*) ' ok?'
   !read (5,'(a1)') ans

   !-----------------------------------------------------------------------
   ! analysis estimate on grid
   !-----------------------------------------------------------------------
   xa= UNDEF
   write (6,*) ' computing analysis on grid x=G*z...', datehr(n)
   write (6,*) ' calcolo G*z...', KKA, KK, II

   res= 0.d0
   do i=1,II ! GG for active gridpoints only
    dx(i)= 0.d0
    xum= 0.d0
    do k=1,KK
     j= jark(k)
     dx(i)= dx(i) +GG(i,j) *az(j)
     xum= xum+ GG(i,j)*zidi(k)
    end do
    res= res+ dx(i)**2
    xa(i)= xb(i)+ dx(i)
    xidi(i)= xum
   end do
   close (3)
   res= sqrt(res/II)

   write (6,*) ' analysis on grid ok, rms(x-xb) ', res
   write (6,'(a,5f14.6)') 'xb ', xb(1), xb(II/3), xb(II/2), xb(2*II/3), xb(II)
   write (6,'(a,5f14.6)') 'xa ', xa(1), xa(II/3), xa(II/2), xa(2*II/3), xa(II)
   write (6,*) ' '
   call flush(6)
   !write (6,*) ' ok?'
   !read (5,'(a1)') ans

  !-----------------------------------------------------------------------
  else if (griglia) then !  KK < KKMIN
   ! UNDEF on grid:
   xb= UNDEF
   xa= UNDEF
   xidi= 0.d0
   write (6,*) '*** zero obs (or fewer than acceptable): ', datehr(n), KK
   call flush(6)
   !write (6,*) ' ok?'
   !read (5,'(a1)') ans
  end if ! griglia .and. KK>=KKMIN

  !-----------------------------------------------------------------------
  ! output obs analysis estimate
  !-----------------------------------------------------------------------
  ! station GrADS output
  !  do k=1,KK
  !   j= jark(k)
  !   stid(k)= atid(j)
  !   s4x(k)= ax(j)
  !   s4y(k)= ay(j)
  !   shgt(k)= ahgt(j)
  !   so(k)= ao(j)
  !   sb(k)= ab(j)
  !   sa(k)= aa(j)
  !   sav(k)= aav(j)
  !  end do
  ! le salvo tutte con gli UNDEF:
  stid= ''
  s4x= UNDEF
  s4y= UNDEF
  shgt= UNDEF
  suhi= UNDEF
  so= UNDEF
  sb= UNDEF
  sa= UNDEF
  sav= UNDEF
  do j=1,KKA
   stid(j)= atid(j)
   s4x(j)= ax(j)
   s4y(j)= ay(j)
   shgt(j)= ahgt(j)
   suhi(j)= auhi(j)
   so(j)= ao(j)
   sb(j)= ab(j)
   sa(j)= aa(j)
   sav(j)= aav(j)
  end do

  write (6,*) 'writing obs analysis...', KKA, FULFIST
  call stzuhiwri(FULFIST, nwrec, KKA, stid, s4x, s4y, shgt, so, sb, sa, sav, suhi)

  write (6,*) ' analysis etc. written on GrADS station file '
  write (6,*) ' '
  call flush(6)
  !write (6,*) ' ok?'
  !read (5,'(a1)') ans

  !-----------------------------------------------------------------------------
  if (griglia) then
   ! output analysis on grid
   i= 0
   do ig=1,IIG
    gxb(ig)= sngl(UNDEF)
    gxa(ig)= sngl(UNDEF)
    gidi(ig)= sngl(UNDEF)
    if (g4msk(ig) > 0.d0) then
     i=i +1
     gxb(ig)= sngl(xb(i))
     gxa(ig)= sngl(xa(i))
     gidi(ig)= sngl(xidi(i))
    end if
   end do
   write (6,*) 'writing grid analysis...', i, FULFIGR
   call griwriUTM(n, FULFIGR, gxb, gxa, gidi)
   write (6,*) ' grid analysis etc. written on binary file '
   write (6,*) ' '
! scrivo su esri-asci per import rasdaman, MS   
   call asciiUTM( FULFIASCIIb, gxb)
   call asciiUTM(FULFIASCIIa, gxa)
   call asciiUTM(FULFIASCIIi, gidi)
   write (6,*) ' grid analysis etc. written on ascii file '
   write (6,*) ' '

   call flush(6)
   !write (6,*) ' ok?'
   !read (5,'(a1)') ans

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
  !write (6,*) ' ok?'
  !read (5,'(a1)') ans

!-----------------------------------------------------------------------
 end do ! end do n
!-----------------------------------------------------------------------

 !-----------------------------------------------------------------------
 ! close CVscore file:
 close (10)
 ! close file CSV:
 close (12)

 !-------------------------------------------------------------------------------
 stop
end program t2m19
!*******************************************************************************
! SUBROUTINES AND FUNCTIONS (here)
! anaguhird
! obsread
! stzuhiwri
!***************************************************************************************************
subroutine anaguhird(STAGFI,stid,sx,sy,slon,slat,shgt,suhi, KKA)
!-------------------------------------------------------------------------------
! read station archive from GrADS "station" file
! compute coordinates UTM->lon,lat
!
! header 28 bytes= 7 *4bytes
! 1(1,2):  8 STID     CHARACTER*8 equiv BID(2)*4
! 2 ( 3): 12 LAT/Y    REAL(4)
! 3 ( 4): 16 LON/X    REAL(4)
! 4 ( 5): 20 TIM      REAL(4)
! 5 ( 6): 24 NLEV     INTEGER(4)
! 6 ( 7): 28 NFLAG    INTEGER(4)
! after the header, NVAR REAL(4) variables: NVAR*4bytes:  REAL(4) PUFF(1:NVARS)
!
! tutte le stazioni in anagrafica
!-------------------------------------------------------------------------------
 use coords
 use ll_utm
 implicit none

 integer, parameter :: KTOT= 1000

! UNDEF here is different from UNDEF in the main code!
 real(4), parameter :: UNDEF= -999.9

 character(200) :: STAGFI
 character(8), dimension(KTOT) ::  stid
 real(8), dimension(KTOT) :: sx, sy, slon, slat, shgt, suhi
 integer KKA

!-----------------------------------------------------------------------
! number of variables:
 integer, parameter :: NVARS= 6
 real(4), dimension(NVARS) :: puff

!-----------------------------------------------------------------------
 real(4) :: tim, r4lox, r4lay
 integer nlev, nflag
 character(8) :: aid

!-----------------------------------------------------------------------
 integer :: irec, j, k
 integer :: lstrim, lstr, iid

!-----------------------------------------------------------------------
 type(GEO_COORD)   :: tLL
 type(UTM_COORD)   :: tUTM  ! tUTM%Zone tUTM%N tUTM%E (c4 r8 r8)

!-----------------------------------------------------------------------
 character(1) :: ans

!===============================================================================
 ! init -> UNDEF
 stid= ' '
 sx  = UNDEF
 sy  = UNDEF
 slon= UNDEF
 slat= UNDEF
 shgt= UNDEF
 suhi= 0.

 !-----------------------------------------------------------------------
 write (6,*) 'anaguhird: file: ',trim(STAGFI)
! write (6,*) ' ok?'
! read (5,'(a1)') ans
 open (8, file= STAGFI, form= 'UNFORMATTED', access= 'DIRECT', recl=4)
 irec= 0
 k= 1

 ! read first header:
 irec= irec +1
 read (8,rec=irec) aid(1:4)
 irec= irec +1
 read (8,rec=irec) aid(5:8)
 irec= irec +1
 read (8,rec=irec) r4lay
 irec= irec +1
 read (8,rec=irec) r4lox
 irec= irec +1
 read (8,rec=irec) tim
 irec= irec +1
 read (8,rec=irec) nlev
 irec= irec +1
 read (8,rec=irec) nflag
 
 !write(6,*) aid,r4lay,r4lox,tim,nlev,nflag

 !-----------------------------------------------------------------------
 ! check for record terminator:
 do while (nlev /= 0)

  ! read variable codes:
  do j=1,NVARS
   irec= irec +1
   read (8,rec=irec) puff(j)
  end do
!  write(6,*) puff
  ! assign variables:
  ! STation ID:
  stid(k)= aid
!  lstr= lstrim(aid)
!  write (6,*) lstr, iid
!  read (aid(1:lstr),*) iid
!  write (stid(k)(3:7),'(i5.5)') aid(3:7)
!  write (6,'(a1,a8,a2,i4)') '[',stid(k),']=',iid
!  if (mod(k,20) == 0) then
!    write (6,*) ' ok?'
!    read (5,'(a1)') ans
!  end if

  sy(k)= dble(r4lay)
  sx(k)= dble(r4lox)

  shgt(k) = dble(puff(1))
  if (puff(2)/=UNDEF) suhi(k) = dble(puff(2))

  !------------------------------------------------------------
  ! codici sensori nel file GrADS:
  !    VARS  6 
  ! 01 elev    0  99  elevazione ASL [m]
  ! 02 wlu     0  99  peso UHI
  ! 03 rete    0  99  codice rete
  ! 04 temp    0  99  termometro
  ! 05 pluvio  0  99  pluviometro
  ! 06 ur      0  99  umidita relativa
  !------------------------------------------------------------

  ! next obs:
  k= k +1

  ! read next header:
  irec= irec +1
  read (8,rec=irec) aid(1:4)
  irec= irec +1
  read (8,rec=irec) aid(5:8)
  irec= irec +1
  read (8,rec=irec) r4lay
  irec= irec +1
  read (8,rec=irec) r4lox
  irec= irec +1
  read (8,rec=irec) tim
  irec= irec +1
  read (8,rec=irec) nlev
  irec= irec +1
  read (8,rec=irec) nflag

  ! go and check for record terminator:
 end do

 ! record terminator found; number of stations:
 KKA= k -1
 close (8)
 write (6,*) 'anaguhird: read KKA=',KKA
! write (6,*) ' ok?'
! read (5,'(a1)') ans

 !-----------------------------------------------------------------------
 do k=1,KKA
 ! coordinate conversion UTM-> lon,lat and variable choice

  tUTM%E = sx(K)
  tUTM%N = sy(k)
  tUTM%Zone = '+32T'
!
  call UTMtoLL(INTERNATIONAL, tUTM, tLL)
!
  slon(k) = tLL%Lon
  slat(k) = tLL%Lat
 end do

!-------------------------------------------------------------------------------
 do k=1, KKA
  write (6,'(a,i5,1x,a8,1x,2f8.2,f8.0,f8.2)') &
&  ' k id lon lat hgt var',k,stid(k),slon(k),slat(k),shgt(k),suhi(k)
!  if (mod(k,20) == 0) then
!   write (6,*) ' ok?'
!   read (5,'(a1)') ans
!  end if
 end do
 write (6,*) 'anaguhird: return'
! write (6,*) ' ok?'
! read (5,'(a1)') ans

!-------------------------------------------------------------------------------
 return
end subroutine anaguhird
!***************************************************************************************************
subroutine obsread(FULFIOB, nrec, KK, stid, yo, nvars, lpos)
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
 character(200) :: FULFIOB
 !integer, save :: nrec
 integer :: nrec
 integer :: KK
 character(8), dimension (KTOT) :: stid
 real(4), dimension (KTOT) :: yo
 integer :: nvars, lpos

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
 integer :: lstr
! functions:
 integer :: lstrim

 character(1) :: ans

!===============================================================================
! temperature
 yomin= -30.d0
 yomax= +45.d0
! pressure
 ! yomin= 500.d0
 ! yomax= 1200.d0
!-----------------------------------------------------------------------
 stid= 'Lo00000 '
 aid= ''

 write (6,*) 'obsread: file :',trim(FULFIOB)
 write (6,*) 'obsread: nrec :',nrec
 write (6,*) 'obsread: nvars:',nvars
 write (6,*) 'obsread: lpos :',lpos
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

! write (6,*) 'ID: ',aid
! write (6,*) 's4lox: ',s4lox
! write (6,*) 's4lay: ',s4lay
! write (6,*) 'tim: ',tim
! write (6,*) 'nlev: ',nlev
! write (6,*) 'nflag: ',nflag
! write (6,*) ' ok?'
! read (5,'(a1)') ans

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
   write (6,'(a,a8,a,f12.3)') 'ID: [' , aid ,'] obs=', puff(lpos)

  ! STation ID:
!  stid(k)= 'Lo00000 '
  stid(k)= aid

!  lstr= lstrim(aid)
!  write (6,'(a1,a4,a1,a4)') ' ',bid(1),':', bid(2)
!  write (6,'(a1,a8,a1,i3)') '[',aid,']', lstr
!  read (aid(1:lstr),*) iid
!  write (stid(k)(3:7),'(i5.5)') iid
  !write (6,*) 'iid=',iid
  !write (6,'(a1,a8,a2,i4)') '[',stid(k),']=', iid

  ! observation:
  yo(k)= puff(lpos)

!  write (6,'(a,4f10.2)') 'obsread: puff', puff(1), puff(2), puff(3), puff(4)
!  write (6,*) 'obsread: puff ', puff
!  write (6,'(a,i5,1x,a8,f8.1)') 'obsread: k id yo', k, stid(k), yo(k)
!  write (6,*) ' ok?'
!  read (5,'(a1)') ans

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
 write (6,*) 'obsread: read KK=',KK

!-----------------------------------------------------------------------
! check for missing data and climatological range
! check on station elevation is done in ANAGREAD

  write (6,*) 'obsread: missing data check'
  write (6,*) '         and climatological range check min-max: ', yomin, yomax
 k= 0
 do while (k < KK)
  k= k +1

  if (yo(k) < yomin .or. yo(k) > yomax) then
   ! output:
   if (yo(k) /= UNDEF) write (6,*) 'obsread: climat check failed id obs: ', stid(k), yo(k)
   KK= KK -1
   do j=k,KK
    stid(j)= stid(j+1)
    yo(j)= yo(j+1)
   end do
   k= k -1
  end if
 end do
 write (6,*) 'obsread: data climatologically checked and'
 write (6,*) '         missing data eliminated, KK=',KK
 write (6,*) ' '

!-----------------------------------------------------------------------
! zero elsewhere
 do j=KK+1,KTOT
  stid(j)= ' '
  yo(j)= 0.d0
 end do

!-----------------------------------------------------------------------
! do k=1,KK
!   write (6,'(a,i5,1x,a8,1x,f20.6)') 'obsread: k id obs', k, stid(k),yo(k)
! end do
! write (6,*) 'obsread: KK=', KK
! write (6,*) ' ok?'
! read (5,'(a1)') ans

!-------------------------------------------------------------------------------
 return
end subroutine obsread
!*******************************************************************************
subroutine stzuhiwri(FULFIST, nrec, KK, stid, s4x, s4y, shgt, so, sb, sa, sav, suhi)
!-------------------------------------------------------------------------------
! write observations on GrADS "station" file
!-------------------------------------------------------------------------------
 implicit none
 integer, Parameter :: KTOT=1000
 real(4), parameter :: UNDEF=-9999.0

 character*200 FULFIST
 !integer, save :: nrec
 integer :: nrec
 integer :: KK
 character(8), dimension (KTOT) :: stid
 real(4), dimension (KTOT) :: s4x, s4y, shgt, so, sb, sa, sav, suhi

 real(4) :: tim
 integer :: nlev, nflag

 integer k
! character ans*1
!===============================================================================
 open (8, file= FULFIST, form= 'UNFORMATTED', access= 'DIRECT', recl=4)
!-------------------------------------------------------------------------------
! this is the good instant:
 tim= 0.
 nflag= 1
 nlev= 1

 do k=1,KK
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

  nrec= nrec +1
  write (8,rec=nrec) shgt(k)
  nrec= nrec +1
  write (8,rec=nrec) so(k)
  nrec= nrec +1
  write (8,rec=nrec) sb(k)
  nrec= nrec +1
  write (8,rec=nrec) sa(k)
  nrec= nrec +1
  write (8,rec=nrec) sav(k)
  nrec= nrec +1
  write (8,rec=nrec) suhi(k)
 end do

!-------------------------------------------------------------------------------
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

!-------------------------------------------------------------------------------

 close (8)
 write (6,*) 'stzuhiwri: KK and last nrec=', KK, nrec

! write (6,*) 'stzuhiwri: return '
! write (6,*) ' ok?'
! read (5,'(a1)') ans

!-------------------------------------------------------------------------------
 return
end subroutine stzuhiwri
!*******************************************************************************
!*******************************************************************************
subroutine griwriUTM(n, FULFIGR, xb, xa, xidi)
!-------------------------------------------------------------------------------
 implicit none
! integer, parameter :: INX=177, INY=174, IM=INX*INY
 integer, parameter :: INX=251, INY=249, IM=INX*INY

 integer :: n
 character(200) :: FULFIGR
 real(4), dimension(IM) :: xb, xa, xidi

 integer :: nrec
! character(1) :: ans

!===============================================================================
 nrec= (n -1) *3

 WRITE (6,*) 'griwri: n IM nrec ', n, IM, nrec

 open (9, file= FULFIGR, form= 'UNFORMATTED', access= 'DIRECT', recl=IM*4)

! xb, xa, xidi
 nrec= nrec +1
 write (9,rec=nrec) xb
 nrec= nrec +1
 write (9,rec=nrec) xa
 nrec= nrec +1
 write (9,rec=nrec) xidi

 close (9)
 write (6,*) 'griwri: last nrec ', nrec
 !write (6,*) 'ok?'
 !read (5,'(a1)') ans

!-------------------------------------------------------------------------------
 return
end subroutine griwriUTM
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
!******************************************************************************
