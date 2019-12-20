! Francesco Uboldi, 2006-2011
! MS, dicembre 2019 - tutto in UTM
!
! anagreadUTM
! gridread: eliminata, sostituita da GetUTMOrography in SpatialStuff.f90
! stzwri
! griwri: ora per griglia UTM 1 km 251X249
! udmat:
!  udsrnv
!  isortin
!  matsortin
!  invmpj
!  invmmj
! FUNCTIONS
! lstrim
! distan
!***************************************************************************************************
subroutine anagreadUTM(STAGFI,stid,sx,sy,slon,slat,shgt, KKA)
!-------------------------------------------------------------------------------
! read station archive from GrADS "station" file 
! compute coordinates GB->lon,lat
!
! header 28 bytes= 7 *4bytes
! 1(1,2):  8 STID     CHARACTER*8 equiv BID(2)*4
! 2 ( 3): 12 LAT/Y    REAL(4)
! 3 ( 4): 16 LON/X    REAL(4)
! 4 ( 5): 20 TIM      REAL(4)
! 5 ( 6): 24 NLEV     INTEGER(4)
! 6 ( 7): 28 NFLAG    INTEGER(4)
! after the header, NVAR REAL(4) variables: NVAR*4bytes:  REAL(4) PUFF(1:NVARS)
!!------------------------------------------------------------
! codici sensori nel file GrADS:
!VARS 6
!elev  0  99  Quota stazione sul livello del mare
!whi   0  99  Indice di urbanita [0-1]
!rete  0  99  ID rete: 1 - RRQA, 2 - CMG, 4 - INM
!temp  0  99  ID sensore T 
!ur    0  99  ID sensore UR 
!prec  0  99  ID sensore PP
!ENDVARS 
!------------------------------------------------------------


! tutte le stazioni in anagrafica
!-------------------------------------------------------------------------------
 use coords
 use ll_utm

 implicit none

 integer, parameter :: KTOT=1000

! UNDEF here is the same as in the main code, both different from UNDEF in obs files!
 real(4), parameter :: UNDEF= -9999.0

 character(200) :: STAGFI
 character(8), dimension(KTOT) ::  stid
 real(8), dimension(KTOT) :: sx, sy, slon, slat, shgt
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
 type(UTM_COORD) :: tUTM

!-----------------------------------------------------------------------
! character(1) :: ans

!===============================================================================
 ! init -> UNDEF
 stid= ' '
 sx  = UNDEF
 sy  = UNDEF
 slon= UNDEF
 slat= UNDEF
 shgt= UNDEF

 !-----------------------------------------------------------------------
 write (6,*) 'anagread: file: ',trim(STAGFI)
 !write (6,*) 'ok?'
 !read (5,'(a1)') ans
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

 !-----------------------------------------------------------------------
 ! check for record terminator:
 do while (nlev /= 0)

  ! read variable codes:
  do j=1,NVARS
   irec= irec +1
   read (8,rec=irec) puff(j)
  end do

  ! assign variables:
  ! STation ID:
  stid(k)= aid
!  lstr= lstrim(aid)
!  read (aid(1:lstr),*) iid
!  write (stid(k)(3:7),'(i5.5)') iid
 ! write (6,'(a1,a8,a2,i4)') '[',stid(k),']=',iid
 ! if (mod(k,20) == 0) then
 !  write (6,*) 'ok?'
 !  read (5,'(a1)') ans
 ! end if

  sy(k)= dble(r4lay)
  sx(k)= dble(r4lox)

  shgt(k) = dble(puff(1))
  ! suhi(k) = dble(puff(29))

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
 write (6,*) 'anagread: read KKA=',KKA
 !write (6,*) 'ok?'
 !read (5,'(a1)') ans

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
! do k=1, KKA
!  write (6,'(a,i5,1x,a8,1x,2f8.2,f8.0)') &
!&  ' k id lon lat hgt var',k,stid(k),slon(k),slat(k),shgt(k)
!  if (mod(k,20) == 0) then
!   write (6,*) 'ok?'
!   read (5,'(a1)') ans
!  end if
! end do
! write (6,*) 'anagread: return'
! write (6,*) ' ok?'
! read (5,'(a1)') ans

!-------------------------------------------------------------------------------
 return
end subroutine anagreadUTM
!***************************************************************************************************
subroutine stzwri(FULFIST, nrec, KK, stid, s4x, s4y, shgt, so, sb, sa, sav)
!-------------------------------------------------------------------------------
! write obs, analysis etc. on GrADS "station" file
!-------------------------------------------------------------------------------
 implicit none
 integer, parameter :: KTOT=1000
 real(4), parameter :: UNDEF=-9999.0

 character(200) :: FULFIST
 ! integer, save :: nrec
 integer :: nrec
 integer :: KK
 character(8), dimension (KTOT) :: stid
 real(4), dimension (KTOT) :: s4x, s4y, shgt, so, sb, sa, sav

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
  ! background:
  nrec= nrec +1
  write (8,rec=nrec) sb(k)
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
 write (6,*) 'stzwri: KK and last nrec=', KK, nrec

 !write (6,*) 'stzwri: return '
 !write (6,*) 'ok?'
 !read (5,'(a1)') ans

!-------------------------------------------------------------------------------
 return
end subroutine stzwri
!*******************************************************************************
subroutine griwri(n, FULFIGR, xb, xa, xidi)
!-------------------------------------------------------------------------------
 implicit none
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
end subroutine griwri
!*******************************************************************************

! START udmat subroutines 
!  udsrnv
!  isortin
!  matsortin
!  invmpj
!  invmmj
!*******************************************************************************
subroutine udsrnv(KTOT, KKP, jarkp, SRNVP, KK, jark, ZS, SRNV)
!-------------------------------------------------------------------------------
! update matrix (S+R)^-1
!-------------------------------------------------------------------------------
! input :
!  KTOT
!  jarkp(KTOT), KKP, SRNVP(KTOT,KTOT)
!  jark(KTOT), KK
!  ZS(KTOT,KTOT)
! output:
!  SRNV(KTOT,KTOT)
! subroutines called:
!  isortin
!  invmmj
!  invmpj
!  matsortin
! Constraints: the code does not check, it must be:
!               KK <= KTOT
!              KKP <= KTOT
!
! Francesco Uboldi, 2006-2011
!-----------------------------------------------------------------------
! steps:
! - sort jark in growing order
! - discard absent stations
! - add new stations at the end
! - sort rows and columns
!-------------------------------------------------------------------------------
 implicit none

 real(8), parameter :: UNDEF= -9999.d0
!-----------------------------------------------------------------------
! input:

! array dimension:
 integer, intent(in) :: KTOT
! previous number of components:
 integer, intent(inout) :: KKP
! previous index arrays:
 integer, intent(inout), dimension(KTOT) :: jarkp
! previous inverse matrix:
 real(8), intent(inout), dimension(KTOT,KTOT) :: SRNVP
! full matrix:
 real(8), intent(in), dimension(KTOT,KTOT) :: ZS
! present number of active elements:
 integer, intent(in) :: KK
! present index array:
 integer, intent(inout), dimension(KTOT) :: jark

!-----------------------------------------------------------------------
! output:

! present inverse matrix:
 real(8), intent(out), dimension(KTOT,KTOT) :: SRNV

!-----------------------------------------------------------------------
! local variables

 integer, allocatable :: mex(:), mad(:)
 integer :: k, k1, k2, l, m, lab, lad, la

 integer inv
 real(8), allocatable :: fid(:,:) 
 real(8) :: fidres, fidmax, afid
! character(1) ans

!===============================================================================
! sort jark in growing order

 write (6,*) 'impongo che jark sia crescente '
 call isortin(jark,KK)
 write (6,*) 'dopo isortin, KK=',KK

!deb
! do k=1,KK
!  write (6,'(a,i4,i4)') 'k jark :', k, jark(k)
!!  write (6,'(a,i4,i4,a2,a8,a1,f6.1)') 'k jark ATID yo:' &
!!&     , k, jark(k), ' [',atid(jark(k)),']', ao(jark(k))
! end do
! write (6,*) ' ok?'
! read (5,'(a1)') ans

!-----------------------------------------------------------------------
! assign SRNVP to SRNV:

 write (6,*) 'assegno SRNVP a SRNV: KKP=', KKP
 SRNV= SRNVP
 write (6,*) ' ho assegnato SRNVP a SRNV ok?'
! read (5,'(a1)') ans

!-----------------------------------------------------------------------
! find stations to discard:

 write (6,*) 'trovo stazioni da togliere'
 allocate(mex(KTOT))
 mex= 0
 lab= 0
 do l=1,KKP
  m=1
  do while (jarkp(l) /= jark(m) .and. m <= KK)
   m= m +1
  end do
  if (m > KK) then
   lab= lab +1
   mex(lab)= jarkp(l)
   write (6,*) 'stazione da togliere:', mex(lab)
  end if
 end do
 !write (6,*) 'lab=',lab,' mex=',(mex(la),la=1,lab)
 call flush (6)
 !write (6,*) ' ok?'
 !read (5,'(a1)') ans

!-----------------------------------------------------------------------
! discard rows and columns corresponding to non-active stations:
 write (6,*) 'tolgo stazioni inattive'
 if (lab < KKP) then
  do la=1,lab
   write (6,*) ' tolgo dalla matrice la:', mex(la)
   call invmmj(ZS, KTOT, KKP, jarkp, SRNVP, mex(la), SRNV, inv)
   ! assegno SRNV a SRNVP:
   SRNVP= SRNV
   write (6,*) 'KKP KK = ', KKP, KK
  end do

 else ! lab >= KKP : le tolgo tutte

  if (lab > KKP) then
   write (6,*) ' !! lab < KKP !! : ', lab, KKP
   write (6,*) 'mi vuoi far togliere piu'' stazioni di quelle che ho ? azzero..'
  end if

  if (lab == KKP) write (6,*) ' ! lab = KKP ! : ', lab, KKP
  jarkp= 0
  SRNVP= UNDEF
  SRNV= UNDEF
  KKP= 0
 end if
 deallocate(mex)
 write (6,*) 'ho finito di togliere'
 call flush (6)
 !write (6,*) ' ok?'
 !read (5,'(a1)') ans

!-----------------------------------------------------------------------
! find stations to add

 write (6,*) 'trovo nuove stazioni da aggiungere'
 allocate(mad(KTOT))
 mad=0
 lad= 0
 do k=1,KK
  k2=1
  do while (jarkp(k2) /= jark(k) .and. k2 <= KKP)
   k2= k2 +1
  end do
  if (k2 > KKP) then
   lad= lad +1
   mad(lad)= jark(k)
   write (6,*) 'stazione da aggiungere:', mad(lad)
  end if
 end do
 !write (6,*) 'lad=',lad,' mad=',(mad(la),la=1,lad)
 !write (6,*) ' ok?'
 !read (5,'(a1)') ans

!-----------------------------------------------------------------------
! add new (rows/columns) stations at the end

 write (6,*) 'aggiungo nuove stazioni alla fine'
 do la=1,lad
  write (6,*) ' aggiungo alla fine la:', mad(la)
  KKP= KKP +1
  jarkp(KKP)= mad(la)
  write (6,*) 'KKP KK = ', KKP, KK
  call invmpj(ZS, KTOT, KKP, jarkp, SRNVP, SRNV, inv)
  ! assegno SRNV a SRNVP
  SRNVP= SRNV
 end do
 deallocate(mad)
 write (6,*) 'ho finito di aggiungere'
 call flush (6)
 !write (6,*) ' ok?'
 !read (5,'(a1)') ans

!-----------------------------------------------------------------------
! riordino righe e colonne della matrice secondo JARK

 write (6,*) 'riordino righe e colonne della matrice'
 call matsortin(SRNV, KTOT, jarkp, KK)

 ! assegno SRNV a SRNVP
 SRNVP= SRNV

 !write (6,*) 'KK= ', KKP, KK
 !do k=1,KK
 ! write (6,'(2(a7,i5))') ' jarkp:', jarkp(k), '  jark=',jark(k)
 !end do
! write (6,*) 'jarkp:', (jarkp(k),k=1,KK)
! write (6,*) 'jARK :', (jark(k),k=1,KK)
 write (6,*) 'KK= ', KKP, KK
 ! write (6,*) ' jarkp e jark sono uguali'
 call flush (6)
 ! write (6,*) ' ok?'
 ! read (5,'(a1)') ans

!-----------------------------------------------------------------------
! controllo matrice inversa

 write (6,*) 'provo l''inversa ZS(jark) SRNV', KK
 allocate(fid(KK,KK))
 do k=1,KK
  do k1=1,KK
   fid(k,k1)= 0.d0
   do k2=1,KK
    fid(k,k1)= fid(k,k1) +ZS(jark(k),jark(k2)) *SRNV(k2,k1)
   end do
  end do
  fid(k,k)= fid(k,k) -1.d0
 end do

 fidres= 0.d0
 fidmax= 0.d0
 do k=1,KK
  do k2=1,KK
   fidres= fidres +fid(k,k2)**2
   afid= abs(fid(k,k2))
   if (afid > fidmax) fidmax= afid
  end do
 end do
 deallocate(fid)
 fidres= sqrt(fidres /(KK *KK))
 write (6,*) 'identity res,max = ',fidres, fidmax
 call flush (6)
 !write (6,*) 'ok?'
 !read (5,'(a1)') ans

!-------------------------------------------------------------------------------
 return
end subroutine udsrnv
!*******************************************************************************
subroutine isortin(iar,NN)
!-------------------------------------------------------------------------------
! sort IAR so that IAR(n) < IAR(n+1)
! I integer array SORT in INcreasing order
!-------------------------------------------------------------------------------
 implicit none
 integer NN
 integer, dimension(NN) :: iar
 integer nsca, n, iaux
!===============================================================================
 nsca= 1
 do while (nsca > 0)
  nsca= 0
  do n=1,NN-1
   if (iar(n) > iar(n+1)) then
    iaux= iar(n)
    iar(n)= iar(n+1)
    iar(n+1)= iaux
    nsca= nsca+ 1
   end if
  end do
 end do
!-------------------------------------------------------------------------------
 return
end subroutine isortin
!*******************************************************************************
subroutine matsortin(xx,LDX,iar,NN)
!-------------------------------------------------------------------------------
! I integer array SORT in INcreasing order
! sort IAR so that IAR(n) < IAR(n+1)
! XX(j,jp)= ANOTHERMATRIX(IAR(j),IAR(jp)) : sort XX too
! Francesco Uboldi, 2006-2011
!-------------------------------------------------------------------------------
 implicit none
 integer LDX, NN
 real(8), dimension(LDX,LDX) :: xx
 integer, dimension(NN) :: iar
 real(8) aux
 integer nsca, n, iaux, j
!===============================================================================
 nsca= 1
 do while (nsca > 0)
  nsca= 0
  do n=1,NN-1
   if (iar(n) > iar(n+1)) then
    iaux= iar(n)
    iar(n)= iar(n+1)
    iar(n+1)= iaux
    do j=1,NN 
     aux= xx(j,n)
     xx(j,n)= xx(j,n+1)
     xx(j,n+1)= aux
    end do
    do j=1,NN 
     aux= xx(n,j)
     xx(n,j)= xx(n+1,j)
     xx(n+1,j)= aux
    end do
    nsca= nsca+ 1
   end if
  end do
 end do
!-------------------------------------------------------------------------------
 return
end subroutine matsortin
!*******************************************************************************
! "sort" subroutines end
! inv_j subroutines start
!*******************************************************************************
subroutine invmpj(aa, LDA, KK, jark, vnm1, vn, inv)
!-------------------------------------------------------------------------------
! all matrices are symmetric
!
! compute inverse of a order KK matrix 
! given the inverse of a submatrix of order KK-1
!
! HERE "J" means the following:
! the matrix of order KK-1 -whose inverse, VNM1, is known- is
! composed of KK-1 row-columns of the (bigger) matrix AA
! these row-columns are: jark(1)..jark(KK-1)
! jark(KK) is the row-column to add, so that the matrix of order KK
! is given by KK row-columns of AA: jark(1)...jark(KK)
! The inverse of this last matrix is calculated
! LDA is the leading dimension of matrices AA, VNM1, VN  (KK <= LDA)
!
! the use of jark results in changes only in the assignment of 
! the auxiliary variables BB(.) and ALP

!
! AA       input : symmetric matrix of order LDA
! LDA      input : its dimension
! KK       input : order of submatrix
! JARK     input : row-columns of AA composing the submatrix of order KK-1
! VMN1     input : inverse of the submatrix of order KK-1
! JARK(KK) input : row-column of AA that composes the submatrix of order KK
!            as its KK row-column
! VN      output : inverse of the submatrix of order KK
! INV     output : flag 0=OK
!
! Francesco Uboldi, 2006-2011
!-------------------------------------------------------------------------------
 implicit none
 integer LDA, KK
 real(8), dimension(LDA,LDA):: aa, vn,  vnm1
 integer, dimension(LDA) :: jark
 integer inv

 real(8), allocatable :: bb(:), dd(:), csi(:)
 real(8) alp, gam
 integer m, n

!===============================================================================
 allocate(bb(KK), dd(KK), csi(KK))

!-----------------------------------------------------------------------
! assign BB, alp
 do n=1,KK-1
  bb(n)= aa(jark(n),jark(KK))
 end do
 alp= aa(jark(KK),jark(KK))

!-----------------------------------------------------------------------
! assign CSI
 do n=1,KK-1
  csi(n)= 0.d0
  do m=1,KK-1
   csi(n)= csi(n) +vnm1(n,m) *bb(m)
  end do
 end do

!-----------------------------------------------------------------------
! assign GAM
 gam= 0.d0
 do n=1,KK-1
  gam= gam +bb(n) *csi(n)
 end do

!-----------------------------------------------------------------------
 if (abs(alp-gam) < 1.d-14) then
  write (6,*) 'invmpj: alp gam diff', alp,gam, (alp-gam)
  write (6,*) 'invmpj: singular matrix', KK
  inv= 1
  deallocate(bb, dd, csi)
  return
 end if
 gam= 1.d0 /(alp -gam)

!-----------------------------------------------------------------------
! assign DD
 do n=1,KK-1
  dd(n)= -gam*csi(n)
 end do

!-----------------------------------------------------------------------
! assign VN first N-1 row, columns
 do n=1,KK-1
  do m=n,KK-1
   vn(n,m)= vnm1(n,m) +gam *csi(n) *csi(m)
   vn(m,n)= vn(n,m)
  end do
 end do

!-----------------------------------------------------------------------
! complete VN
 do n=1,KK-1
  vn(n,KK)= dd(n)
  vn(KK,n)= dd(n)
 end do

 vn(KK,KK)= gam
 deallocate(bb, dd, csi)

!-------------------------------------------------------------------------------
 return
end subroutine invmpj
!*******************************************************************************
subroutine invmmj(aa, LDM, MM, jark, bbm, mex, aanv, inv)
!-------------------------------------------------------------------------------
! all matrices are symmetric
!
! compute inverse of an order MM-1 submatrix of an order MM matrix
! given the inverse of the order MM matrix
!
! HERE "J" means the following:
! the matrix of order MM -whose inverse, BBM, is known- is
! composed of MM row-columns of the (bigger) matrix AA
! these row-columns are: jark(1)..jark(MM)
! MEX is the row-column to exclude: it is excluded from jark, too,
! so that, on output, the matrix of order MM-1
! is given by MM-1 row-columns of AA: jark(1)...jark(MM-1)
! The inverse of this last matrix is calculated
! LDM is the leading dimension of matrices AA, BBM, AANV  (MM <= LDM)
!
! AA       input : symmetric matrix of order LDM
! LDM      input : its dimension
! MM       input : order of submatrix
! JARK     input : row-columns of AA composing the submatrix of order MM
! BBM      input : inverse of the submatrix of order MM
! MEX      input : row-column of AA to be excluded
! AANV    output : inverse of the submatrix of order MM-1
! INV     output : flag 0=OK
!
! Francesco Uboldi, 2006-2011
!-------------------------------------------------------------------------------
! INV=  0 : inversion OK
! INV= -1 : inversion FAILED
! INV=  1 : inversion done, but better to check (small denominator)
! INV= 10 : MEX was not found in JARK: AANV=BBM no row-column was excluded
!-------------------------------------------------------------------------------
 implicit none
 integer LDM, MM
 real(8), dimension(LDM,LDM) :: aa, bbm, aanv
 integer, dimension(LDM) :: jark
 integer mex, inv
 integer MM1, m, l, k, lm, km, jmex
!===============================================================================
 inv= 0
!-----------------------------------------------------------------------
 if (MM <= 1) then
  write (6,*) 'invmm1: grunt! MM=', MM
  inv= -1
  return
 end if

!-----------------------------------------------------------------------
 m=1
 do while (mex /= jark(m) .and. m <= MM)
  m= m +1
 end do
 if (m > MM) then
  write (6,*) 'invmm1: mex is different from all jark elements:'
  aanv= bbm
  write (6,*) 'invmm1: no row-columns were excluded'
  inv= 10
  return
 end if

!-----------------------------------------------------------------------
! else:
 jmex= m

!-----------------------------------------------------------------------
 if (bbm(jmex,jmex) == 0.d0) then
  write (6,*) 'invmm1: I cannot invert this! bbm(jmex,jmex)=' &
& , bbm(jmex,jmex)
  inv= -1
  return
 end if

!-----------------------------------------------------------------------
 if (abs(bbm(jmex,jmex)) <= 1.d-13) then
  write (6,*) 'invmm1: I''m not sure I can do this..' &
& , ' bbm(jmex,jmex)=', bbm(jmex,jmex)
  write (6,*) 'however..'
  inv= 1
 end if

!-----------------------------------------------------------------------
 MM1= MM -1

!-----------------------------------------------------------------------
 ! write (6,*) 'invmmj: mex jmex',mex, jmex
 do l=1,MM1
  lm= l
  if (l >= jmex) lm= l+1
  do k=1,MM1
   km= k
   if (k >= jmex) km= k+1
   aanv(l,k)= bbm(lm,km)
  end do
  aanv(MM,l)= bbm(jmex,lm)
  ! aanv(l,MM)= bbm(jmex,lm)
  aanv(l,MM)= bbm(lm,jmex)
 end do
 aanv(MM,MM)= bbm(jmex,jmex)

!-----------------------------------------------------------------------
 do l=1,MM1
  do k=1,MM1
   aanv(l,k)= aanv(l,k)- aanv(MM,l)*aanv(MM,k)/aanv(MM,MM)
  end do
 end do

! exclude mex from jark:
 do l= jmex, MM1
  jark(l)= jark(l+1)
 end do

! memory:
 jark(MM)= mex

 MM= MM1

!-------------------------------------------------------------------------------
 return
end subroutine invmmj
!*******************************************************************************
! udmat subroutines end

!*******************************************************************************
integer function lstrim(buff)
!-------------------------------------------------------------------------------
! eliminate final blanks and final not printable characters
! - buff (ch): string (input)
! lstrim = lenght without final blanks
!-------------------------------------------------------------------------------
 character buff*(*)
!===============================================================================
 do k=len(buff),1,-1
  if (buff(k:k) > ' ' .and. buff(k:k) < '~') then
   lstrim= k
   return
  end if
 end do
 lstrim= 0
!-------------------------------------------------------------------------------
 return
end function lstrim
!*******************************************************************************
real(8) function distan(theta1, phi1, theta2, phi2)
                      ! lon1    lat1  lon2    lat2
!-------------------------------------------------------------------------------
! a(rad)= a(deg)* rtd
!-------------------------------------------------------------------------------
 implicit none
 ! real(8), parameter :: UNDEF= -9999.d0, RT=6371.d0
 real(8), parameter :: UNDEF= -9999.d0, RT=6378.d0

 real(8) :: phi1, phi2, theta1, theta2

 real(8) :: hpi, rtd, gam, alpha, rphi1, rphi2, rtheta1, rtheta2
!===============================================================================

 if (phi1 == UNDEF .or.phi2 == UNDEF .or. theta1 == UNDEF .or. theta2 == UNDEF) then
  distan= UNDEF
  return
 end if

 if (phi1 == phi2 .and. theta1 == theta2) then
  distan= 0.d0
  return
 end if

!-------------------------------------------------------------------------------
! gam= cosd(phi1) *cosd(phi2) *cosd(theta1 -theta2) +sind(phi1) *sind(phi2)
!-------------------------------------------------------------------------------

 hpi= asin(1.d0)
 rtd= hpi /90.d0
 rphi1= rtd *phi1
 rphi2= rtd *phi2
 rtheta1= rtd *theta1
 rtheta2= rtd *theta2

 gam= cos(rphi1) *cos(rphi2) *cos(rtheta1 -rtheta2) + sin(rphi1) *sin(rphi2)
 alpha= acos(gam)

 distan= RT* alpha

!-------------------------------------------------------------------------------
 return
end function distan
!*******************************************************************************
