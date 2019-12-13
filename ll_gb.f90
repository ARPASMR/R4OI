module GAUSS_BOAGA

  use COORDS

  private CONT, EVAL, N, RHO_PHI

  real(8), parameter  :: CONT = 0.9996d0  ! Contraction coefficient

contains

  ! Ancillary (private) routines

  function EVAL(P, imin, imax, X) result(Y)
    implicit none
    real(8), dimension(*), intent(in)  :: P
    real(8), intent(in)      :: X
    integer, intent(in)      :: imin, imax
    real(8)            :: Y
    integer            :: i1,i
    i1 = imax - imin + 1
    Y = P(i1)
    do i=i1,1,-1
      Y = Y*X + P(i)
    end do
    if(imin>0) then
      Y = Y*X**imin
    end if
  end function EVAL

  function N(PHI,A,E) result(RetVal)
    real(8), intent(in)  :: PHI
    real(8), intent(in)  :: A
    real(8), intent(in)  :: E
    real(8)    :: RetVal
    RetVal = A/dsqrt(1.0d0-(E*dsin(PHI))**2)
  end function N

  function RHO_INT(PHI,A,E2,AE,BE,CE,DE,EE,FE) result(RetVal)
    real(8), intent(in)  :: PHI
    real(8), intent(in)  :: A
    real(8), intent(in)  :: E2
    real(8), intent(in)  :: AE
    real(8), intent(in)  :: BE
    real(8), intent(in)  :: CE
    real(8), intent(in)  :: DE
    real(8), intent(in)  :: EE
    real(8), intent(in)  :: FE
    real(8)    :: RetVal
    RetVal = A*(1.0d0-E2)* &
      ( &
        AE*PHI &
        -BE*dsin( 2.0d0*PHI) &
        +CE*dsin( 4.0d0*PHI) &
        -DE*dsin( 6.0d0*PHI) &
        +EE*dsin( 8.0d0*PHI) &
        -FE*dsin(10.0d0*PHI) &
      )
  end function RHO_INT

  ! Actual conversions

  function LL_TO_GB(L) result(RSL)

    implicit none

    ! Call list
    type (GEO_COORD), intent(in)  :: L
    type (GAUSS_COORD)      :: RSL

    ! Locals
    real(8)  :: A, E, E2
    real(8)  :: Pi
    real(8)  :: PHI
    real(8)  :: LAMBDA
    real(8)  :: ETA2
    real(8)  :: NG
    real(8)  :: EG
    real(8)  :: LAMBDA2
    real(8)  :: ORIG
    real(8)  :: E0
    real(8)  :: AE
    real(8)  :: BE
    real(8)  :: CE
    real(8)  :: DE
    real(8)  :: EE
    real(8)  :: FE
    real(8)  :: T1
    real(8)  :: U0
    real(8)  :: U1
    real(8)  :: Q1
    real(8)  :: T2
    real(8)  :: V0
    real(8)  :: V1
    real(8)  :: V2
    real(8)  :: Q2

    ! Constants
    real(8), dimension(0:5), parameter    :: AV = &
& (/ 1.0d0, 3.0d0/4.0d0, 45.0d0/64.0d0, 175.0d0/256.0d0, 11025.0d0/16384.0d0, 43659.0d0/65536.0d0 /)
    real(8), dimension(1:5), parameter    :: BV = &
& (/ 3.0d0/4.0d0, 15.0d0/16.0d0, 525.0d0/512.0d0, 2205.0d0/2048.0d0, 72765.0d0/65536.0d0 /)
    real(8), dimension(2:5), parameter    :: CV = &
& (/ 15.0d0/64.0d0, 105.0d0/256.0d0, 2205.0d0/4096.0d0, 10395.0d0/16384.0d0 /)
    real(8), dimension(3:5), parameter    :: DV = &
& (/ 35.0d0/512.0d0, 315.0d0/2048.0d0, 31185.0d0/131072.0d0 /)
    real(8), dimension(4:5), parameter    :: EV = &
& (/ 315.0d0/16384.0d0, 3465.0d0/65536.0d0 /)
    real(8), dimension(5:5), parameter    :: FV = &
& (/ 693.0d0/131072.0d0 /)

    ! Fix ellipsoid parameters
    A = DATUMS(INTERNATIONAL)%A
    E2 = DATUMS(INTERNATIONAL)%E2
    E = dsqrt(E2)
    write (6,*) 'LL_to_GB: a e2 e', A, E2, E

    ! Polynomial evaluations
    AE = EVAL(AV,0,5,E2)
    BE = EVAL(BV,1,5,E2)/2.0d0
    CE = EVAL(CV,2,5,E2)/4.0d0
    DE = EVAL(DV,3,5,E2)/6.0d0
    EE = EVAL(EV,4,5,E2)/8.0d0
    FE = EVAL(FV,5,5,E2)/10.0d0

    ! Fuse determination
    if(L%LON<=12.d0 .and. L%LON>=6.0) then
      RSL%S = WEST
      ORIG = 9.0d0
      E0 = 1500000.0d0
    elseif(L%LON<=18.0d0 .and. L%LON>=12.0) then
      RSL%S = EAST
      ORIG = 15.0d0
      E0 = 2520000.0d0
    else
      RSL%S = WEST
      ORIG = 9.0d0
      E0 = 1500000.0d0
    end if

    ! Preliminary values
    Pi = datan(1.d0)*4.d0
    PHI = L%LAT * Pi / 180.0d0
    LAMBDA = (L%LON - ORIG) * Pi / 180.0d0
    LAMBDA2 = LAMBDA**2
    ETA2 = 0.00676817d0*dcos(PHI)**2
    T1 = LAMBDA2*dsin(PHI)*dcos(PHI)/2.0d0
    U0 = 1.0d0
    U1 = LAMBDA2*(5.0d0-dtan(PHI)**2+9.0d0*ETA2+4.0d0*ETA2**2)/12.0d0
    Q1 = U0 + U1*dcos(PHI)**2
    T2 = LAMBDA*dcos(PHI)
    V0 = 1.0d0
    V1 = LAMBDA2*(1.0d0-dtan(PHI)**2+ETA2)/6.0d0
    V2 = (LAMBDA2**2) * (5.0d0-18.0d0*dtan(PHI)**2+dtan(PHI)**4)/120.0d0
    Q2 = V0 + V1*dcos(PHI)**2 + V2*dcos(PHI)**4
    NG = RHO_INT(PHI,A,E2,AE,BE,CE,DE,EE,FE) + N(PHI,A,E)*T1*Q1
    EG = N(PHI,A,E)*T2*Q2

    ! Return results
    RSL%N = CONT*NG
    RSL%E = CONT*EG + E0

  end function LL_TO_GB

  function GB_TO_LL(G) result(GG)

    implicit none

    ! Call list
    type (GAUSS_COORD), intent(in)  :: G
    type (GEO_COORD)      :: GG

    ! Locals
    real(8)  :: Pi
    real(8)  :: Pi4
    real(8)  :: A
    real(8)  :: E2
    real(8)  :: X
    real(8)  :: Y
    real(8)  :: NF
    real(8)  :: EF
    real(8)  :: PHI0
    real(8)  :: LAMBDA0
    real(8)  :: AA
    real(8)  :: BB
    real(8)  :: ALPHA
    real(8)  :: BETA
    real(8)  :: NU
    real(8)  :: PSI
    real(8)  :: PHI
    real(8)  :: LAMBDA
    real(8)  :: M0
    real(8)  :: M1
    real(8)  :: M2
    real(8)  :: M3
    real(8)  :: LM
    real(8), dimension(0:4), parameter  :: MP = (/ 1.0d0, 1.0d0/4.0d0, 7.0d0/64.0d0, 15.0d0/256.0d0, 579.0d0/16384.0d0 /)
    real(8), dimension(1:4), parameter  :: AL = (/ 1.0d0/4.0d0, 3.0d0/64.0d0, 5.0d0/256.0d0, 175.0d0/16384.0d0 /)
    real(8), dimension(2:4), parameter  :: BL = (/ 1.0d0/24.0d0, 35.0d0/768.0d0, 71.0d0/1536.0d0 /)
    real(8), dimension(3:4), parameter  :: CL = (/ 127.0d0/5808.0d0, 9803.0d0/185856.0d0 /)
    real(8), dimension(4:4), parameter  :: DL = (/ 33239.0d0/1626240.0d0 /)
    real(8), dimension(0:3), parameter  :: AF = (/ 1.0d0, 1.0d0, 1.0d0, 1.0d0 /)
    real(8), dimension(1:3), parameter  :: BF = (/ 7.0d0, 17.0d0, 30.0d0 /)
    real(8), dimension(2:3), parameter  :: CF = (/ 224.0d0, 889.0d0 /)
    real(8), dimension(3:3), parameter  :: DF = (/ 4279.0d0 /)
    real(8), dimension(0:3)    :: POL
    real(8)  :: L0
    real(8)  :: L1
    real(8)  :: L2
    real(8)  :: L3
    real(8)  :: M
    real(8)  :: P0
    real(8)  :: P1
    real(8)  :: P2
    real(8)  :: P3

    ! Initialisation
    Pi4 = datan(1.d0)
    Pi  = Pi4*4.0d0
    A = DATUMS(INTERNATIONAL)%A
    E2 = DATUMS(INTERNATIONAL)%E2
    L0 = EVAL(AL,1,4,E2)
    L1 = EVAL(BL,2,4,E2)
    L2 = EVAL(CL,3,4,E2)
    L3 = EVAL(DL,4,4,E2)
    M  = EVAL(MP,0,4,E2)/A
    P0 = EVAL(AF,0,3,E2)
    P1 = EVAL(BF,1,3,E2)/6.0d0
    P2 = EVAL(CF,2,3,E2)/120.0d0
    P3 = EVAL(DF,3,3,E2)/1260.0d0
    POL = (/ P0, -P1, P2, -P3 /)

    ! Deducing false origin based on fuse
    if(G%S==WEST) then
      NF = 0.0d0
      EF = 1500000.0d0
      PHI0 = 0.0d0
      LAMBDA0 = 9.0d0
    else
      NF = 0.0d0
      EF = 2520000.0d0
      PHI0 = 0.0d0
      LAMBDA0 = 15.0d0
    end if

    ! Computing X and Y as local fuse Gauss coordinates
    X = M*(G%N - NF)/0.9996d0
    Y = M*(G%E - EF)/0.9996d0

    ! Computes auxiliary values
    AA = dsin(X)*dcosh(Y)
    BB = dcos(X)*dsinh(Y)
    ALPHA = dasin(dsin(X)/dcosh(Y))
    BETA  = datan(dsinh(Y)/dcos(X))
    M0 = AA
    M1 = AA*(AA**2-3.0d0*BB**2)
    M2 = AA*(AA**4-10.0d0*(AA*BB)**2+5.0d0*BB**4)
    M3 = AA*(AA**6-21.0d0*AA**4*BB**2+35.0d0*AA**2*BB**4-7.0d0*BB**6)
    LM = -L0*M0+L1*M1-L2*M2+L3*M3
    NU = dlog(dtan(Pi4+(ALPHA/2.0d0))) + LM
    PSI = 2.0d0*(datan(dexp(1.0d0)**NU)-Pi4)

    ! Latitude in radian form
    PHI = PSI+(E2/2.0d0)*dsin(2.0d0*PSI)*EVAL(POL,0,3,dsin(PSI)**2)

    ! Longitude in radian form
    M0 =  BB
    M1 = -BB*(BB**2-3.0d0*AA**2)
    M2 =  BB*(BB**4-10.0d0*(BB*AA)**2+5.0d0*AA**4)
    M3 = -BB*(BB**6-21.0d0*BB**4*AA**2+35.0d0*BB**2*AA**4-7.0d0*AA**6)
    LM = -L0*M0 + L1*M1 - L2*M2 + L3*M3
    LAMBDA = BETA + LM

    ! Converts latitude and longitude to final form
    GG%LAT = PHI*180.0d0/Pi + PHI0
    GG%LON = LAMBDA*180.0d0/Pi + LAMBDA0

  end function GB_TO_LL

end module GAUSS_BOAGA
