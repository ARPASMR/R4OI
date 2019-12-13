module LatLon

  use COORDS  ! Import ellipsoid data

contains

  function CONVERT1(G, IN_ELL) result(XYZ)
    implicit none
    ! Call list
    type (GEO_COORD), intent(in)  :: G
    integer, intent(in)      :: IN_ELL
    type (GEOCENTRIC)      :: XYZ
    ! Locals
    real(8)  :: PHI
    real(8)  :: LAMBDA
    real(8)  :: SPH
    real(8)  :: CPH
    real(8)  :: A
    real(8)  :: E2
    real(8)  :: FCT
    real(8)  :: Pi
    ! Auxiliary
    Pi = DATAN(1.d0)*4.d0
    PHI = G%LAT * Pi/180.d0
    LAMBDA = G%LON * Pi/180.d0
    SPH = dsin(PHI)
    CPH = dcos(PHI)
    A  = DATUMS(IN_ELL)%A
    E2 = DATUMS(IN_ELL)%E2
    FCT = A/dsqrt(1.d0-E2*SPH**2)
    ! Actual conversion
    XYZ%X = FCT*CPH*dcos(LAMBDA)
    XYZ%Y = FCT*CPH*dsin(LAMBDA)
    XYZ%Z = FCT*(1.d0-E2)*SPH
  end function CONVERT1

  function CONVERT2(G, OUT_ELL) result(Q)
    ! Call list
    type (GEOCENTRIC), intent(in)  :: G
    integer, intent(in)      :: OUT_ELL
    type (GEO_COORD)      :: Q
    ! Locals
    real(8)  :: A
    real(8)  :: E2
    real(8)  :: EP
    real(8)  :: B
    real(8)  :: P
    real(8)  :: PSI
    real(8)  :: PHI
    real(8)  :: LAMBDA
    real(8)  :: Pi
    ! Auxiliary
    Pi  = DATAN(1.d0)*4.d0
    A   = DATUMS(OUT_ELL)%A
    E2  = DATUMS(OUT_ELL)%E2
    EP  = E2/(1.d0-E2)
    B   = A*dsqrt(1.d0-E2)
    P   = dsqrt(G%X**2+G%Y**2)
    PSI = datan2(G%Z*A,P*B)
    LAMBDA = datan2(G%Y,G%X)
    PHI = datan((G%Z+EP*B*dsin(PSI)**3)/(P-E2*A*dcos(PSI)**3))
    Q%LAT = PHI*180.d0/Pi
    Q%LON = LAMBDA*180.d0/Pi
  end function CONVERT2

  function CONVERT(G, IN_ELL, OUT_ELL) result(H)
    ! Call list
    type (GEO_COORD), intent(in)  :: G
    integer, intent(in)      :: IN_ELL
    integer, intent(in)      :: OUT_ELL
    type (GEO_COORD)      :: H
    ! Locals
    type (GEOCENTRIC)      :: XYZ
    ! Conversions
    XYZ = CONVERT1(G,IN_ELL)
    H = CONVERT2(XYZ,OUT_ELL)
  end function CONVERT

end module LatLon

