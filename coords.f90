module COORDS

  type Ellipsoid
    integer    :: id
    character(len=32)  :: ellipsoidName
    real(kind=8)    :: A
    real(kind=8)    :: E2
  end type Ellipsoid

    integer, parameter  :: AIRY                     =  1
    integer, parameter  :: AUSTRALIAN_NATIONAL      =  2
    integer, parameter  :: BESSEL_1841              =  3
    integer, parameter  :: BESSEL_1841_NAMBIA       =  4
    integer, parameter  :: CLARKE_1866              =  5
    integer, parameter  :: CLARKE_1880              =  6
    integer, parameter  :: EVEREST                  =  7
    integer, parameter  :: FISCHER_1960_MERCURY     =  8
    integer, parameter  :: FISCHER_1968             =  9
    integer, parameter  :: GRS_1967                 = 10
    integer, parameter  :: GRS_1980                 = 11
    integer, parameter  :: HELMERT_1906             = 12
    integer, parameter  :: HOUGH                    = 13
    integer, parameter  :: INTERNATIONAL            = 14
    integer, parameter  :: KRASSOVSKY               = 15
    integer, parameter  :: MODIFIED_AIRY            = 16
    integer, parameter  :: MODIFIED_EVEREST         = 17
  integer, parameter  :: MODIFIED_FISCHER_1960  = 18
    integer, parameter  :: SOUTH_AMERICAN_1969      = 19
    integer, parameter  :: WGS_60                   = 20
    integer, parameter  :: WGS_66                   = 21
    integer, parameter  :: WGS_72                   = 22
    integer, parameter  :: WGS_84                   = 23

  type (Ellipsoid), dimension(0:23), parameter  ::  DATUMS = (/      &
    Ellipsoid( -1, "Placeholder", 0._8, 0._8),              &
    Ellipsoid( 1, "Airy", 6377563._8, 0.00667054_8),          &
    Ellipsoid( 2, "Australian National", 6378160._8, 0.006694542_8),  &
    Ellipsoid( 3, "Bessel 1841", 6377397._8, 0.006674372_8),      &
    Ellipsoid( 4, "Bessel 1841 (Nambia) ", 6377484._8, 0.006674372_8),  &
    Ellipsoid( 5, "Clarke 1866", 6378206._8, 0.006768658_8),      &
    Ellipsoid( 6, "Clarke 1880", 6378249._8, 0.006803511_8),      &
    Ellipsoid( 7, "Everest", 6377276._8, 0.006637847_8),        &
    Ellipsoid( 8, "Fischer 1960 (Mercury) ", 6378166._8, 0.006693422_8),&
    Ellipsoid( 9, "Fischer 1968", 6378150._8, 0.006693422_8),      &
    Ellipsoid( 10, "GRS 1967", 6378160._8, 0.006694605_8),        &
    Ellipsoid( 11, "GRS 1980", 6378137._8, 0.00669438_8),        &
    Ellipsoid( 12, "Helmert 1906", 6378200._8, 0.006693422_8),      &
    Ellipsoid( 13, "Hough", 6378270._8, 0.00672267_8),          &
    Ellipsoid( 14, "International", 6378388._8, 0.00672267_8),      &
    Ellipsoid( 15, "Krassovsky", 6378245._8, 0.006693422_8),      &
    Ellipsoid( 16, "Modified Airy", 6377340._8, 0.00667054_8),      &
    Ellipsoid( 17, "Modified Everest", 6377304._8, 0.006637847_8),    &
    Ellipsoid( 18, "Modified Fischer 1960", 6378155._8, 0.006693422_8),  &
    Ellipsoid( 19, "South American 1969", 6378160._8, 0.006694542_8),  &
    Ellipsoid( 20, "WGS 60", 6378165._8, 0.006693422_8),        &
    Ellipsoid( 21, "WGS 66", 6378145._8, 0.006694542_8),        &
    Ellipsoid( 22, "WGS-72", 6378135._8, 0.006694318_8),        &
    Ellipsoid( 23, "WGS-84", 6378137._8, 0.00669438_8)          &
  /)

  ! Latitude/Longitude coordinates
  type GEO_COORD
    real(8)  :: Lat
    real(8)  :: Lon
  end type GEO_COORD

    ! Lokalmodell rotated geographic coordinates
  type LAMI_COORD
    real(8)  :: Lat
    real(8)  :: Lon
  end type LAMI_COORD

  ! Geocentric metric coordinates
  type GEOCENTRIC
    real(8)  :: X
    real(8)  :: Y
    real(8)  :: Z
  end type GEOCENTRIC

  ! UTM coordinates
  type UTM_COORD
    character(len=4)  :: Zone
    real(8)    :: N
    real(8)    :: E
  end type UTM_COORD

  ! Gauss-Boaga coordinates
  integer, parameter  :: WEST = 1
  integer, parameter  :: EAST = 2
  type GAUSS_COORD
    integer  :: S
    real(8)  :: N
    real(8)  :: E
  end type GAUSS_COORD

end module COORDS
