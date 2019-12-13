module ll_utm

  ! Equations from USGS Bulletin 1532 
  ! East Longitudes are positive, West longitudes are negative.
  ! North latitudes are positive, South latitudes are negative
  ! Lat and Long are in decimal degrees
  ! Written by Chuck Gantz- chuck.gantz@globalstar.com (C++ version)
  ! Translated to Fortran 90 by M.Favaron (Servizi Territorio coop. a r.l.)

  use COORDS  ! Import ellipsoid and coordinate definitions

contains

  subroutine LLtoUTM(ReferenceEllipsoid, G, UTM)

    integer, intent(in)      :: ReferenceEllipsoid
    type (GEO_COORD), intent(in)  :: G
    type (UTM_COORD)      :: UTM

    ! Locals
    real(kind=8)    ::  UTMNorthing
    real(kind=8)    ::  UTMEasting
    character(len=4)  ::  UTMZone
    real(kind=8)    ::  Lat
    real(kind=8)    ::  Long
    real(kind=8)    ::  a, eccSquared, k0, LongOrigin, eccPrimeSquared
    real(kind=8)    ::  N, T, C, Aa, M
    real(kind=8)    ::  LongTemp, LatRad, LongRad, LongOriginRad
    integer    ::  ZoneNumber

    ! Parameters
    real(kind=8)    ::  PI
    real(kind=8)    ::  FOURTHPI
    real(kind=8)    ::  deg2rad
    real(kind=8)    ::  rad2deg

    ! Initial values
    FourthPi = datan(1.d0)
    Pi = FourthPi*4.d0
    deg2rad = Pi/180.d0
    rad2deg = 180.d0/Pi

    ! Initial transfer
    Lat = G%Lat
    Long = G%Lon

    ! Actual conversions
    a = DATUMS(ReferenceEllipsoid)%A
    eccSquared = DATUMS(ReferenceEllipsoid)%E2
    k0 = 0.9996_8
    LongTemp = Long-int((Long+180._8)/360._8)*360._8
    LatRad = Lat*deg2rad
    LongRad = LongTemp*deg2rad
    ZoneNumber = int((LongTemp + 180._8)/6._8) + 1
    if(Lat>=56._8 .and. Lat<64._8 .and. LongTemp>=3._8 .and. LongTemp<12._8) ZoneNumber = 32
    ! Special zones for Svalbard
    if(Lat>=72._8 .and. Lat<84._8) then
      if(LongTemp>=0._8 .and. LongTemp<9._8 ) then
        ZoneNumber = 31
      elseif(LongTemp>=9._8 .and. LongTemp<21._8) then
        ZoneNumber = 33
      elseif(LongTemp>=21._8 .and. LongTemp<33._8) then
        ZoneNumber = 35
      elseif(LongTemp>=33._8 .and. LongTemp<42._8) then
        ZoneNumber = 37
      end if
    end if
!
! blindato se c'e' il +: '+32T'
      if (UTM%Zone(1:1) == '+') then
       read (UTM%Zone,'(i3)') ZoneNumber
       write (6,*) 'LLtoUTM: UTM%Zone=[', UTM%Zone,']'
       write (6,*) 'LLtoUTM: locked ZoneNumber=', ZoneNumber
      end if

!
    LongOrigin = (ZoneNumber - 1)*6._8 - 180._8 + 3._8
    LongOriginRad = LongOrigin * deg2rad
    write(UTMZone,"(i3,a1)") ZoneNumber, UTMLetterDesignator(Lat)
    eccPrimeSquared = eccSquared/(1._8-eccSquared)
    ! Original: N = a/sqrt(1._8-eccSquared*sin(LatRad)*sin(LatRad))
    N = a/sqrt(1._8-eccSquared*sin(LatRad)**2)
    ! Original: T = tan(LatRad)*tan(LatRad)
    T = tan(LatRad)**2
    ! Original: C = eccPrimeSquared*cos(LatRad)*cos(LatRad)
    C = eccPrimeSquared*cos(LatRad)**2
    Aa = cos(LatRad)*(LongRad-LongOriginRad)
    M = a*((1._8 - eccSquared/4._8 - 3._8*eccSquared**2/64._8 - 5._8*eccSquared**3/256._8)*LatRad        &
      - (3._8*eccSquared/8._8 + 3._8*eccSquared**2/32._8 + 45._8*eccSquared**3/1024._8)*sin(2._8*LatRad)  &
      + (15._8*eccSquared**2/256._8 + 45._8*eccSquared**3/1024._8)*sin(4._8*LatRad)            &
      - (35._8*eccSquared**3/3072._8)*sin(6._8*LatRad))
    UTMEasting = (k0*N*(Aa+(1._8-T+C)*Aa**3/6._8              &
      + (5._8-18._8*T+T**2+72._8*C-58._8*eccPrimeSquared)*Aa**5/120._8)  &
      + 500000.0_8)
    UTMNorthing = (k0*(M+N*tan(LatRad)*(Aa**2/2._8+(5._8-T+9._8*C+4._8*C**2)*Aa**4/24._8  &
       + (61._8-58._8*T+T**2+600._8*C-330._8*eccPrimeSquared)*Aa**6/720._8)))
    if(Lat<0._8) UTMNorthing = UTMNorthing + 10000000.0_8

    ! Transfer
    UTM%Zone = UTMZone
    UTM%N    = UTMNorthing
    UTM%E    = UTMEasting

  end subroutine LLtoUTM



  function UTMLetterDesignator(Lat) result(LetterDesignator)
    real(kind=8), intent(in)  ::  Lat
    character(len=1)      ::  LetterDesignator

    if((84._8>=Lat) .and. (Lat>=72._8)) then
      LetterDesignator = 'X'
    elseif((72._8>Lat) .and. (Lat>=64._8)) then
      LetterDesignator = 'W'
    elseif((64._8>Lat) .and. (Lat>=56._8)) then
      LetterDesignator = 'V'
    elseif((56._8>Lat) .and. (Lat>=48._8)) then
      LetterDesignator = 'U'
    elseif((48._8>Lat) .and. (Lat>=40._8)) then
      LetterDesignator = 'T'
    elseif((40._8>Lat) .and. (Lat>=32._8)) then
      LetterDesignator = 'S'
    elseif((32._8>Lat) .and. (Lat>=24._8)) then
      LetterDesignator = 'R'
    elseif((24._8>Lat) .and. (Lat>=16._8)) then
      LetterDesignator = 'Q'
    elseif((16._8>Lat) .and. (Lat>=8._8)) then
      LetterDesignator = 'P'
    elseif(( 8._8>Lat) .and. (Lat>=0._8)) then
      LetterDesignator = 'N'
    elseif(( 0._8>Lat) .and. (Lat>=-8._8)) then
      LetterDesignator = 'M'
    elseif((-8._8>Lat) .and. (Lat>=-16._8)) then
      LetterDesignator = 'L'
    elseif((-16._8>Lat) .and. (Lat>=-24._8)) then
      LetterDesignator = 'K'
    elseif((-24._8>Lat) .and. (Lat>=-32._8)) then
      LetterDesignator = 'J'
    elseif((-32._8>Lat) .and. (Lat>=-40._8)) then
      LetterDesignator = 'H'
    elseif((-40._8>Lat) .and. (Lat>=-48._8)) then
      LetterDesignator = 'G'
    elseif((-48._8>Lat) .and. (Lat>=-56._8)) then
      LetterDesignator = 'F'
    elseif((-56._8>Lat) .and. (Lat>=-64._8)) then
      LetterDesignator = 'E'
    elseif((-64._8>Lat) .and. (Lat>=-72._8)) then
      LetterDesignator = 'D'
    elseif((-72._8>Lat) .and. (Lat>=-80._8)) then
      LetterDesignator = 'C'
    else
      LetterDesignator = 'Z' ! This is here as an error flag to show that the Latitude is outside the UTM limits
    end if

  end function UTMLetterDesignator



  subroutine UTMtoLL(ReferenceEllipsoid, UTM, G)

    implicit none

    ! Call list
    integer, intent(in)      :: ReferenceEllipsoid
    type (UTM_COORD), intent(in)  :: UTM
    type (GEO_COORD)      :: G

    ! Locals
    real(kind=8)    ::  UTMNorthing
    real(kind=8)    ::  UTMEasting
    character(len=4)  ::  UTMZone
    real(kind=8)    ::  Lat
    real(kind=8)    ::  Long

    ! Locals
    real(kind=8)    ::  k0 = 0.9996_8
    real(kind=8)    ::  a, eccSquared, eccPrimeSquared, e1
    real(kind=8)    ::  N1, T1, C1, R1, D, M
    real(kind=8)    ::  LongOrigin
    real(kind=8)    ::  mu, phi1, phi1Rad
    real(kind=8)    ::  x, y
    integer        ::  ZoneNumber
    character(len=1)  ::  ZoneLetter
    integer        ::  NorthernHemisphere ! 1 for northern hemisphere, 0 for southern

    ! Parameters
    real(kind=8), parameter  ::  PI = 3.14159265_8
    real(kind=8), parameter  ::  FOURTHPI = PI / 4._8
    real(kind=8), parameter  ::  deg2rad = PI / 180._8
    real(kind=8), parameter  ::  rad2deg = 180.0_8 / PI

    ! Transfer
    UTMZone = UTM%Zone
    UTMNorthing = UTM%N
    UTMEasting  = UTM%E

    ! Conversion
    a = DATUMS(ReferenceEllipsoid)%A
    eccSquared = DATUMS(ReferenceEllipsoid)%E2
    e1 = (1._8-sqrt(1._8-eccSquared))/(1._8+sqrt(1._8-eccSquared))
    x = UTMEasting - 500000.0_8 ! remove 500,000 meter offset for longitude
    y = UTMNorthing
    read(UTMZone,"(i3,a1)") ZoneNumber, ZoneLetter
    if(ichar(ZoneLetter)>=ichar('N')) then
      NorthernHemisphere = 1  ! point is in northern hemisphere
    else
      NorthernHemisphere = 0  ! point is in southern hemisphere
      y = y - 10000000.0_8  ! remove 10,000,000 meter offset used for southern hemisphere
    end if
    LongOrigin = (ZoneNumber - 1)*6._8 - 180._8 + 3._8  ! +3 puts origin in middle of zone
    eccPrimeSquared = (eccSquared)/(1._8-eccSquared)
    M = y / k0
    mu = M/(a*(1._8-eccSquared/4._8-3._8*eccSquared**2/64._8-5._8*eccSquared**3/256._8))
    phi1Rad = mu  + (3._8*e1/2._8-27._8*e1**3/32._8)*sin(2._8*mu)  &
      + (21._8*e1**2/16._8-55._8*e1**4/32._8)*sin(4._8*mu)      &
      +(151._8*e1**3/96._8)*sin(6._8*mu)
    phi1 = phi1Rad*rad2deg
    N1 = a/sqrt(1._8-eccSquared*sin(phi1Rad)**2)
    T1 = tan(phi1Rad)**2
    C1 = eccPrimeSquared*cos(phi1Rad)**2
    R1 = a*(1._8-eccSquared)/(1._8-eccSquared*sin(phi1Rad)**2)**1.5
    D = x/(N1*k0)
    Lat = phi1Rad - (N1*tan(phi1Rad)/R1)*(D**2/2._8-(5._8+3._8*T1+10._8*C1-4._8*C1**2-9._8*eccPrimeSquared)*D**4/24._8  &
      +(61._8+90._8*T1+298._8*C1+45._8*T1**2-252._8*eccPrimeSquared-3._8*C1**2)*D**6/720._8)
    Lat = Lat * rad2deg
    Long = (D-(1._8+2._8*T1+C1)*D**3/6._8+(5._8-2._8*C1+28._8*T1-3._8*C1**2+8._8*eccPrimeSquared+24._8*T1**2)  &
      *D**5/120._8)/cos(phi1Rad)
    Long = LongOrigin + Long * rad2deg

    ! Transfer
    G%Lat = Lat
    G%Lon = Long

  end subroutine UTMtoLL

end module ll_utm
