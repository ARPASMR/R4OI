! Module "SpatialStuff" - Routines supporting the encoding and decoding of spatial
!                         data according to (and enforcing of) Lombardy conventions.
!
! By: Mauri Favaron (in reality, really a pun: I've inherited almost all of this code,
!                    and my contribution has been mainly in tidying it. The real authors
!                    are giants like Cristian Lussana, Maria Ranci, Francesco Uboldi,
!					 Marta Salvati and Umberto Pellegrini. I have not had a way to credit
!					 any individual contributions, but anyway, a VERY BIG THANK YOU ALL!
!					 Francesco especially, whose "impromptu" on code is pervasive.)
! 
MODULE SpatialStuff

	USE coords
	USE ll_utm
	USE gauss_boaga
	
	IMPLICIT NONE
	
	! I'm defining all "PRIVATE", so any name defined in this module will
	! not be seen outside of it unless explicitly named "PUBLIC". This, to
	! help preventing name-space cluttering (I don't know whether you have
	! my same problem, but as programs grow big my mono-neuron tends to
	! forget all names looking unimportant, regardless they really are).
	PRIVATE
	
	! Public interface: Here is where names accessible to everyone are
	PUBLIC	:: GetUTMOrography
	PUBLIC	:: GetGBOrography
	PUBLIC	:: SampleField
	
CONTAINS

	! *****************************************
	! * Support for standard orographic files *
	! * (also known as Digital Terrain Model, *
	! * or DTM for short.                     *
	! *****************************************
	
	! Standard Lombardy orographic files come in two forms, depending on reference system.
	! One is the UTM (fuse 32), and the other is Gauss-Boaga (Western fuse). Both forms are
	! stored as binary files, composed by a sequence of equal size "fields". Each field
	! contains values sampled over a regularly-shaped grid, which is rectangular in its
	! reference system. (Of course UTM32 and GB-W are £rectangular in different ways": if
	! you project a rectangular grid in the first to the other reference, you obtain a grid
	! which is not rectangular any longer (although on a region small as Lombardy you would
	! hardly notice a difference, if you watch a map of it on an A4 sheet.
	!
	! It is worth noticing the actual correspondence between terrain coordinates and within-field
	! positional index is *not* stored into the file, but is known in advance (and hard-wired into
	! DTM read routines GetUTMOrography and GetGBOrography).
	!
	! To date, the field contents is a bit different in the two references. The following table applies.
	!
	!  _____________________________________________________________________________________________
	! | Field #  | Gauss-Boaga                            | UTM                                     |
	! |----------|----------------------------------------|-----------------------------------------|
	! |     1    |                                The DTM, in meters.                               |
	! |----------|----------------------------------------|-----------------------------------------|
	! |     2    |       A boolean mask, indicating a grid point belongs to Lombardy or not.        |
	! |----------|----------------------------------------|-----------------------------------------|
	! |     3    |         Slope along 'x' axis           |           "Urban" boolean mask          |
	! |----------|----------------------------------------|-----------------------------------------|
	! |     4    |         Slope along 'y' axis           |                  n/a                    |
	! |----------|----------------------------------------|-----------------------------------------|
	! |     5    |         "Urban" boolean mask           |                  n/a                    |
	! |__________|________________________________________|_________________________________________|
	!
	! Another important difference is, the two grid sizes differ in node number and spacing:
	!
	!   - Gauss-Boaga DTM: 1500m spacing, Nx = 177, Ny = 174.
	!
	!   - UTM DTM:         1000m spacing, Nx = 251, Ny = 249.
	!
	! These differences are made for the most part transparent by the two reading routines.
	! The slopes along X and Y are ignored.
	
	! A word of caution: I've written this  ignored.routine after one of Francesco. In so doing I've modified
	! theoriginal quite conservatively, but for a single point: I added a field to Lombardy standard UTM
	! orography, on last position (record no.3), holding the urban incidence matrix. This is where
	! 'g4uhi' comes from.
	SUBROUTINE GetUTMOrography(sDTM, ii, g4x, g4y, g4lon, g4lat, g4hgt, g4msk, g4uhi)
	
		! Standard size of Lombardy UTM orography data & mask
		INTEGER, PARAMETER		:: INX=251, INY=249, IM=INX*INY
		
		! Routine arguments
		CHARACTER(LEN=*), INTENT(IN)		:: sDTM										! Data file
		INTEGER, INTENT(OUT)				:: ii											! Data set size
		REAL(4), DIMENSION(IM), INTENT(OUT)	:: g4x, g4y, g4lon, g4lat, g4hgt, g4msk, g4uhi	! UTM and corresponding LL data, as read from standard Lombardy file
		
		! Local variables
		REAL(8) 			:: xmin, ymin, dx, dy
		INTEGER				:: ix, jy, i
		TYPE(GEO_COORD)		:: tLL   ! tLL%Lat  tLL%Lon        (r8 r8)
		TYPE(UTM_COORD)		:: tUTM  ! tUTM%Zone tUTM%N tUTM%E (c4 r8 r8)
		INTEGER				:: iErrCode

		! Make number of data visible
		ii = IM

		! Set grid SW point in UTM and define spacing information
		xmin=  436.d3
		ymin= 4918.d3
		dx  =    1.d3
		dy  =    1.d3
		
		! Generate all grid points in UTM
		DO i = 1, ii
			jy= (i-1) /inx +1
			ix= i -(jy-1) *inx
			tUTM%E = xmin+ (ix -1) *dx 
			tUTM%N = ymin+ (jy -1) *dy
			tUTM%Zone = '+32T'
			CALL UTMtoLL(INTERNATIONAL, tUTM, tLL)
			g4x(i)= sngl(tUTM%E)
			g4y(i)= sngl(tUTM%N)
			g4lon(i)= tLL%Lon
			g4lat(i)= tLL%Lat
		END DO
		
		! Read orographic data and dispatch them to the appropriate variables
		OPEN(3, FILE= sDTM, FORM= 'UNFORMATTED', ACCESS='DIRECT', RECL=4*IM)
		READ(3,rec=1) g4hgt
		READ(3,rec=2) g4msk
		READ(3,rec=3,IOSTAT=iErrCode) g4uhi
		CLOSE(3)
		
	END SUBROUTINE GetUTMOrography
	

	SUBROUTINE GetGBOrography(sDTM, ii, g4x, g4y, g4lon, g4lat, g4hgt, g4msk, g4uhi)
	
		! Standard size of Lombardy UTM orography data & mask
		INTEGER, PARAMETER		:: INX=177, INY=174, IM=INX*INY
		
		! Routine arguments
		CHARACTER(LEN=*), INTENT(IN)		:: sDTM										! Data file
		INTEGER, INTENT(OUT)				:: ii											! Data set size
		REAL(4), DIMENSION(IM), INTENT(OUT)	:: g4x, g4y, g4lon, g4lat, g4hgt, g4msk, g4uhi	! Gauss-Boaga and corresponding LL data, as read from standard Lombardy file

		
		! Local variables
		REAL(8) 			:: xmin, ymin, dx, dy
		INTEGER				:: ix, jy, i
		TYPE(GEO_COORD)		:: tLL   ! tLL%Lat  tLL%Lon        (r8 r8)
		TYPE(GAUSS_COORD)	:: tGB   ! tGB%Zone tUTM%N tUTM%E (c4 r8 r8)

		! Make number of data visible
		ii = IM

		! Set grid SW point in UTM and define spacing information
		xmin= 1436301.375d0
		ymin= 4916704.500d0
		dx  =    1500.d0
		dy  =    1500.d0
		dy  =    1500.d0
		dy  =    1500.d0
		
		! Generate all grid points in UTM
		DO i = 1, ii
			jy    = (i-1) /inx +1
			ix    = i -(jy-1) *inx
			tGB%E = xmin + (ix - 1)*dx
			tGB%N = ymin + (jy - 1)*dy
			tGB%S = WEST
			tLL   = gb_to_ll(tGB)
			g4x(i)= sngl(tGB%E)
			g4y(i)= sngl(tGB%N)
			g4lon(i)= tLL%Lon
			g4lat(i)= tLL%Lat
		END DO
		
		! Read orographic data and dispatch them to the appropriate variables
		OPEN(3, FILE= sDTM, FORM= 'UNFORMATTED', ACCESS='DIRECT', RECL=4*IM)
		READ(3,rec=1) g4hgt
		READ(3,rec=2) g4msk
		READ(3,rec=6) g4uhi
		CLOSE (3)
		
	END SUBROUTINE GetGBOrography
	

	! Resamples a field on a set of points based on closeness (no interpolation performed).
	!
	! A warning: this routine is *very* slow, mainly by the way I've written it. Sure is very
	! safe. But if you want not to waste five good minutes of your life (on a dual-core
	! Intel Atom), you may imagine using something lighter.
	SUBROUTINE SampleField(gLatIn, gLonIn, gFieldIn, gLatOut, gLonOut, gFieldOut)
	
		! Routine arguments
		REAL(4), DIMENSION(:), INTENT(IN)	:: gLatIn
		REAL(4), DIMENSION(:), INTENT(IN)	:: gLonIn
		REAL(4), DIMENSION(:), INTENT(IN)	:: gFieldIn
		REAL(4), DIMENSION(:), INTENT(IN)	:: gLatOut
		REAL(4), DIMENSION(:), INTENT(IN)	:: gLonOut
		REAL(4), DIMENSION(:), INTENT(OUT)	:: gFieldOut
		
		! Locals
		INTEGER	:: iNumDataIn
		INTEGER	:: iNumDataOut
		INTEGER	:: iPos
		INTEGER	:: i, j
		REAL	:: rMinDistance, rDistance
		
		! Build output field by locating the closest original points' values
		iNumDataIn  = SIZE(gLatIn)
		iNumDataOut = SIZE(gLatOut)
		DO i = 1, iNumDataOut
			rMinDistance = HUGE(rMinDistance)
			DO j = 1, iNumDataIn
				rDistance    = SQRT((gLatIn(j)-gLatOut(i))**2+(gLonIn(j)-gLonOut(i))**2)
				IF(rDistance < rMinDistance) THEN
					rMinDistance = rDistance
					iPos         = j
				END IF
			END DO
			gFieldOut(i) = gFieldIn(iPos)
		END DO
		
	END SUBROUTINE SampleField
	
END MODULE SpatialStuff


! Module test procedure (Mauri: maybe worth a separate file?)
!PROGRAM TestSpatialStuff
!
!	USE SpatialStuff
!	
!	IMPLICIT NONE
!	
!	INTEGER, PARAMETER			:: INX=177, INY=174	! Gauss-Boaga grid
!	INTEGER, PARAMETER			:: JNX=251, JNY=249	! UTM grid
!	REAL(4), DIMENSION(INX*INY)	:: g41x, g41y, g41lon, g41lat, g41hgt, g41msk, g41uhi
!	REAL(4), DIMENSION(JNX*JNY)	:: g42x, g42y, g42lon, g42lat, g42hgt, g42msk, g42uhi
!	INTEGER						:: iRetCode, ii
!	
!	CALL GetGBOrography("./topography_1500.dat", ii, g41x, g41y, g41lon, g41lat, g41hgt, g41msk, g41uhi)
!	PRINT *,ii
!	PRINT *,'Data for GB file read'
!	PRINT *,"Hgt (m): ", MINVAL(g41hgt), MAXVAL(g41hgt)
!	PRINT *,"Urb (m): ", MINVAL(g41uhi), MAXVAL(g41uhi)
!	PRINT *,'Num.urb. ', COUNT(g41uhi > 0)
!	
!	CALL GetUTMOrography("./oromaskutm.dat", ii, g42x, g42y, g42lon, g42lat, g42hgt, g42msk, g42uhi)
!	PRINT *,'Data for UTM file read'
!	PRINT *,"Hgt (m): ", MINVAL(g42hgt), MAXVAL(g42hgt)
!	
!	CALL SampleField(g41lat, g41lon, g41uhi, g42lat, g42lon, g42uhi)
!	PRINT *,"Urb (m): ", MINVAL(g42uhi), MAXVAL(g42uhi)
!	
!	OPEN(10, FILE="./oromaskutm_mauri.dat", FORM= 'UNFORMATTED', ACCESS='DIRECT', RECL=4*JNX*JNY)
!	WRITE(10,REC=1) g42hgt
!	WRITE(10,REC=2) g42msk
!	WRITE(10,REC=3) g42uhi
!	CLOSE(10)
!
!	CALL GetUTMOrography("./oromaskutm_mauri.dat", ii, g42x, g42y, g42lon, g42lat, g42hgt, g42msk, g42uhi)
!	PRINT *,'Data for UTM file read'
!	PRINT *,"Hgt (m): ", MINVAL(g42hgt), MAXVAL(g42hgt)
!	PRINT *,"Urb (m): ", MINVAL(g42uhi), MAXVAL(g42uhi)
!	
!END PROGRAM TestSpatialStuff

