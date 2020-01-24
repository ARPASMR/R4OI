all:
	gfortran --verbose -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/ coords.f90 ll_gb.f90 ll_utm.f90 SpatialStuff.f90 subs19.f90 subsf89.f t2m19.f90  -o t2m19
	gfortran --verbose -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7/ coords.f90 ll_gb.f90 ll_utm.f90 SpatialStuff.f90 subs19.f90 subsf89.f rhtd19.f90 -o rhtd19
	gfortran --verbose -L/usr/lib/gcc/x86_64-redhat-linux/4.4.7 coords.f90 ll_gb.f90 ll_utm.f90 SpatialStuff.f90 subs19.f90 subsf89.f plzln19.f90 -o plzln19

clean:
	rm -f *.o
	rm -f *.mod
	rm -f t2m19
	rm -f rhtd19
	rm -f plzln19

