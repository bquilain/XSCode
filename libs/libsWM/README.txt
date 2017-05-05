In this folder we find the libraries for the WM software. 

LAST UPDATE : Matthieu L. 2016/10/10

This version has been put to the semi-offline machine before T2K run8.
It corresponds (it should) to ingrid git commit 8ae03c

********  FOLDER INFO  ********

src contains the .cc and .h files

lib_5.22 contains the *.o *.so *Dict.cc *Dict.h *Dict.o *LinkDef.h
produced with ROOT v5r22p00n00

lib_5.34  contains the *.o *.so *Dict.cc *Dict.h *Dict.o *LinkDef.h
produced with ROOT v5r34p09n03


********   IMPORTANT   ********

1- use the libraries with the ROOT version corresponding to your
ROOTSYS variable

Notice: * MonteCarlo requires the version 5.22
	* T2KReWeight (if you use it at some point) will require a
	version > 5.28
		      

2- before you make in the lib_5.22 or lib_5.34 repository, ensure your
ROOTSYS is set to the appropriate ROOT version


3- if you want to compile with your own ROOT version, that's easy. Do
the following:
	mkdir lib_`version`
	cp lib_5.22/Makefile lib_`version`/
	cp lib_5.22/gen_root_so.sh lib_`version`/
	cd lib_`version`
Check your ROOTSYS and make!