LIBS    = -largeArrayDims -lgomp -lm
INCLUDE =
FILES   = sectorDOUBLE.cpp sectorSINGLE.cpp

make:
	mex $(LIBS) $(CFLAGS) $(INCLUDE) scan_convert.cpp $(FILES)
