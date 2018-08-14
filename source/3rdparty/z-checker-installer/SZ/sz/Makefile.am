#AM_CFLAGS = -I./include -I../zlib
#LDFLAGS=-fPIC -shared

AUTOMAKE_OPTIONS=foreign
if FORTRAN
include_HEADERS=include/CompressElement.h include/DynamicDoubleArray.h include/rw.h\
		include/dictionary.h include/DynamicFloatArray.h include/VarSet.h include/sz.h include/Huffman.h\
		include/DynamicByteArray.h include/DynamicIntArray.h include/iniparser.h include/TightDataPointStorageD.h include/TightDataPointStorageF.h sz.mod rw.mod

lib_LTLIBRARIES=libsz.la
libsz_la_CFLAGS=-g -I./include -I../zlib/
libsz_la_LDFLAGS = -version-info  1:4:0
libsz_la_LIDADD=../zlib/.libs/libzlib.a
libsz_la_SOURCES=src/ByteToolkit.c src/dataCompression.c src/DynamicIntArray.c src/iniparser.c src/szf.c \
		src/CompressElement.c src/DynamicByteArray.c src/rw.c src/TightDataPointStorageD.c src/TightDataPointStorageF.c \
		src/conf.c src/DynamicDoubleArray.c src/rwf.c src/TypeManager.c \
		src/dictionary.c src/DynamicFloatArray.c src/VarSet.c src/test_zlib.c src/Huffman.c src/sz_float.c src/sz_double.c src/sz.c \
		src/TightDataPointStorageF_pwr.c src/TightDataPointStorageD_pwr.c src/sz_float_pwr.c src/sz_double_pwr.c \
		src/sz_interface.F90 src/rw_interface.F90
libsz_la_LINK=$(AM_V_CC)$(LIBTOOL) --tag=FC --mode=link $(FCLD) $(libsz_la_CFLAGS) -O3 $(libsz_la_LDFLAGS) -o $(lib_LTLIBRARIES)
else
include_HEADERS=include/CompressElement.h include/DynamicDoubleArray.h include/rw.h\
		include/dictionary.h include/DynamicFloatArray.h include/VarSet.h include/sz.h include/Huffman.h\
		include/DynamicByteArray.h include/DynamicIntArray.h include/iniparser.h include/TightDataPointStorageD.h include/TightDataPointStorageF.h

lib_LTLIBRARIES=libsz.la
libsz_la_CFLAGS=-g -I./include -I../zlib/
libsz_la_LDFLAGS = -version-info  1:4:0
libsz_la_LIDADD=../zlib/.libs/libzlib.a
libsz_la_SOURCES=src/ByteToolkit.c src/dataCompression.c src/DynamicIntArray.c src/iniparser.c\
		src/CompressElement.c src/DynamicByteArray.c src/rw.c src/TightDataPointStorageD.c src/TightDataPointStorageF.c \
		src/conf.c src/DynamicDoubleArray.c src/TypeManager.c \
		src/dictionary.c src/DynamicFloatArray.c src/VarSet.c src/test_zlib.c src/Huffman.c src/sz_float.c src/sz_double.c src/sz.c \
		src/TightDataPointStorageF_pwr.c src/TightDataPointStorageD_pwr.c src/sz_float_pwr.c src/sz_double_pwr.c
libsz_la_LINK=$(AM_V_CC)$(LIBTOOL) --tag=CC --mode=link $(CCLD) $(libsz_la_CFLAGS) -O3 $(libsz_la_LDFLAGS) -o $(lib_LTLIBRARIES)
endif
