AUTOMAKE_OPTIONS=foreign
AM_CFLAGS = -I../sz/include -I../zlib $(GSL_CFLAGS) $(GSL_STATIC)
AM_FCFLAGS = -I../sz/include -I../sz -I../zlib $(GSL_CFLAGS) $(GSL_STATIC)

if FORTRAN
bin_PROGRAMS=testint_compress testint_decompress testfloat_compress testfloat_decompress testdouble_compress testdouble_decompress testfloat_batch_compress testdouble_batch_compress sz testdouble_batch_f testdouble_compress_f testdouble_decompress_f testfloat_CompDecomp_subblock testdouble_CompDecomp_subblock
else
bin_PROGRAMS=testint_compress testint_decompress testfloat_compress testfloat_decompress testdouble_compress testdouble_decompress testfloat_batch_compress testdouble_batch_compress sz testfloat_CompDecomp_subblock testdouble_CompDecomp_subblock
endif

#if GSL
#bin_PROGRAMS += sz_gsl
#endif
	
testint_compress_SOURCES=testint_compress.c
testint_compress_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a -lm
testint_decompress_SOURCES=testint_decompress.c
testint_decompress_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a -lm
testfloat_compress_SOURCES=testfloat_compress.c
testfloat_compress_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a -lm
testfloat_decompress_SOURCES=testfloat_decompress.c
testfloat_decompress_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a -lm
testdouble_compress_SOURCES=testdouble_compress.c
testdouble_compress_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a -lm
testdouble_decompress_SOURCES=testdouble_decompress.c
testdouble_decompress_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a -lm
if FORTRAN
testdouble_batch_f_SOURCES=testdouble_batch_f.f90
testdouble_batch_f_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a -lm
testdouble_compress_f_SOURCES=testdouble_compress_f.f90
testdouble_compress_f_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a -lm
testdouble_decompress_f_SOURCES=testdouble_decompress_f.f90
testdouble_decompress_f_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a -lm
endif
testfloat_batch_compress_SOURCES=testfloat_batch_compress.c
testfloat_batch_compress_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a -lm
testdouble_batch_compress_SOURCES=testdouble_batch_compress.c
testdouble_batch_compress_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a -lm
testfloat_CompDecomp_subblock_SOURCES=testfloat_CompDecomp_subblock.c
testfloat_CompDecomp_subblock_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a -lm
testdouble_CompDecomp_subblock_SOURCES=testdouble_CompDecomp_subblock.c
testdouble_CompDecomp_subblock_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a -lm
if GSL
sz_SOURCES=sz_gsl.c
sz_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a $(GSL_LIBS) $(GSL_HDR) $(GSL_LIB) -lm
else
sz_SOURCES=sz.c
sz_LDADD=../sz/.libs/libSZ.a ../zlib/.libs/libzlib.a -lm
endif
