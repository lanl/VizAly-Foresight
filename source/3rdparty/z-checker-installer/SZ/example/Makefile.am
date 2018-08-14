AUTOMAKE_OPTIONS=foreign
AM_CFLAGS = -lm -I../sz/include -I../zlib
AM_FCFLAGS = -I../sz/include -I../sz -I../zlib
if FORTRAN
bin_PROGRAMS=testfloat_compress testfloat_decompress testdouble_compress testdouble_decompress testfloat_batch_compress testdouble_batch_compress testdouble_batch_f testdouble_compress_f testdouble_decompress_f
else
bin_PROGRAMS=testfloat_compress testfloat_decompress testdouble_compress testdouble_decompress testfloat_batch_compress testdouble_batch_compress
endif
testfloat_compress_SOURCES=testfloat_compress.c
testfloat_compress_LDADD=../sz/.libs/libsz.a ../zlib/.libs/libzlib.a
testfloat_decompress_SOURCES=testfloat_decompress.c
testfloat_decompress_LDADD=../sz/.libs/libsz.a ../zlib/.libs/libzlib.a
testdouble_compress_SOURCES=testdouble_compress.c
testdouble_compress_LDADD=../sz/.libs/libsz.a ../zlib/.libs/libzlib.a
testdouble_decompress_SOURCES=testdouble_decompress.c
testdouble_decompress_LDADD=../sz/.libs/libsz.a ../zlib/.libs/libzlib.a
if FORTRAN
testdouble_batch_f_SOURCES=testdouble_batch_f.f90
testdouble_batch_f_LDADD=../sz/.libs/libsz.a ../zlib/.libs/libzlib.a
testdouble_compress_f_SOURCES=testdouble_compress_f.f90
testdouble_compress_f_LDADD=../sz/.libs/libsz.a ../zlib/.libs/libzlib.a
testdouble_decompress_f_SOURCES=testdouble_decompress_f.f90
testdouble_decompress_f_LDADD=../sz/.libs/libsz.a ../zlib/.libs/libzlib.a
endif
testfloat_batch_compress_SOURCES=testfloat_batch_compress.c
testfloat_batch_compress_LDADD=../sz/.libs/libsz.a ../zlib/.libs/libzlib.a
testdouble_batch_compress_SOURCES=testdouble_batch_compress.c
testdouble_batch_compress_LDADD=../sz/.libs/libsz.a ../zlib/.libs/libzlib.a
