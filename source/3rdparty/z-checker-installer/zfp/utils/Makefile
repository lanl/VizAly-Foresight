include ../Config

ZCPATH = ../../Z-checker/zc-install
ZCFLAG  = -I$(ZCPATH)/include $(ZCPATH)/lib/libzc.a

TARGET = ../bin/zfp
TARGET2 = ../bin/zfp-zc

all: $(TARGET) $(TARGET2)

$(TARGET): zfp.c ../lib/libzfp.a
	mkdir -p ../bin
	$(CC) $(CFLAGS) zfp.c -L../lib -lzfp -lm -o $(TARGET) 

$(TARGET2): zfp-zc.c ../lib/libzfp.a
	mkdir -p ../bin
	$(CC) $(CFLAGS) zfp-zc.c -L../lib -lzfp -lm -o $(TARGET2) $(ZCFLAG)

clean:
	rm -f $(TARGET) $(TARGET2) fields.o
