CFLAGS= -O3 -pthread -DLOG_USE_COLOR

LIBS= -lfftw3_threads -lfftw3 -lm -lrt

INCLUDE= 

OBJS=utils/log.o src/service/raw_samples/samples_server.o 

TARGET=server

all: 
	build

.c.o:
	gcc $(CFLAGS) $(INCLUDE) -c -o $@ $<

build:	$(OBJS)
	gcc -o $(TARGET) $(CFLAGS) $(OBJS) $(LIBS)

clean:
	-rm -rf $(OBJS) $(TARGET)

.PHONY: all build clean