CFLAGS= -pthread -DLOG_USE_COLOR -Og

LIBS= -lfftw3_threads -lfftw3 -lm -lrt -liniparser 

INCLUDE= 

OBJS=src/clear_freq_search.o src/ini_parser.o src/misc_read_writes.o src/log.o src/clear_frequency_server.o 

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
