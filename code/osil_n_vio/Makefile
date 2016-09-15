CC = gcc
LIBS = -lm
CFLAGS = -Wall -O2 -march=native

EXECUTABLES = osil vio vio1 osiljob viojob osilplot vioplot

.PHONY: all clean

all: $(EXECUTABLES)

osil: osil.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

vio: vio.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

vio1: vio_1.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

osiljob: osiljob.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

viojob: viojob.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

osilplot: osilplot.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

vioplot: vioplot.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o $(EXECUTABLES)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
