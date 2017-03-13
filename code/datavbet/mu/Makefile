CC = gcc
LIBS = -lm
CFLAGS = -Wall -O2 -march=native

EXECUTABLES = vbet vbetjob

.PHONY: all clean

all: $(EXECUTABLES)

vbet: vbet.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

vbetjob: vbetjob.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o $(EXECUTABLES)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
