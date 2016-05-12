CC = gcc
LIBS = -lm
CFLAGS = -Wall -O2 -march=native

EXECUTABLES = djob dist

.PHONY: all clean

all: $(EXECUTABLES)

dist: dist.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

djob: djob.o
	$(CC) $(CFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o $(EXECUTABLES)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
