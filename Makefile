CC = g++
CFLAGS = -Wall -Wextra

all: nbody

nbody: nbody.cpp
	$(CC) $(CFLAGS) -g -o $@ $^ -lm

clean:
	rm -f nbody

