CC = gcc
SIM5LIB = ../../lib

CFLAGS = -I$(SIM5LIB) -L$(SIM5LIB) -Wall -O3 -w -fgnu89-inline
LFLAGS = $(CFLAGS) -lm


default: kerr-orbits

%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

clean:
	rm -f *.o


kerr-orbits-src = \
    $(SIM5LIB)/sim5lib.c \
    kerr-orbits.c \

kerr-orbits-obj = $(kerr-orbits-src:.c=.o)

kerr-orbits: $(kerr-orbits-obj)
	$(CC) $(kerr-orbits-obj) -o $@ $(LFLAGS) 


