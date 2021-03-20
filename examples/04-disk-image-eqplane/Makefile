CC = gcc
SIM5LIB = ../../src

CFLAGS = -I$(SIM5LIB) -L$(SIM5LIB) -Wall -O3 -w -fgnu89-inline
LFLAGS = $(CFLAGS) -lm


default: disk-image

%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

clean:
	rm -f *.o


disk-image-src = \
    $(SIM5LIB)/sim5lib.c \
    disk-image.c \

disk-image-obj = $(disk-image-src:.c=.o)

disk-image: $(disk-image-obj)
	$(CC) $(disk-image-obj) -o $@ $(LFLAGS) 


