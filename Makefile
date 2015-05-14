CFLAGS = -fPIC -Ilib -Llib
LFLAGS = -lm -lgsl -lgslcblas

default: sim5lib

clean:
	rm -f *.o lib/*.o bin/*

#%.o: %.c
#	gcc -o $@ $(CFLAGS) -c $<

sim5lib: 
	@mkdir bin
	gcc -Wall -O3 -shared -o bin/sim5lib.so $(CFLAGS) $(LFLAGS) lib/sim5lib.c


sim5lib-dbg: clean
	gcc -Wall -O3 -shared -o bin/sim5lib.so $(CFLAGS) -g -pg $(LFLAGS) lib/sim5lib.c

sim5lib-python: sim5lib
	swig -python -w314 lib/sim5lib.swig
	gcc -c sim5lib.c sim5lib_wrap.c $(CFLAGS) -Ibin -Lbin -I/usr/include/python2.7 $(LFLAGS)
	gcc -shared sim5lib.o sim5lib_wrap.o -o bin/_sim5lib_module.so
