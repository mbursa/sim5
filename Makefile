CFLAGS = -Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas -O3 -fPIC -Isrc -Lsrc -std=gnu11 -fgnu89-inline
LFLAGS = -lm

CC=gcc

default: lib

clean:
	@rm -f *.o src/*.o bin/*
	@rm -f src/*_wrap.c

#%.o: %.c
#	$(CC) -o $@ $(CFLAGS) -c $<

all: lib python export
	@echo "done"

lib: lib-clean
	$(CC) -c src/sim5lib.c -o src/sim5lib.o $(CFLAGS) $(LFLAGS)

cuda: lib-clean
	nvcc -arch=sm_35 -Isrc -Lsrc -O3 -dc src/sim5lib.cu

lib-clean:
	@echo "Cleaning..."
	@rm -f src/*.o

python: lib
	swig -python -w314 src/sim5lib.swig
	mv src/sim5lib.py lib/sim5lib.py
	sed -i "s/'_sim5lib'/'sim5lib'/g" lib/sim5lib.py
	sed -i "s/_sim5lib/sim5lib/g" src/sim5lib_wrap.c
	cat py/sim5*.py >> bin/sim5lib.py
	$(CC) -c src/sim5lib_wrap.c -o src/sim5lib_wrap.o $(CFLAGS) -I/usr/include/python2.7 $(LFLAGS)
	$(CC) -shared src/sim5lib.o src/sim5lib_wrap.o $(CFLAGS) $(LFLAGS) -o lib/sim5lib.so
	rm -f src/*_wrap.*


export:
	echo 'Compiling for export'
	@rm -f bin/sim5lib.h
	@for i in `cat src/sim5lib.h | grep -e '^\#include ".*.h"$$' | sed -n  's/.*"\(.*\)"/\1/p'`; do cat src/$$i >> bin/sim5lib.h; done
	@sed -i 's/[ \t]*$$//;/^[ \t]*\/\//d;s/\/\/.*$$//' bin/sim5lib.h
	@rm -f bin/sim5lib.c
	@echo "#include \"sim5lib.h\"" > bin/sim5lib.c
	@for i in `cat src/sim5lib.c | grep -e '^\#include ".*.c"$$' | sed -n  's/.*"\(.*\)"/\1/p'`; do cat src/$$i >> bin/sim5lib.c; done
	@sed -i 's/[ \t]*$$//;/^\/\//d;s/\/\/.*$$//;s/    / /' bin/sim5lib.c
	$(CC) -c bin/sim5lib.c -o bin/sim5lib.o $(CFLAGS) $(LFLAGS)
	@rm bin/sim5lib.o


debug: lib
	$(CC) -shared -o bin/sim5lib.so $(CFLAGS) -O0 -g -pg $(LFLAGS) src/sim5lib.c


test: lib
	@rm -f bin/sim5lib-tests
	$(CC) -c src/sim5unittests.c -o src/sim5unittests.o $(CFLAGS) $(LFLAGS)
	$(CC) src/sim5unittests.o src/sim5lib.o -o bin/sim5lib-tests $(CFLAGS) $(LFLAGS)
	if [ -e bin/sim5lib-tests ]; then bin/sim5lib-tests; fi

.PHONY: doc

doc:
	doxygen doc/doxygen.cfg
	xsltproc doc/doxygen.xsl doc/xml/index.xml > doc/sim5lib-doc.md
