CFLAGS = -Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas -O3 -fPIC -Ilib -Llib
#LFLAGS = -lm -lgsl -lgslcblas
LFLAGS = -lm

default: lib

clean:
	@rm -f *.o lib/*.o bin/*
	@rm -f lib/lib/*_wrap.c

#%.o: %.c
#	gcc -o $@ $(CFLAGS) -c $<

all: lib python export
	@echo "done"

lib: lib-clean
	gcc -c lib/sim5lib.c -o lib/sim5lib.o $(CFLAGS) $(LFLAGS)

cuda: lib-clean
	nvcc -arch=sm_35 -Ilib -Llib -O3 -dc lib/sim5lib.cu

lib-clean:
	@echo "Cleaning..."
	@rm -f lib/*.o

python: lib
	swig -python -w314 lib/sim5lib.swig
	mv lib/sim5lib_module.py bin/sim5lib_module.py
	sed -i "s/'_sim5lib_module'/'sim5lib'/g" bin/sim5lib_module.py
	sed -i "s/_sim5lib_module/sim5lib/g" lib/sim5lib_wrap.c
	gcc -c lib/sim5lib_wrap.c -o lib/sim5lib_wrap.o $(CFLAGS) -I/usr/include/python2.7 $(LFLAGS)
	gcc -shared lib/sim5lib.o lib/sim5lib_wrap.o $(CFLAGS) $(LFLAGS) -o bin/sim5lib.so
#	rm -f ./*.o lib/lib/*_wrap.c


export: 
	echo 'Compiling for export'
	@rm -f bin/sim5lib.h
	@for i in `cat lib/sim5lib.h | grep -e '^\#include ".*.h"$$' | sed -n  's/.*"\(.*\)"/\1/p'`; do cat lib/$$i >> bin/sim5lib.h; done
	@sed -i 's/[ \t]*$$//;/^\/\//d;s/\/\/.*$$//' bin/sim5lib.h
	@rm -f bin/sim5lib.c
	@echo "#include \"sim5lib.h\"" > bin/sim5lib.c
	@for i in `cat lib/sim5lib.c | grep -e '^\#include ".*.c"$$' | sed -n  's/.*"\(.*\)"/\1/p'`; do cat lib/$$i >> bin/sim5lib.c; done
	@sed -i 's/[ \t]*$$//;/^\/\//d;s/\/\/.*$$//;s/    / /' bin/sim5lib.c
	gcc -c bin/sim5lib.c -o bin/sim5lib.o $(CFLAGS) $(LFLAGS)
	@rm bin/sim5lib.o


debug: lib
	gcc -shared -o bin/sim5lib.so $(CFLAGS) -O0 -g -pg $(LFLAGS) lib/sim5lib.c


test: lib
	@rm -f bin/sim5lib-tests
	gcc -c lib/sim5unittests.c -o lib/sim5unittests.o $(CFLAGS) $(LFLAGS) 
	gcc lib/sim5unittests.o lib/sim5lib.o -o bin/sim5lib-tests $(CFLAGS) $(LFLAGS) 
	if [ -e bin/sim5lib-tests ]; then bin/sim5lib-tests; fi

.PHONY: doc

doc:
	doxygen doxygen.cfg
	xsltproc doxygen.xsl doc/xml/index.xml > doc/sim5lib-doc.md
