CFLAGS = -Wall -Wextra -Wno-unused-parameter -Wno-unknown-pragmas -O3 -fPIC -Isrc -Lsrc -std=gnu11 -fgnu89-inline
LFLAGS = -lm

CC=gcc

# prerequisities
EXISTS_SWIG := $(shell command -v swig 2> /dev/null)
EXISTS_PYDEV := $(shell command -v python3-config 2> /dev/null)
EXISTS_DOXYGEN := $(shell command -v doxygen 2> /dev/null)
EXISTS_XSLTPROC := $(shell command -v xsltproc 2> /dev/null)
EXISTS_NVCC := $(shell command -v nvcc 2> /dev/null)
PYDEV_INC := $(shell python3-config --includes)


default: all

all: lib export
	@echo "done"

clean:
	@rm -f *.o src/*.o lib/*
	@rm -f src/*_wrap.c

#%.o: %.c
#	$(CC) -o $@ $(CFLAGS) -c $<


lib: lib-clean
	@[ -f src/sim5config.h ] || cp src/sim5config.h.default src/sim5config.h
	$(CC) -c src/sim5lib.c -o src/sim5lib.o $(CFLAGS) $(LFLAGS)

cuda: lib-clean nvcc-check
	nvcc -arch=sm_35 -Isrc -Lsrc -O3 -dc src/sim5lib.cu

lib-clean:
	@echo "Cleaning..."
	@rm -f src/*.o

python: lib
ifndef EXISTS_SWIG
	$(error "FAILED PREREQUISITY: SWIG has not been found on the system PATH, install it with apt install swig")
endif
ifndef EXISTS_PYDEV
	$(error "FAILED PREREQUISITY: header files for Python are not installed, install it with apt install python3-dev")
endif
	@mkdir -p lib
	swig -python -py3 -w314 -w301 src/sim5lib.swig
	@mv src/sim5lib.py python/sim5lib.py
	@sed -i "s/'_sim5lib'/'sim5lib'/g" python/sim5lib.py
	@sed -i "s/_sim5lib/sim5lib/g" src/sim5lib_wrap.c
#	cat python/sim5*.py >> lib/sim5lib.py
	$(CC) -c src/sim5lib_wrap.c -o src/sim5lib_wrap.o $(CFLAGS) $(PYDEV_INC) $(LFLAGS) -w
	$(CC) -shared src/sim5lib.o src/sim5lib_wrap.o $(CFLAGS) $(LFLAGS) -o lib/sim5lib.so
	patch python/sim5lib.py python/sim5lib.py.patch
	@rm -f src/*_wrap.*

export:
	echo 'Compiling for export'
	@mkdir -p lib
	@rm -f lib/sim5lib.h
	@for i in `cat src/sim5lib.h | grep -e '^\#include ".*.h"$$' | sed -n  's/.*"\(.*\)"/\1/p'`; do cat src/$$i >> lib/sim5lib.h; done
	@sed -i 's/[ \t]*$$//;/^[ \t]*\/\//d;s/\/\/.*$$//' lib/sim5lib.h
	@rm -f lib/sim5lib.c
	@echo "#include \"sim5lib.h\"" > lib/sim5lib.c
	@for i in `cat src/sim5lib.c | grep -e '^\#include ".*.c"$$' | sed -n  's/.*"\(.*\)"/\1/p'`; do cat src/$$i >> lib/sim5lib.c; done
	@sed -i 's/[ \t]*$$//;/^\/\//d;s/\/\/.*$$//;s/    / /' lib/sim5lib.c
	$(CC) -c lib/sim5lib.c -o lib/sim5lib.o $(CFLAGS) $(LFLAGS)


debug: lib
	$(CC) -shared -o lib/sim5lib.so $(CFLAGS) -O0 -g -pg $(LFLAGS) src/sim5lib.c


test: lib
	@rm -f bin/sim5lib-tests
	$(CC) -c src/sim5unittests.c -o src/sim5unittests.o $(CFLAGS) $(LFLAGS)
	$(CC) src/sim5unittests.o src/sim5lib.o -o bin/sim5lib-tests $(CFLAGS) $(LFLAGS)
	if [ -e bin/sim5lib-tests ]; then bin/sim5lib-tests; fi

.PHONY: doc python


doc: doxygen-check
	doxygen doc/doxygen.cfg
	xsltproc doc/doxygen.xsl doc/xml/index.xml > doc/sim5lib-doc.md


#nvcc-check:
#ifndef EXISTS_NVCC
#    $(error "FAILED PREREQUISITY: NVCC has not been found on the system PATH, install the NVIDIA CUDA SDK")
#endif

