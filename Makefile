CFLAGS = -Wextra -O3 -fPIC -Ilib -Llib
LFLAGS = -lm -lgsl -lgslcblas

default: sim5lib

clean:
	@rm -f *.o lib/*.o bin/*
	@rm -f lib/lib/*_wrap.c

#%.o: %.c
#	gcc -o $@ $(CFLAGS) -c $<

sim5lib: sim5lib-c sim5lib-python
	@echo "done"

sim5lib-c: 
	@rm -f lib/*.o
	gcc -c lib/sim5lib.c -o lib/sim5lib.o $(CFLAGS) $(LFLAGS)


sim5lib-python: sim5lib-c
	swig -python -w314 lib/sim5lib.swig
	mv lib/sim5lib_module.py bin/sim5lib_module.py
	sed -i "s/'_sim5lib_module'/'sim5lib'/g" bin/sim5lib_module.py
	sed -i "s/_sim5lib_module/sim5lib/g" lib/sim5lib_wrap.c
	gcc -c lib/sim5lib_wrap.c -o lib/sim5lib_wrap.o $(CFLAGS) -I/usr/include/python2.7 $(LFLAGS)
	gcc -shared lib/sim5lib.o lib/sim5lib_wrap.o $(CFLAGS) $(LFLAGS) -o bin/sim5lib.so
#	rm -f ./*.o lib/lib/*_wrap.c


sim5lib-export: 
	echo Compiling for export
	@rm -f lib/sim5lib-full.h
	@for i in `cat lib/sim5lib.h | grep -e '^\#include ".*.h"$$' | sed -n  's/.*"\(.*\)"/\1/p'`; do cat lib/$$i >> lib/sim5lib-full.h; done
	@sed -i 's/[ \t]*$$//;/^\/\//d;s/\/\/.*$$//' lib/sim5lib-full.h
	@rm -f lib/sim5lib-full.c
	@echo "#include \"sim5lib-full.h\"" > lib/sim5lib-full.c
	@for i in `cat lib/sim5lib.c | grep -e '^\#include ".*.c"$$' | sed -n  's/.*"\(.*\)"/\1/p'`; do cat lib/$$i >> lib/sim5lib-full.c; done
	@sed -i 's/[ \t]*$$//;/^\/\//d;s/\/\/.*$$//;s/    / /' lib/sim5lib-full.c
	gcc -c lib/sim5lib-full.c -o lib/sim5lib-full.o $(CFLAGS) $(LFLAGS)


sim5lib-dbg: sim5lib-c
	gcc -shared -o bin/sim5lib.so $(CFLAGS) -O0 -g -pg $(LFLAGS) lib/sim5lib.c


sim5lib-tests: sim5lib-c
	gcc -c lib/sim5unittests.c -o lib/sim5unittests.o $(CFLAGS) $(LFLAGS) 
	gcc lib/sim5unittests.o lib/sim5lib.o -o bin/sim5lib-tests $(CFLAGS) $(LFLAGS) 

