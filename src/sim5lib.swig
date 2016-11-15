%module sim5lib_module
%{
    #define DEVICEFUNC

    /* the resulting C file should be built as a python extension */
    //#define SWIG_FILE_WITH_INIT
    /*  Includes the header in the wrapper code */
    #include "sim5lib.h"
%}

%include cpointer.i
%include carrays.i
%array_class(int, intArray);
%array_class(double, doubleArray);
%pointer_class(int, intp);
%pointer_class(double, doublep);
//%pointer_functions(double ,double_ptr);  // this vill allow to access content of general 'double *' variables
%array_functions(double, double_array);  // will allow to access content of general 'double *' arrays

// Some array helpers
%inline %{
    // create a connection array helper
    // usage in python: G = connection()
    double (*connection())[4][4] {
        return (double (*)[4][4]) calloc(64,sizeof(double));
    }
%}

%include "sim5config.h"
%include "sim5const.h"
%include "sim5elliptic.h"
%include "sim5kerr.h"
%include "sim5raytrace.h"
%include "sim5kerr-geod.h"
%include "sim5disk-nt.h"

// vector() routine has been added to sim5kerr.c
%pythoncode %{
def sim5vector(components):
    v  = doubleArray(4)
    v[0] = components[0]
    v[1] = components[1]
    v[2] = components[2]
    v[3] = components[3]
    return v
%}
