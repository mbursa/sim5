/*
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "sim5utils.h"
#include "sim5math.h"
*/

void gprintf(FILE* file, const char *template, ...)
//**************************************************
{
    #ifndef CUDA
    va_list ap;

    va_start (ap, template);
    vfprintf (stderr, template, ap);
    va_end (ap);

    va_start (ap, template);
    vfprintf (file, template, ap);
    va_end (ap);
    #endif
}


void warning(const char *template, ...)
//**************************************************
{
    #ifndef CUDA
    va_list ap;
    fprintf (stderr, "ERROR: ");
    va_start (ap, template);
    vfprintf (stderr, template, ap);
    va_end (ap);
	fprintf(stderr, "\n");
    #endif
}


void error(const char *template, ...)
//**************************************************
{
    #ifndef CUDA
    va_list ap;
    fprintf (stderr, "ERROR: ");
    va_start (ap, template);
    vfprintf (stderr, template, ap);
    va_end (ap);
	fprintf(stderr, "\n");
	exit(-1);
    #endif
}


void sort_array(double *array, int N)
//**************************************************
{
    int compare_func(const void *x, const void *y) {
      return (*(double*)x - *(double*)y);
    }
    qsort(array, N, sizeof(double), compare_func);
}


void sort_array_f(float *array, int N)
//**************************************************
{
    int compare_func(const void *x, const void *y) {
      return (*(float*)x - *(float*)y);
    }
    qsort(array, N, sizeof(float), compare_func);
}


void* array_alloc(size_t capacity, size_t element_size)
{
    size_t size_total = element_size*capacity    // space allocated for elements
                      + sizeof(size_t)           // space for info about size of stored elements 
                      + sizeof(size_t)           // space for info about allocated capacity
                      + sizeof(size_t);          // space for info about stored element count
    
    void* ptr = malloc(size_total);
    memset(ptr, '\0', size_total);

    // store element size    
    *(size_t*)(ptr) = element_size;
    ptr += sizeof(size_t);

    // store allocated capacity
    *(size_t*)(ptr) = capacity;
    ptr += sizeof(size_t);

    // store element count
    *(size_t*)(ptr) = 0;
    ptr += sizeof(size_t);

    return ptr;
}


void array_free(void* ptr)
{
    ptr -= 3*sizeof(size_t);
    free(ptr);
}


void* array_realloc(void* array, size_t new_capacity)
{
    size_t* count = (size_t*)(array-1*sizeof(size_t));
    size_t* capa  = (size_t*)(array-2*sizeof(size_t));
    size_t* esize = (size_t*)(array-3*sizeof(size_t));
    
    size_t size_total = new_capacity * (*esize)    // space allocated for elements
                      + sizeof(size_t)             // space for info about size of stored elements 
                      + sizeof(size_t)             // space for info about allocated capacity
                      + sizeof(size_t);            // space for info about stored element count
    
    (*capa)  = new_capacity;
    (*count) = min(*count, new_capacity);

    void* src_ptr = (array-3*sizeof(size_t));
    src_ptr = realloc(src_ptr, size_total);
    return src_ptr+3*sizeof(size_t);
}


inline 
long array_count(void* array)
{
    array -= 1*sizeof(size_t);
    return *(size_t*)(array);
}


inline 
long array_capa(void* array)
{
    array -= 2*sizeof(size_t);
    return *(size_t*)(array);
}


inline 
size_t array_esize(void* array)
{
    array -= 3*sizeof(size_t);
    return *(size_t*)(array);
}


void array_push(void** array_ptr, const void* data)
{
    void* array = *array_ptr;

    size_t* count = (size_t*)(array-1*sizeof(size_t));
    size_t* capa  = (size_t*)(array-2*sizeof(size_t));
    size_t* esize = (size_t*)(array-3*sizeof(size_t));

    if (*count+1 > *capa) {
        array = *array_ptr = array_realloc(array, 2*(*capa));
        count = (size_t*)(array-1*sizeof(size_t));
        capa  = (size_t*)(array-2*sizeof(size_t));
        esize = (size_t*)(array-3*sizeof(size_t));
    }

    void* dptr = array + (*esize)*(*count);
    memcpy(dptr, data, *esize);
    (*count)++;
}


inline
void array_push_int(void** array_ptr, const int data)
{
    array_push(array_ptr, &data);
}


inline
void array_push_long(void** array_ptr, const long data)
{
    array_push(array_ptr, &data);
}


inline
void array_push_double(void** array_ptr, const double data)
{
    array_push(array_ptr, &data);
}


int  array_exists(void* array, const void* data)
{
    size_t count = *(size_t*)(array-1*sizeof(size_t));
    size_t esize = *(size_t*)(array-3*sizeof(size_t));
    
    size_t i;
    for (i=0; i<count; i++) {
        if (memcmp(array, data, esize) == 0) return 1;
        array += esize;
    }

    return 0;
}

inline
void array_push_if_not_exists(void** array_ptr, const void* data)
{
    if (!array_exists(*array_ptr, data)) array_push(array_ptr, data);
}


void array_reverse(void* array)
{
    long count = *(size_t*)(array-1*sizeof(size_t));
    long esize = *(size_t*)(array-3*sizeof(size_t));
    
    void* frst_ptr = array;
    void* last_ptr = array + esize*(count-1);
    void* tmpdata1 = malloc(esize);
    void* tmpdata2 = malloc(esize);
    
    while(frst_ptr < last_ptr) {
        memcpy(tmpdata1, frst_ptr, esize);
        memcpy(tmpdata2, last_ptr, esize);
        memcpy(frst_ptr, tmpdata2, esize);
        memcpy(last_ptr, tmpdata1, esize);
        frst_ptr += esize;
        last_ptr -= esize;
    }
    
    free(tmpdata1);
    free(tmpdata2);
}


