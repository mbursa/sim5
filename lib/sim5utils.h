#ifndef _SIM5UTILS_H_
#define _SIM5UTILS_H_

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

void gprintf(FILE* file, const char *template, ...);

void error(const char *template, ...);

void warning(const char *template, ...);

void sort_array(double *array, int N);
void sort_array_f(float *array, int N);
// sorts array of numbers


void* array_alloc(size_t capacity, size_t element_size);
void array_free(void* array);
void* array_realloc(void* array, size_t new_capacity);
long array_count(void* arrry);
long array_capa(void* array);
size_t array_esize(void* array);
void array_push(void** array_ptr, const void* data);
void array_push_int(void** array_ptr, const int data);
void array_push_long(void** array_ptr, const long data);
void array_push_double(void** array_ptr, const double data);
int  array_exists(void* array, const void* data);
void array_push_if_not_exists(void** array_ptr, const void* data);
void array_reverse(void* array);

char* key_value_get(const char *string, const char *key);

#endif 
