#ifndef _SIM5UTILS_H_
#define _SIM5UTILS_H_

#ifdef __cplusplus
extern "C" {
#endif

void gprintf(FILE* file, const char *templatex, ...);

void error(const char *templatex, ...);
void warning(const char *templatex, ...);

void sort_array(double *array, int N);
void sort_array_f(float *array, int N);
// sorts array of numbers

#ifndef CUDA
void* sim5_array_alloc(size_t capacity, size_t element_size);
void sim5_array_free(void* array);
void* sim5_array_realloc(void* array, size_t new_capacity);
long sim5_array_count(void* arrry);
long sim5_array_capa(void* array);
size_t sim5_array_esize(void* array);
void sim5_array_push(void** array_ptr, const void* data);
void sim5_array_push_int(void** array_ptr, const int data);
void sim5_array_push_long(void** array_ptr, const long data);
void sim5_array_push_double(void** array_ptr, const double data);
void sim5_array_push_float(void** array_ptr, const float data);
int  sim5_array_exists(void* array, const void* data);
void sim5_array_push_if_not_exists(void** array_ptr, const void* data);
void sim5_array_reverse(void* array);
#endif

char* key_value_get(const char *string, const char *key);

void backtrace();

#ifdef __cplusplus
}
#endif


#endif 
