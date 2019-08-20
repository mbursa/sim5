//! \file sim5utils.c
//! Utility routines.



void gprintf(FILE* file, const char *templatex, ...)
//**************************************************
// prints message both to file and stderr
{
    #ifndef CUDA
    va_list ap;

    // print to stderr
    va_start (ap, templatex);
    vfprintf (stderr, templatex, ap);
    va_end (ap);

    // print to file
    va_start (ap, templatex);
    vfprintf (file, templatex, ap);
    va_end (ap);
    #endif
}


void warning(const char *templatex, ...)
//**************************************************
// prints warning to stderr
{
    #ifndef CUDA
    va_list ap;
    fprintf (stderr, "WRN: ");
    va_start (ap, templatex);
    vfprintf (stderr, templatex, ap);
    va_end (ap);
	fprintf(stderr, "\n");
    #endif
}


void error(const char *templatex, ...)
//**************************************************
// prints error to stderr
{
    #ifndef CUDA
    va_list ap;
    fprintf (stderr, "ERROR: ");
    va_start (ap, templatex);
    vfprintf (stderr, templatex, ap);
    va_end (ap);
    fprintf(stderr, "\n");
    //exit(-1);
    #endif
}




int sort_array_d_compare_func(const void *x, const void *y) {
    return (*(double*)x - *(double*)y);
}

void sort_array(double *array, int N)
//**************************************************
{
    qsort(array, N, sizeof(double), sort_array_d_compare_func);
}



int sort_array_f_compare_func(const void *x, const void *y) {
    return (*(float*)x - *(float*)y);
}

void sort_array_f(float *array, int N)
//**************************************************
{
    qsort(array, N, sizeof(float), sort_array_f_compare_func);
}


#ifndef CUDA
void* sim5_array_alloc(size_t capacity, size_t element_size)
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


void sim5_array_free(void* ptr)
{
    ptr -= 3*sizeof(size_t);
    free(ptr);
}


void* sim5_array_realloc(void* array, size_t new_capacity)
{
    size_t* count = (size_t*)(array-1*sizeof(size_t));
    size_t* capa  = (size_t*)(array-2*sizeof(size_t));
    size_t* esize = (size_t*)(array-3*sizeof(size_t));
    
    size_t size_total = new_capacity * (*esize)    // space allocated for elements
                      + sizeof(size_t)             // space for info about size of stored elements 
                      + sizeof(size_t)             // space for info about allocated capacity
                      + sizeof(size_t);            // space for info about stored element count
    
    (*capa)  = new_capacity;
    (*count) = fmin(*count, new_capacity);

    void* src_ptr = (array-3*sizeof(size_t));
    src_ptr = realloc(src_ptr, size_total);
    return src_ptr+3*sizeof(size_t);
}


inline 
long sim5_array_count(void* array)
{
    array -= 1*sizeof(size_t);
    return *(size_t*)(array);
}


inline 
long sim5_array_capa(void* array)
{
    array -= 2*sizeof(size_t);
    return *(size_t*)(array);
}


inline 
size_t sim5_array_esize(void* array)
{
    array -= 3*sizeof(size_t);
    return *(size_t*)(array);
}


void sim5_array_push(void** array_ptr, const void* data)
{
    void* array = *array_ptr;

    size_t* count = (size_t*)(array-1*sizeof(size_t));
    size_t* capa  = (size_t*)(array-2*sizeof(size_t));
    size_t* esize = (size_t*)(array-3*sizeof(size_t));

    if (*count+1 > *capa) {
        array = *array_ptr = sim5_array_realloc(array, 2*(*capa));
        count = (size_t*)(array-1*sizeof(size_t));
        capa  = (size_t*)(array-2*sizeof(size_t));
        esize = (size_t*)(array-3*sizeof(size_t));
    }

    void* dptr = array + (*esize)*(*count);
    memcpy(dptr, data, *esize);
    (*count)++;
}


inline
void sim5_array_push_int(void** array_ptr, const int data)
{
    sim5_array_push(array_ptr, &data);
}


inline
void sim5_array_push_long(void** array_ptr, const long data)
{
    sim5_array_push(array_ptr, &data);
}


inline
void sim5_array_push_float(void** array_ptr, const float data)
{
    sim5_array_push(array_ptr, &data);
}


inline
void sim5_array_push_double(void** array_ptr, const double data)
{
    sim5_array_push(array_ptr, &data);
}


int  sim5_array_exists(void* array, const void* data)
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
void sim5_array_push_if_not_exists(void** array_ptr, const void* data)
{
    if (!sim5_array_exists(*array_ptr, data)) sim5_array_push(array_ptr, data);
}


void sim5_array_reverse(void* array)
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
#endif


char* key_value_get(const char *string, const char *key) 
{
    static char value[256];
    char* p;
    
    if ((!string) || (!key)) return NULL;    

    // make local copy of the string (segfault comes othervise)
    char* string_copy = (char*)malloc(strlen(string)+1);
    strcpy (string_copy, string);

    for (p = strtok(string_copy,","); p != NULL; p = strtok(NULL, ",")) {
        if ((strncmp(p, key, strlen(key))==0) && (p[strlen(key)]=='=')) {
            p += strlen(key) + 1;
            strncpy(value, p, sizeof(value));
            free(string_copy);
            return value;
        }
    }

    // if no match    
    free(string_copy);
    return NULL;
}




#define MAXTOKENS       1024
#define MAXLINE         8096     // fgets buff
#define MINLEN          3        // skip lines shorter as

// split string into tokens, return token array
char **split(char *string, char *delim) {
    char **tokens = NULL;
    char *working = NULL;
    char *token = NULL;
    int idx = 0;

    tokens  = (char**)malloc(sizeof(char *) * MAXTOKENS);
    if(tokens == NULL) return NULL;
    working = (char*)malloc(sizeof(char) * strlen(string) + 1);
    if(working == NULL) return NULL;

    // to make sure, copy string to a safe place
    strcpy(working, string);
    for(idx = 0; idx < MAXTOKENS; idx++) tokens[idx] = NULL;

    token = strtok(working, delim);
    idx = 0;

    // always keep the last entry NULL termindated
    while((idx < (MAXTOKENS - 1)) && (token != NULL)) {
        tokens[idx] = (char*)malloc(sizeof(char) * strlen(token) + 1);
        if(tokens[idx] != NULL) {
            strcpy(tokens[idx], token);
            idx++;
            token = strtok(NULL, delim);
        }
    }

    free(working);
    return tokens;
}



long getlinecount(FILE* f) {
    long fpos = ftell(f);
    char line[8192];
    long result = 0;
    fseek(f,0,SEEK_SET);
    while (fgets(line, 8192, f) != NULL) result++; 
    fseek(f,fpos,SEEK_SET);
    return result;
}






// Call this function to get a backtrace.
void backtrace() {
/*
    unw_cursor_t cursor;
    unw_context_t context;

    unw_getcontext(&context);
    unw_init_local(&cursor, &context);

    int n=0;
    while ( unw_step(&cursor) ) {
        unw_word_t ip, sp, off;

        unw_get_reg(&cursor, UNW_REG_IP, &ip);
        unw_get_reg(&cursor, UNW_REG_SP, &sp);

        char symbol[256] = {"<unknown>"};
        char *name = symbol;

        off = 0;
        //if ( !unw_get_proc_name(&cursor, symbol, sizeof(symbol), &off) ) {
        //  int status;
        //  if ( (name = abi::__cxa_demangle(symbol, NULL, NULL, &status)) == 0 )
        //    name = symbol;
        //}

        printf("#%-2d 0x%016" PRIxPTR " sp=0x%016" PRIxPTR " %s + 0x%" PRIxPTR "\n", ++n, ip, sp, name, off);

        if (name != symbol) free(name);
    }
*/
}

    
