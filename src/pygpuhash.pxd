
cdef extern from "common.h":
    ctypedef unsigned long long int KEY_T
    ctypedef unsigned long long int * KEY_PTR
    ctypedef unsigned long long int VALUE_T
    ctypedef unsigned long long int * VALUE_PTR

cdef extern from "gpuhash.h":
    void createHashTable(
            KEY_PTR d_keys,
            VALUE_PTR d_values,
            unsigned int length,
            KEY_PTR * d_TK,
            VALUE_PTR * d_TV,
            unsigned int * tableLength,
            unsigned int ** d_bucketSize,
            unsigned int * bucketCount);