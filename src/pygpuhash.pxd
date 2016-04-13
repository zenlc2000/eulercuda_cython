
cdef extern from "common.h":
	ctypedef int KEY_T
	ctypedef int * KEY_PTR
	ctypedef int VALUE_T
	ctypedef int * VALUE_PTR

cdef extern from "gpuhash.h":
    void createHashTable(KEY_PTR d_keys,
                VALUE_PTR d_values,
                unsigned int length,
                KEY_PTR * d_TK,
                VALUE_PTR * d_TV,
                unsigned int * tableLength,
                unsigned int ** d_bucketSize,
                unsigned int * bucketCount);