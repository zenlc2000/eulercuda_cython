
cimport cython
cimport pygpuhash
#from pygpuhash cimport createHashTable

#	createHashTable(d_kmerKeys, d_kmerValues, kmerCount, &d_TK, &d_TV,
#			&tableLength, &d_bucketSeed, &bucketCount);


def create_hash_table(
        KEY_T[::1] d_keys,
        VALUE_T[::1] d_values,
        unsigned int length):

    cdef KEY_T *d_TK
    cdef VALUE_T *d_TV
    cdef unsigned int tableLength
    cdef unsigned int *d_bucketSize
    cdef unsigned int  bucketCount

    return pygpuhash.createHashTable(
        &d_keys[0],
        &d_values[0],
        length,
        &d_TK,
        &d_TV,
        &tableLength,
        &d_bucketSize,
        &bucketCount)