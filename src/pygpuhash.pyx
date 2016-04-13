

from pygpuhash cimport createHashTable


def create_hash_table(KEY_T[::1] d_keys,
                      VALUE_T[::1] d_values,
                      unsigned int length,
                      KEY_T[::1] * d_TK,
                      VALUE_T[::1] * d_TV,
                      unsigned int[::1] tableLength,
                      unsigned int[::1] * d_bucketSize,
                      unsigned int[::1] bucketCount):
    return createHashTable(&d_keys[0], &d_values[0], length, &d_TK[0],
                           &d_TV[0], &tableLength[0], d_bucketSize[0],&bucketCount[0])