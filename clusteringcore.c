/*
 * Copyright 2014 Uri Laserson
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// clusteringcore.c
// Defines clusteringcore python extension module

#include <Python.h>
#include <numpy/arrayobject.h>

// DEBUG
//#include <stdio.h>

int intarraymin( int *data, int n ) {
    int i, minval ;
    for ( minval = data[0], i = 1 ; i < n ; i++ ) {
        if ( data[i] < minval ) {
            minval = data[i] ;
        }
    }
    return minval ;
}

static PyObject *clusteringcore_levenshtein( PyObject *self, PyObject *args ){
    char *seq1, *seq2 ;
    int len1, len2 ;
    npy_intp dim[2] ;
    int i, j ;
    int cost, best ;
    int ext[3] ;
    PyArrayObject *scores = NULL;
    
    // get sequence args
    if ( !PyArg_ParseTuple(args,"s#s#",
                            &seq1, &len1,
                            &seq2, &len2) ) {
        return NULL ;
    }
    
    // check for trivial case
    if ( len1 == 0 || len2 == 0 ) {
        return Py_BuildValue( "i", (len1 < len2 ? len2 : len1) ) ;
    }
    
    // allocate and initialize score matrix F
    dim[0] = len1+1 ;
    dim[1] = len2+1 ;
    scores = (PyArrayObject *)PyArray_ZEROS( 2, dim, NPY_INT, 0 ) ;
    if (scores == NULL) return NULL ;
    
    for ( i = 0, j = 0 ; i <= len1 ; i++ ) {
        *((int*)PyArray_GETPTR2(scores,i,j)) = i ;
    }
    for ( i = 0, j = 0 ; j <= len2 ; j++ ) {
        *((int*)PyArray_GETPTR2(scores,i,j)) = j ;
    }
    
    // compute DP score matrix
    for ( i = 1 ; i <= len1 ; i++ ) {
        for ( j = 1 ; j <= len2 ; j++ ) {
            cost = (seq1[i-1] == seq2[j-1]) ? 0 : 1 ;
            
            ext[0] = *((int*)PyArray_GETPTR2(scores,i-1,j-1)) + cost ;
            ext[1] = *((int*)PyArray_GETPTR2(scores,i-1,j))   + 1 ;
            ext[2] = *((int*)PyArray_GETPTR2(scores,i,j-1))   + 1 ;
            
            *((int*)PyArray_GETPTR2(scores,i,j)) = intarraymin(ext,3) ;
        }
    }
    
    best = *((int*)PyArray_GETPTR2(scores,len1,len2)) ;
    Py_DECREF(scores) ;
    return Py_BuildValue( "i", best ) ;
}

static PyMethodDef clusteringcoremethods[] = {
    {"levenshtein", clusteringcore_levenshtein, METH_VARARGS}
} ;

PyMODINIT_FUNC initclusteringcore() {
    Py_InitModule( "clusteringcore", clusteringcoremethods ) ;
    import_array() ;
}
