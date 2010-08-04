#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdlib.h>

// DEBUG
//#include <stdio.h>

int arrayargmax( double *data, int n ) {
    int maxidx, curidx ;
    for ( maxidx = 0, curidx = 1 ; curidx < n ; curidx++ ) {
        if ( data[curidx] > data[maxidx] ) {
            maxidx = curidx ;
        }
    }
    return maxidx ;
}

void alignNW( PyArrayObject *M, PyArrayObject *Ix, PyArrayObject *Iy, PyArrayObject *BT,
              const char *seq1, const int len1, const char *seq2, const int len2) {
    
    // define parameters
    double match     =  0.5  ;
    double mismatch  = -0.75 ;
    double gapopen   = -2.0  ;
    double gapextend = -1.5  ;
    
    int nrows = len1 + 1 ;
    int ncols = len2 + 1 ;
    int i, j ;
    
    double *Mij, *Mim1jm1, *Mim1j, *Mijm1 ;
    double *Ixij, *Ixim1jm1, *Ixim1j ;
    double *Iyij, *Iyim1jm1, *Iyijm1 ;
    int *BTij ;
    
    double s ;
    double ext[3] ;
    double IxGO, IxGE, IyGO, IyGE ;
    int best ;
    
    for ( i = 1 ; i < nrows ; i++ ) {
        for ( j = 1 ; j < ncols ; j++ ) {
            
            Mij  = (double*)PyArray_GETPTR2(M, i,j) ;
            Ixij = (double*)PyArray_GETPTR2(Ix,i,j) ;
            Iyij = (double*)PyArray_GETPTR2(Iy,i,j) ;
            Mim1jm1  = (double*)PyArray_GETPTR2(M, i-1,j-1) ;
            Ixim1jm1 = (double*)PyArray_GETPTR2(Ix,i-1,j-1) ;
            Iyim1jm1 = (double*)PyArray_GETPTR2(Iy,i-1,j-1) ;
            Mim1j = (double*)PyArray_GETPTR2(M,i-1,j) ;
            Mijm1 = (double*)PyArray_GETPTR2(M,i,j-1) ;
            Ixim1j = (double*)PyArray_GETPTR2(Ix,i-1,j) ;
            Iyijm1 = (double*)PyArray_GETPTR2(Iy,i,j-1) ;
            BTij = (int*)PyArray_GETPTR2(BT,i,j) ;
            
            s = (seq1[i-1] == seq2[j-1]) ? match : mismatch ;
            ext[0] = *Mim1jm1  + s ;
            ext[1] = *Ixim1jm1 + s ;
            ext[2] = *Iyim1jm1 + s ;
            IxGO = *Mim1j  + gapopen   ;
            IxGE = *Ixim1j + gapextend ;
            IyGO = *Mijm1  + gapopen   ;
            IyGE = *Iyijm1 + gapextend ;
            
            best = arrayargmax(ext,3) ;
            
            *Mij = ext[best] ;
            *Ixij = (IxGO >= IxGE) ? IxGO : IxGE ;
            *Iyij = (IyGO >= IyGE) ? IyGO : IyGE ;
            
            *BTij = best ;  // 0 = (i-1,j-1) ; 1 = (i-1,j) ; 2 = (i,j-1)
            
            // DEBUG
            //printf("char1: %c\tchar2: %c\tmatch? %d\ts: %f\tM: %f\tIx: %f\tBT: %i\n",seq1[i-1],seq2[j-1],seq1[i-1]==seq2[j-1],s,*Mij,*Ixij,*BTij) ;
        }
    }
    return ;
}

void alignSW( PyArrayObject *F, PyArrayObject *BT,
              const char *seq1, const int len1, const char *seq2, const int len2) {
    
    // define parameters
    double match     =  0.5  ;
    double mismatch  = -0.75 ;
    double gapextend = -1.5  ;
    
    int nrows = len1 + 1 ;
    int ncols = len2 + 1 ;
    int i, j ;
    
    double *Fij, *Fim1jm1, *Fim1j, *Fijm1 ;
    int *BTij ;
    
    double s ;
    double ext[4] ;
    int best ;
    
    for ( i = 1 ; i < nrows ; i++ ) {
        for ( j = 1 ; j < ncols ; j++ ) {
            
            Fij     = (double*)PyArray_GETPTR2(F,i,j) ;
            Fim1jm1 = (double*)PyArray_GETPTR2(F,i-1,j-1) ;
            Fim1j   = (double*)PyArray_GETPTR2(F,i-1,j) ;
            Fijm1   = (double*)PyArray_GETPTR2(F,i,j-1) ;
            BTij    = (int*)PyArray_GETPTR2(BT,i,j) ;
            
            s = (seq1[i-1] == seq2[j-1]) ? match : mismatch ;
            ext[0] = *Fim1jm1  + s ;
            ext[1] = *Fim1j + gapextend ;
            ext[2] = *Fijm1 + gapextend ;
            ext[3] = 0 ;
            
            best = arrayargmax(ext,4) ;
            
            *Fij = ext[best] ;
            *BTij = best ;  // 0 = (i-1,j-1) ; 1 = (i-1,j) ; 2 = (i,j-1) ; 3 = END (0)
            
            // DEBUG
            //printf("char1: %c\tchar2: %c\tmatch? %d\ts: %f\tF: %f\tBT: %i\n",seq1[i-1],seq2[j-1],seq1[i-1]==seq2[j-1],s,*Fij,*BTij) ;
        }
    }
    return ;
}

static PyObject *alignmentcore_alignNW( PyObject *self, PyObject *args ) {
    char *seq1, *seq2 ;
    int len1, len2 ;
    PyArrayObject *M, *Ix, *Iy, *BT ;
    
    // get sequence args
    if ( !PyArg_ParseTuple(args,"OOOOs#s#",
                            &M,
                            &Ix,
                            &Iy,
                            &BT,
                            &seq1, &len1,
                            &seq2, &len2) ) {
        return NULL ;
    }
    
    // call function
    alignNW( M, Ix, Iy, BT, seq1, len1, seq2, len2 ) ;
    
    return Py_BuildValue( "d", 0.0 ) ;
}

// wrapper functions

static PyObject *alignmentcore_alignSW( PyObject *self, PyObject *args ) {
    char *seq1, *seq2 ;
    int len1, len2 ;
    PyArrayObject *F, *BT ;
    
    // get sequence args
    if ( !PyArg_ParseTuple(args,"OOs#s#",
                            &F,
                            &BT,
                            &seq1, &len1,
                            &seq2, &len2) ) {
        return NULL ;
    }
    
    // call function
    alignSW( F, BT, seq1, len1, seq2, len2 ) ;
    
    return Py_BuildValue( "d", 0.0 ) ;
}

static PyMethodDef alignmentcoremethods[] = {
    {"alignNW", alignmentcore_alignNW, METH_VARARGS},
    {"alignSW", alignmentcore_alignSW, METH_VARARGS}
} ;

void initalignmentcore() {
    Py_InitModule( "alignmentcore", alignmentcoremethods ) ;
    import_array();
}
