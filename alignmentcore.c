#include "kbandalign.h"
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
			
			*BTij = best ;	// 0 = (i-1,j-1) ; 1 = (i-1,j) ; 2 = (i,j-1)
			
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
			*BTij = best ;	// 0 = (i-1,j-1) ; 1 = (i-1,j) ; 2 = (i,j-1) ; 3 = END (0)
			
			// DEBUG
			//printf("char1: %c\tchar2: %c\tmatch? %d\ts: %f\tF: %f\tBT: %i\n",seq1[i-1],seq2[j-1],seq1[i-1]==seq2[j-1],s,*Fij,*BTij) ;
		}
	}
	return ;
}

static PyObject *alignmentcore_kalign( PyObject *self, PyObject *args ){
    
    char *seq1, *seq2;

    if( !PyArg_ParseTuple(args, "ss",
                            &seq1,
                            &seq2) ){
        return NULL;
    }

    KAlignment *algn = new KAlignment(seq1,seq2);
    algn->align();

    int similarity = algn->globalScore();
    delete algn;

    return Py_BuildValue( "i", similarity ); 

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




















// ALTERNATE IMPLEMENTATION of alignNW
/*
void alignNW( PyArrayObject *M, PyArrayObject *Ix, PyArrayObject *Iy, PyArrayObject *BT,
			  const char *seq1, const int len1, const char *seq2, const int len2) {
	
	// extract proper pointers etc.
	char *Mdata  = M->data  ;
	char *Ixdata = Ix->data ;
	char *Iydata = Iy->data ;
	char *BTdata = BT->data ;
	
	npy_intp *Mstrides  = M->strides  ;
	npy_intp *Ixstrides = Ix->strides ;
	npy_intp *Iystrides = Iy->strides ;
	npy_intp *BTstrides = BT->strides ;
	
	// define parameters
	double match     =  0.5  ;
	double mismatch  = -0.75 ;
	double gapopen   = -2.0  ;
	double gapextend = -1.5  ;
	
	int nrows = len1 + 1 ;
	int ncols = len2 + 1 ;
	int i, j ;
	int im1, jm1 ;
	
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
			
			
			//Mij  = (double*)PyArray_GETPTR2(M, i,j) ;
			//#define PyArray_BYTES(obj) (((PyArrayObject *)(obj))->data)
			//#define PyArray_STRIDES(obj) (((PyArrayObject *)(obj))->strides)
			//#define PyArray_GETPTR2(obj, i, j) ((void *)(PyArray_BYTES(obj) +             \
			//                                            (i)*PyArray_STRIDES(obj)[0] +     \
			//                                            (j)*PyArray_STRIDES(obj)[1]))
			//
			//Mim1jm1  = (double*)PyArray_GETPTR2(M, i-1,j-1) ;
			
			
			
			im1 = i-1 ;
			jm1 = j-1 ;
			
			Mij  = (double*)(Mdata  + i*Mstrides[0]  + j*Mstrides[1])  ;
			Ixij = (double*)(Ixdata + i*Ixstrides[0] + j*Ixstrides[1]) ;
			Iyij = (double*)(Iydata + i*Iystrides[0] + j*Iystrides[1]) ;
			Mim1jm1  = (double*)(Mdata  + im1*Mstrides[0]  + jm1*Mstrides[1])  ;
			Ixim1jm1 = (double*)(Ixdata + im1*Ixstrides[0] + jm1*Ixstrides[1]) ;
			Iyim1jm1 = (double*)(Iydata + im1*Iystrides[0] + jm1*Iystrides[1]) ;
			Mim1j = (double*)(Mdata  + im1*Mstrides[0]  + j*Mstrides[1]) ;
			Mijm1 = (double*)(Mdata  + i*Mstrides[0]  + jm1*Mstrides[1]) ;
			Ixim1j = (double*)(Ixdata + im1*Ixstrides[0] + j*Ixstrides[1]) ;
			Iyijm1 = (double*)(Iydata + i*Iystrides[0] + jm1*Iystrides[1]) ;
			BTij = (int*)(BTdata  + i*BTstrides[0]  + j*BTstrides[1]) ;
			
			s = (seq1[i-1] == seq2[j-1]) ? match : mismatch ;
			ext[0] = *((double*)(Mdata  + im1*Mstrides[0]  + jm1*Mstrides[1]))  + s ;
			ext[1] = *((double*)(Ixdata + im1*Ixstrides[0] + jm1*Ixstrides[1])) + s ;
			ext[2] = *((double*)(Iydata + im1*Iystrides[0] + jm1*Iystrides[1])) + s ;
			IxGO = *Mim1j  + gapopen   ;
			IxGE = *Ixim1j + gapextend ;
			IyGO = *Mijm1  + gapopen   ;
			IyGE = *Iyijm1 + gapextend ;
			
			best = arrayargmax(ext,3) ;
			
			*Mij = ext[best] ;
			*Ixij = (IxGO >= IxGE) ? IxGO : IxGE ;
			*Iyij = (IyGO >= IyGE) ? IyGO : IyGE ;
			
			*BTij = best ;	// 0 = (i-1,j-1) ; 1 = (i-1,j) ; 2 = (i,j-1)
			
			// DEBUG
			//printf("char1: %c\tchar2: %c\tmatch? %d\ts: %f\tM: %f\tIx: %f\tBT: %i\n",seq1[i-1],seq2[j-1],seq1[i-1]==seq2[j-1],s,*Mij,*Ixij,*BTij) ;
		}
	}
	return ;
}
*/
