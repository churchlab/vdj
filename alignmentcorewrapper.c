// alignmentcorewrapper.c
// Defines alignmentcore python extension module

#include <Python.h>
#include <numpy/arrayobject.h>

// alignmentcore.alignNW(M,Ix,Iy,BT,seq1,seq2)

// extern declaration
extern void alignNW( PyArrayObject *M, PyArrayObject *Ix, PyArrayObject *Iy, PyArrayObject *BT,
					 const char *seq1, const int len1, const char *seq2, const int len2 ) ;

extern void alignSW( PyArrayObject *F, PyArrayObject *BT,
					 const char *seq1, const int len1, const char *seq2, const int len2 ) ;

static PyObject *alignmentcore_alignNW( PyObject *self, PyObject *args ) {
	const char *seq1, *seq2 ;
	const int len1, len2 ;
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

static PyObject *alignmentcore_alignSW( PyObject *self, PyObject *args ) {
	const char *seq1, *seq2 ;
	const int len1, len2 ;
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
}