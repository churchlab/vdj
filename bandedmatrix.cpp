#include <stdio.h>

#include "bandedmatrix.h"

using namespace std;

template <class S>
MatrixRow<S>::MatrixRow(int offset, int length, S* row){
    _offset = offset;
    _length = length;
    _dat = row;
}

template <class S>
S& MatrixRow<S>::operator[](const int ii){
    if( ii - _offset >= 0 && ii - _offset < _length){
        return _dat[ii - _offset]; 
    } else {
        raise(SIGSEGV);
    }
}

template <class T>
BandedMatrix<T>::BandedMatrix(int x, int y){
    _numCols = x;
    _numRows = y;
    _left  = (_numCols < _numRows ? _numRows - _numCols : 1);
    _right = (_numRows < _numCols ? _numCols - _numRows : 1);
    
    _cols = vector<T*>(_numCols);
    _startPos = vector<int>(_numCols);
    resize(false, false);
    _dat = (T*) malloc(sizeof(T) * _size);
    makeCols();
}

template <class T>
MatrixRow<T> BandedMatrix<T>::operator[](const int ii){
    int offset = _startPos[ii];
    int length = (int) (ii == _numRows - 1 ? (_dat + _size) - _cols[ii] : _cols[ii + 1] - _cols[ii]);
    
    if( ii >= 0 && ii < _numCols ){
        return MatrixRow<T>(offset, length, _cols[ii]);
    } else {
        raise(SIGSEGV);
    }
}

template <class T>
void BandedMatrix<T>::makeCols(){

    int colLength = 1 + _left;
    T* index = _dat; 
    int effectiveRight = 0;
    int effectiveLeft  = _left;
    int colStart = 0;
   for(int ii = 0; ii < _numCols; ++ii ){
        _cols[ii] = index;
        _startPos[ii] = colStart;
        
        if( ii + _left + 1 >= _numRows && effectiveLeft >= 0 ){
            effectiveLeft--;
        }

        if( ii < _right ){
            effectiveRight++;
        } else {
            colStart++;
        }

        colLength = (ii < _numRows ? 1 : 0) 
            + effectiveRight 
            + effectiveLeft;
        index += colLength;
   }
}

template <class T>
int BandedMatrix<T>::resize(bool gLeft, bool gRight){
    
    if( gLeft ){
        _left = min(_left * 2, _numRows - 1);
    }

    if( gRight ){
        _right = min(_right * 2, _numCols - 1);
    }

    /*
     * Lop off the symmetric component of the matrix and
     * compute the size of the area to the left of the diagonal
     * and to the right of it.
     */
    int upperTriangle = _numCols - 1 - _right;
    int lowerTriangle = _numRows - 1 - _left;

    _size = 
        (_numCols * _numRows) -
        (upperTriangle * (upperTriangle + 1)) / 2 -
        (lowerTriangle * (lowerTriangle + 1)) / 2;
    assert(_size > 0);   
    return _size;
}

template <class T>
T BandedMatrix<T>::set(int x, int y, T v){
    printf("Setting (%d, %d) at ", x, y);
    T* ptr = getPtr(x,y);
    printf("%p\n", ptr);
    *ptr = v;
    return *ptr;
}

template <class T>
T BandedMatrix<T>::get(int x, int y){
    printf("Getting (%d, %d) at ", x, y);
    T* ptr = (getPtr(x,y));
    printf("%p\n", ptr);
    return *ptr;
}

template <class T>
T* BandedMatrix<T>::getPtr(int x, int y){
    T* col_ptr = _cols[x];
    int offset = _startPos[x];
    int index  = x - offset;

    if( index < 0 ){
        return NULL;    
    }
    
    if( index > _size ){
        return NULL;
    }

    if( x > _numCols || x < 0 ){
        return NULL;
    }

    if( y > _numRows || y < 0 ){
        return NULL;
    }

    return (col_ptr + index);

}

int main(int argc, void** argv){

    BandedMatrix<int> m = BandedMatrix<int>(5,5);
    printf("Set %d\n", m.set(4,4,1));
    printf("Got %d\n", m.get(4,4));
    printf("%d\n", m[4][4]);
    m[4][4] = -1;
    printf("%d\n", m[4][4]);
    return 0; 
}
