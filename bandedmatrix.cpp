#include <stdio.h>

#include "bandedmatrix.h"

using namespace std;

/*
 * BandedMatrixIterator
 */

template <class T>
BandedMatrixIterator<T>::BandedMatrixIterator(){
    _size = 0;
    _index = 0;
    _x = 0;
    _y = 0;
}

template <class T>
BandedMatrixIterator<T>::BandedMatrixIterator(const BandedMatrixIterator<T> &itr){
    _cols = itr._cols;
    _currCol = itr._currCol;
    _size = itr._size;
    _index = itr._index;
    _x = itr._x;
    _y = itr._y;
    _dat = itr._dat;
}

template <class T>
BandedMatrixIterator<T>::BandedMatrixIterator(vector<MatrixCol<T> > &cols, int size){
    _currCol = cols[0];
    _dat = _currCol._dat;
    _cols = cols;
    _x = 0;
    _y = 0;
    _size = size;

    
}

template <class T>
BandedMatrixIterator<T>::~BandedMatrixIterator(){

}

template <class T>
bool BandedMatrixIterator<T>::operator==(BandedMatrixIterator<T> &itr){

    return( _dat == itr._dat && _index == itr._index );

}

template <class T>
bool BandedMatrixIterator<T>::operator!=(BandedMatrixIterator<T> &itr){
    return !( _dat == itr._dat || _index != itr._index );
}

template <class T>
bool BandedMatrixIterator<T>::hasNext(){
    return _index < _size - 1;
}

template <class T>
bool BandedMatrixIterator<T>::hasPrev(){
    return _index > 0;
}

template <class T>
BandedMatrixIterator<T>& BandedMatrixIterator<T>::operator++() {
    _y++;
    _index++;
    if( _y - _currCol._offset - 1 > _currCol._length ){
        _x++;
        _currCol = _cols[_x];
        _y = _currCol._offset;
    }
    return *this;
}

template <class T>
BandedMatrixIterator<T>& BandedMatrixIterator<T>::operator--() {
    _y--;
    _index--;
    if( _y - _currCol._offset < 0 ){
        _x--;
        _currCol = _cols[_x];
        _y = _currCol._offset;
    }
    return *this;
}


template <class T>
T& BandedMatrixIterator<T>::operator*() const {
    printf("%p + %d = %p\n", _dat, _index, _dat + _index);
    return _dat[_index];
}

template <class T>
pair<int, int> BandedMatrixIterator<T>::coordinates(){
    return pair<int, int>(_x,_y);
}

/*
 * MatrixCol: Convenient addressing of matrix rows
 */
template <class S>
MatrixCol<S>::MatrixCol(){
    _offset = 0;
    _length = 0;
    _dat = NULL;
}

template <class S>
MatrixCol<S>::MatrixCol(int offset, int length, S* row){
    _offset = offset;
    _length = length;
    _dat = row;
}

template <class S>
S& MatrixCol<S>::operator[](const int ii){
    if( ii - _offset >= 0 && ii - _offset < _length){
        return _dat[ii - _offset]; 
    } else {
        raise(SIGSEGV);
    }
}

template <class T>
BandedMatrix<T>::BandedMatrix(){
    _numCols = 0;
    _numRows = 0;
    _left    = 0;
    _right   = 0;
    _cols    = vector<MatrixCol<T> >(0);
    _startPos = vector<int>(0);
    _size    = 0;
}

template <class T>
BandedMatrix<T>::BandedMatrix(int x, int y){
    _numCols = x;
    _numRows = y;
    _left  = (_numCols < _numRows ? _numRows - _numCols : 1);
    _right = (_numRows < _numCols ? _numCols - _numRows : 1);
    
    _cols = vector<MatrixCol<T> >(_numCols);
    _startPos = vector<int>(_numCols);
    resize(false, false);
    _dat = (T*) malloc(sizeof(T) * _size);
    makeCols();
    _begin = BandedMatrixIterator<T>(_cols, _size);
    _end   = BandedMatrixIterator<T>(_cols, _size);
}

template <class T>
MatrixCol<T> BandedMatrix<T>::operator[](const int ii){
    return _cols[ii];
}

template<class T>
bool BandedMatrix<T>::inbounds(int x, int y){
    if( x < 0 || x >= _numCols || y < 0 || y >= _numRows ){
        return false;
    }

    MatrixCol<T> column = _cols[x];
    
    if( y - column._offset < 0 || y - column._offset >= column._length ){
        return false;
    }

    return true;

}

template <class T>
void BandedMatrix<T>::makeCols(){
    int colLength = 1 + _left;
    T* index = _dat; 
    int effectiveRight = 0;
    int effectiveLeft  = _left;
    int colStart = 0;
   for(int ii = 0; ii < _numCols; ++ii ){
        _cols[ii] = MatrixCol<T>(colStart, colLength, index);
        
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
const BandedMatrixIterator<T> BandedMatrix<T>::begin(){
    return _begin;
}

template <class T>
const BandedMatrixIterator<T> BandedMatrix<T>::end(){
    return _end;
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
int BandedMatrix<T>::size(){
    return _size;
}

BandedAligner::BandedAligner(string ref){
    _ref = ref;
    _initialized = false;
}

void BandedAligner::initialize(string test){
    _test = test;
    _initialized = true;
    _matrix = BandedMatrix<int>((int)_ref.size(), (int)_test.size());

}

int BandedAligner::align(){
    
    BandedMatrixIterator<int> itr = _matrix.begin();
    int ii = 0;
    int min = -16777216;
    int x;
    int y;
    for(; ii < _matrix.size(); ++ii){

        pair<int, int> coords = itr.coordinates();
        x = coords.first;
        y = coords.second;

        int diag = min;
        int up = min;
        int left = min;

        if( x > 0 && y > 0 ){
            diag = _matrix[x - 1][y - 1] +
                (_ref[x] == _test[y] ? 1 : -1);
        }
        
        if( _matrix.inbounds(x - 1, y) ){
           left = _matrix[x-1][y] - 1; 
        }

        if( _matrix.inbounds(x, y - 1) ){
            up = _matrix[x][y-1] - 1;
        }

        _matrix[x][y] = max(diag, max(left, up));
        ++itr;
    }
    return _matrix[x][y];
}

int main(int argc, void **argv){

    BandedAligner alg("AAAAAAAAAAAA");
       alg.initialize("AAAATTTTAAAA");
    printf("%d\n", alg.align());
    return 0;
}

