#include <stdio.h>

#include "bandedalignment.h"

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
void BandedMatrixIterator<T>::test(){
    int ii = 1;
    while( hasNext() ){
        ++(*this);
        _index++;
        ii++;
    }

    printf("Expected %d cells, got %d.\n", _size, ii);

    while( hasPrev() ){
        _index--;
        --(*this);
    }

}

template <class T>
BandedMatrixIterator<T>& BandedMatrixIterator<T>::operator++() {
    _y++;
    _index++;
    if( _y - _currCol._offset >= _currCol._length ){
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

    assert(_x >= 0 && _y >= 0);
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
    _dat     = NULL;
}

template <class T>
BandedMatrix<T>::BandedMatrix(int x, int y){
    _numCols = x;
    _numRows = y;
    _left  = (_numCols < _numRows ? _numRows - _numCols : 1);
    _right = (_numRows < _numCols ? _numCols - _numRows : 1);

    _cols = vector<MatrixCol<T> >(_numCols);
    _startPos = vector<int>(_numCols);
    _dat = NULL;
    resize(false, false); // makeCols() is called from resize
    }

template <class T>
BandedMatrix<T>::~BandedMatrix(){
    if( _dat ){ free(_dat); }
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
    printf("Building the columns... ");
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
    printf("done!\n");
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
        printf("Grew the left...\n");
        _left = min(_left * 2, _numRows - 1);
    }

    if( gRight ){
        printf("Grew the right...\n");
        _right = min(_right * 2, _numCols - 1);
    }

    /*
     * Lop off the symmetric component of the matrix and
     * compute the size of the area to the left of the diagonal
     * and to the right of it.
     */
    int upperTriangle = _numCols - 1 - _right;
    int lowerTriangle = _numRows - 1 - _left;
    printf("%d -> ", _size);
    _size = 
        (_numCols * _numRows) -
        (upperTriangle * (upperTriangle + 1)) / 2 -
        (lowerTriangle * (lowerTriangle + 1)) / 2;
    assert(_size > 0);  
    printf("%d\n", _size);
    if( _dat ){ free(_dat); }
    _dat = (T*) malloc(sizeof(T) * _size);

    makeCols();
    _begin = BandedMatrixIterator<T>(_cols, _size);
    _end   = BandedMatrixIterator<T>(_cols, _size);

    //_begin.test();
    return _size;
}

template <class T>
int BandedMatrix<T>::size(){
    return _size;
}

BandedAligner::BandedAligner(string ref){
    _ref = ref;
    _initialized = false;
    _aligned = false;
}

void BandedAligner::initialize(string test){
    _test = test;
    _initialized = true;
    _aligned = false;
    _matrix = BandedMatrix<int>((int)_ref.size(), (int)_test.size());
    setBounds(0,0,0);
}

int BandedAligner::step(){

    BandedMatrixIterator<int> itr = _matrix.begin();
    int ii = 0;
    int min = -8388608;
    int x = 0;
    int y = 0;
    int prevX = -1;
    int upperScore = min;
    int lowerScore = min;
    int columnBest = min;
    int score = min;
    int maxX = _matrix.cols();
    int maxY = _matrix.rows();

    for(; ii < _matrix.size(); ++ii){

        assert(x >= 0 && x < maxX);
        assert(y >= 0 && y < maxY);
        if( x > prevX ){
            lowerScore = score;
        }

        pair<int, int> coords = itr.coordinates();
        x = coords.first;
        y = coords.second;

        int diag = min;
        int up = min;
        int left = min;

        if( x > 0 && y > 0 ){
            diag = _matrix[x - 1][y - 1] +
                (_ref[x] == _test[y] ? 1 : -1);
        } else {
            diag = (_ref[x] == _test[y] ? 1 : - 1) - x - y;
        }

        if( _matrix.inbounds(x - 1, y) ){
            left = _matrix[x-1][y] - 1; 
        }

        if( _matrix.inbounds(x, y - 1) ){
            up = _matrix[x][y-1] - 1;
        }

        score = max(diag, max(left, up));
        _matrix[x][y] = score;

        columnBest = max(columnBest, score);

        // we've rolled over to a new column
        if( x > prevX ){
            if(x != 0 &&
               x != maxX &&
               y != 0 &&
               y != maxY ){
                bool growLeft, growRight;

                if( upperScore >= columnBest ){
                    printf("Upper score %d beats %d\n", upperScore, columnBest);
                    growRight = true;
                }

                if( lowerScore >= columnBest ){
                    printf("Lower score %d beats %d\n", lowerScore, columnBest);
                    growLeft = true;
                }

                if( growRight || growLeft ){
                    printf("Resizing...\n");
                    _matrix.resize(growLeft, growRight);

                    return min;
                }
            }
            //reset our variables
            upperScore = score;
            prevX = x;
            columnBest = min;
        }

        printf("[%d][%d]: %d (d: %d u: %d l: %d)\n", x, y, _matrix[x][y], diag, up, left);
        ++itr;
    }
    _aligned = true;
    setBounds(x,y, score);
    return _matrix[x][y];
}

int BandedAligner::align(){
    int score = 0;
    while( !_aligned ){
        score = step();
    }
    return score;
}

void BandedAligner::setBounds(int x, int y, int score){

    int maxX = _matrix.cols();
    int maxY = _matrix.rows();
    int diagonal = min(maxX - x, maxY - y);
    int horizontal = max(maxX - x, maxY - y);
    _upperBound = score + diagonal - (horizontal - diagonal);
    _lowerBound = score - (maxX - x) - (maxY - y);
}

int main(int argc, void **argv){

    BandedAligner alg("AAAAAAAAAAAA");
       alg.initialize("GGGGGGGGGAAA");
    printf("%d\n", alg.step());
    printf("%d\n", alg.step());
    printf("%d\n", alg.step());
    printf("%d\n", alg.step());
    printf("%d\n", alg.step());
    printf("%d\n", alg.step());
    printf("%d\n", alg.step());
    return 0;
}

