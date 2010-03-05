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
    _maxX = 0;
    _maxY = 0;
}

template <class T>
BandedMatrixIterator<T>::BandedMatrixIterator(const BandedMatrixIterator<T> &itr){
    _cols = itr._cols;
    _currCol = itr._currCol;
    _size = itr._size;
    _index = itr._index;
    _x = itr._x;
    _y = itr._y;
    _maxX = itr._maxX;
    _maxY = itr._maxY;
    _dat = itr._dat;
}

template <class T>
BandedMatrixIterator<T>::BandedMatrixIterator(vector<MatrixCol<T> > &cols, int size, int maxRows, int maxCols){
    _currCol = cols[0];
    _dat = _currCol._dat;
    _cols = cols;
    _x = 0;
    _y = 0;
    _size = size;
    _maxX = maxCols - 1;
    _maxY = maxRows - 1;
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
    return _dat[_index];
}

template <class T>
pair<int, int> BandedMatrixIterator<T>::coordinates(){
    return pair<int, int>(_x,_y);
}

template <class T>
bool BandedMatrixIterator<T>::boundary(){
  
    bool upBound = false;
    bool leftBound = false;

    if( _x > 0 && _y > 0 ){
        upBound = ( 0 == _y - _currCol._offset );
    
        MatrixCol<T> prevCol = _cols[_x - 1];
        int poffset = prevCol._offset;
        int plen    = prevCol._length;
        leftBound = (poffset + plen - 1 < _y);
    }

    return (upBound || leftBound);
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
    int colLength = 1 + _left;
    T* index = _dat; 
    int allowableRight = 0;
    int allowableLeft  = _left;
    int colStart = 0;
    for(int ii = 0; ii < _numCols; ++ii){
        
        int d = 1; //(ii < _numRows - 1 ? 1 : 0);

        allowableLeft  = _left  + min(0, _numRows - ii - 1 - _left);
        allowableRight = _right + min(0, ii - _right);
        colStart = ii - allowableRight;
/*
        assert(allowableLeft >= 0);
        assert(allowableRight >= 0);
*/
        colLength = d
            + allowableLeft 
            + allowableRight;
        index += colLength;
       
        //printf("Column: %d Allowable Right: %d Allowable Left: %d Length: %d Offset: %d\n", ii, allowableRight, allowableLeft, colLength, colStart);

        _cols[ii] = MatrixCol<T>(colStart, colLength, index);
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
    if( _dat ){ free(_dat); }
    _dat = (T*) malloc(sizeof(T) * _size);

    makeCols();
    _begin = BandedMatrixIterator<T>(_cols, _size, _numRows, _numCols);
    //TODO set this iterator to the end.
    _end   = BandedMatrixIterator<T>(_cols, _size, _numRows, _numCols);

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

BandedAligner::~BandedAligner(){
}

void BandedAligner::initialize(string test){
    _test = test;
    _initialized = true;
    _aligned = false;
    _matrix = new BandedMatrix<int>((int)_ref.size(), (int)_test.size());
    setBounds(0,0,0);
}

pair<string, string> BandedAligner::getBacktrace(){

    stack<Direction> backtrace;

    if( !_aligned ){
        align();
    }

    int x = _matrix->cols() - 1;
    int y = _matrix->rows() - 1;

    while( x > 0 && y > 0 ){
        assert( _matrix->inbounds(x, y - 1 ));
        assert( _matrix->inbounds(x - 1, y - 1 ));
        assert( _matrix->inbounds(x - 1, y ));

        int up   = (*_matrix)[x][y-1];
        int left = (*_matrix)[x-1][y];
        int diag = (*_matrix)[x-1][y-1];

        printf("BT: [%d][%d] up: %d left: %d diag: %d\t", x, y, up, left, diag);    
        if( diag >= up && diag >= left )
        {
            backtrace.push(Diagonal);
            x--;
            y--;
            printf("DIAGONAL\n");
        } else if( up >= left ){
            backtrace.push(Up);
            y--;
            printf("UP\n");
        } else {
            backtrace.push(Left);
            x--;
            printf("LEFT\n");
        }
    }
    
    //It was just too tempting to write this:
    while(x-->0){
        backtrace.push(Left);
    }

    while(y-->0){
        backtrace.push(Up);
    }
    
    string refAlign;
    string testAlign;

    x = 0;
    y = 0;
    while(!backtrace.empty()){
        Direction d = backtrace.top();
        backtrace.pop();

        switch(d){
            case Diagonal:
                refAlign  += _ref[x];
                testAlign += _test[y];
                x++;
                y++;
                break;
            case Up:
                refAlign  += "-";
                testAlign += _test[y];
                y++;
                break;
            case Left:
                refAlign  += _ref[x];
                testAlign += "-";
                x++;
                break;
        }
    }
   
    printf("R: %s\nT: %s\n", refAlign.c_str(), testAlign.c_str());
    return pair<string, string>(refAlign, testAlign);


}

int BandedAligner::step(){

    BandedMatrixIterator<int> itr = _matrix->begin();
    int ii = 0;
    const int min = -8388608;
    int x = 0;
    int y = 0;
    int prevX = -1;
    int upperScore = min;
    int lowerScore = min;
    int columnBest = 0;
    int score = min;
    int maxX = _matrix->cols();
    int maxY = _matrix->rows();
    pair<int, int> coords; 
    
    bool upperBoundary = false;
    bool lowerBoundary = (_matrix->left() + 1 == maxY ? false : true);
    
    int left = min;
    int diag = min;
    int up = min;

    for(; ii < _matrix->size(); ++ii){
        x = itr.coordinates().first;

        if( x > prevX ){
            lowerScore = score;

            printf("Column rolled over. Upper: %d Lower: %d Best: %d\n", upperScore, lowerScore, columnBest);
            bool growRight = false;
            bool growLeft  = false;

            if( upperBoundary && upperScore >= columnBest ){
                growRight = true;
            }

            if( lowerBoundary && lowerScore >= columnBest ){
                growLeft = true;
            }

            if( growLeft || growRight ){
                printf("Growing Left[%s] Right[%s]\n\n\n", (growLeft?"X":" "), (growRight?"X":" "));
                _matrix->resize(growLeft, growRight);
                return min;
            }
        }

        y = itr.coordinates().second;
        
        assert(x >= 0 && x < maxX);
        assert(y >= 0 && y < maxY);

        diag = min;
        up = min;
        left = min;

        if( x > 0 && y > 0 ){
            diag = (*_matrix)[x - 1][y - 1] +
                (_ref[x] == _test[y] ? 1 : -1);
        } else {
            diag = (_ref[x] == _test[y] ? 1 : - 1) - x - y;
        }

        if( _matrix->inbounds(x - 1, y) ){
            left = (*_matrix)[x-1][y] - 1; 
        }

        if( _matrix->inbounds(x, y - 1) ){
            up = (*_matrix)[x][y-1] - 1;
        }

        score = max(diag, max(left, up));
        (*_matrix)[x][y] = score;

        if( x > prevX ){
            upperScore = score;
            upperBoundary = itr.boundary();
            columnBest = score;
            prevX = x;
        }

        columnBest = max(columnBest, score);
        lowerBoundary = itr.boundary();
        printf("[%d][%d]: %d (d: %d u: %d l: %d) %s\n", x, y, (*_matrix)[x][y], diag, up, left, (itr.boundary() ? "BOUNDARY" : ""));
        ++itr;
    }

    bool growLeft = false;
    bool growRight = false;
    if( min == left ){
        growLeft = true;
    }
    if( min == up ){
        growRight = true;
    }

    if( growLeft || growRight ){
        _matrix->resize(growLeft, growRight);
        return min;
    }
    _aligned = true;
    setBounds(x,y, score);
    _score = score;
    return (*_matrix)[x][y];
}

int BandedAligner::align(){
    int score = 0;
    while( !_aligned ){
        score = step();
    }
    return score;
}

void BandedAligner::setBounds(int x, int y, int score){

    int maxX = _matrix->cols();
    int maxY = _matrix->rows();
    int diagonal = min(maxX - x, maxY - y);
    int horizontal = max(maxX - x, maxY - y);
    _upperBound = score + diagonal - (horizontal - diagonal);
    _lowerBound = score - (maxX - x) - (maxY - y);
}


int main(int argc, char **argv){

    BandedAligner algn("AAAAAAAAAAAAAA");
       algn.initialize("AAATTTTTTAAAAA");
       algn.align();
       algn.getBacktrace();
       return 0;

}
