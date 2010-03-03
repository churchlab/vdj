#include <assert.h>
#include <signal.h>

#include <utility>
#include <vector>

template <class S>
class MatrixRow
{
    public:
        MatrixRow(int, int, S*);
        S& operator[](const int);
    private:
        int _offset;
        int _length;
        S*  _dat;
};

template <class T>
class BandedMatrix
{

    public:
        BandedMatrix(int, int);
        MatrixRow<T> operator[](const int);
        T get(int, int);
        T set(int, int, T);
        int resize(bool, bool);
    private:
        T* _dat; //used shared swap space later;
        T  _val;
        std::vector<T*> _cols;
        std::vector<int> _startPos; 
        void makeCols();
        T* getPtr(int, int);

        int _numRows;
        int _numCols;

        int _size;
        int _left;
        int _right;
};

