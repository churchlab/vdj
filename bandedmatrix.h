#include <assert.h>
#include <signal.h>

#include <utility>
#include <vector>

template <class T>
class MatrixCol;

template <class T>
class BandedMatrixIterator
{
    public:
        BandedMatrixIterator();
        BandedMatrixIterator(const BandedMatrixIterator<T> &itr);
        BandedMatrixIterator(std::vector<MatrixCol<T> > &itr, int size);
        ~BandedMatrixIterator();
        //FIXME these operators are broken. fix em.
        bool operator== (BandedMatrixIterator<T> &itr); 
        bool operator!= (BandedMatrixIterator<T> &itr); 
        
        bool hasNext();
        bool hasPrev();

        //FIXME these operators only work as prefix ops.
        BandedMatrixIterator<T>& operator++();
        BandedMatrixIterator<T>& operator--();
        T& operator*() const;
        std::pair<int, int> coordinates();
    private:
        T* _dat;
        std::vector<MatrixCol <T> > _cols;
        MatrixCol<T> _currCol;
        int _size;
        int _x;
        int _y;
        int _index;
};

template <class S>
class MatrixCol
{
    public:
        MatrixCol();
        MatrixCol(int, int, S*);
        S& operator[](const int);
    protected:
        friend class BandedMatrixIterator<S>;
        S* _dat;
        int _offset;
        int _length;
};

template <class T>
class BandedMatrix
{

    public:
        BandedMatrix(int, int);
        MatrixCol<T> operator[](const int);
        const BandedMatrixIterator<T> begin();
        const BandedMatrixIterator<T> end();
        bool inbounds(int, int);
        int resize(bool, bool);
    private:
        T* _dat; //TODO use shared swap space later;
        T  _val;
        std::vector<MatrixCol<T> > _cols;
        std::vector<int> _startPos; 
        void makeCols();
        T* getPtr(int, int);

        int _numRows;
        int _numCols;

        int _size;
        int _left;
        int _right;

        BandedMatrixIterator<T> _begin;
        BandedMatrixIterator<T> _end;
};

