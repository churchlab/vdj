#include "bandedalignment.h"

using namespace std;

BandedAlignment::BandedAlignment(string ref){
    _ref = ref;
    _initialized = false;
}

void BandedAlignment::initialize(string test){
    _test = test;
    _initialized = true;
    _matrix = BandedMatrix<int>((int)_ref.size(), (int)_test.size());

}

int BandedAlignment::align(){
    
    BandedMatrixIterator<int> itr = _matrix.begin();
    int ii = 0;
    int min = -16777216;
    for(; ii < _matrix.size(); ++ii){

        pair<int, int> coords = itr.coordinates();
        int x = coords.first;
        int y = coords.second;

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

}
