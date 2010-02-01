// k-band alignment
//   Sequence alignment using a diagonal band.
//   Where n,m are inputs and e is the edit distance:
//     O(e^2) + O(n) + O(m)

#include "kalign.h"

using namespace std;

inline short int Alignment::getScore(
        unsigned short int x,
        unsigned short int y,
        unsigned short int index){

    short int diag, top, prv;
    short int di, ti, li;
    unsigned short int overflow  = 0;
    unsigned short int underflow = 0;

    // Get the overflow for the previous row
    if( right + y - max_cols > 0 ){
        overflow = right + y - max_cols;
    }

    // How much did the previous row underflow by?
    if( y - left < 1 ){
        underflow = left - y;
    }

    // If the cell to our left exists, get the value
    if( y - x < left && x > 0 ) {
        li = index - 1;
        prv = dp_matrix[li] + gap;
    } else {
        li = -1;
        prv = MIN_INT;
    }

    // Do the same with the cell above us
    if( x - y < right && y != 0 ){
        ti = index - ( left - underflow) - ( right - overflow);
        top = dp_matrix[ti] + GAP;
    } else {
        ti = -1;
        top = MIN_INT;
    }

    // and finally the diagonal
    if( y > 0 && x > 0 ){
        di = index - ( left - underflow) - ( right - overflow) - 1;
        diag = dp_matrix[di];
        diag += ( row_seq[y] == col_seq[x] ? match : mismatch );
    } else {
        diag = gap * x + gap * y + ( row_seq[y] == col_seq[x] ? match : mismatch );
        di = -1;
    }
    // take the max of the three neighbors
    dp_matrix[index] = (diag > top ? (diag > prv ? diag : prv) : (top > prv ? top : prv));

    //printf("x:%3d y:%3d [%3d] %3d %3d (%3d %3d %3d [%c, %c] %s)\n", x, y, index, diag, dp_matrix[index], di, ti, li, row_seq[y], col_seq[x], col_seq);
    return dp_matrix[index];
}

short int Alignment::globalScore(){
    if( aligned ){
        return score;
    } else {
        return lBound;
    }   
}

Alignment::Alignment( string seqA, string seqB ){
    dp_matrix = NULL; 
    max_rows = seqA.size();
    max_cols = seqB.size();

    if( max_rows < max_cols ){
        Alignment( seqB, seqA );
        return;
    }

    row_seq = (char*) malloc(sizeof(char) * max_rows);
    col_seq = (char*) malloc(sizeof(char) * max_cols);

    strcpy(row_seq, seqA.c_str());
    strcpy(col_seq, seqB.c_str());

    right = 1;
    left  = (max_cols < max_rows ? max_rows - max_cols : 1); // ensure that we reach the bottom corner

    uBound = setUpperBound(0,0,0);
    lBound = setLowerBound(0,0,0);
    score = lBound;

    aligned = false;
    failure = false;
    resize(1,1);
}

Alignment::~Alignment(){
    if( dp_matrix ){ free(dp_matrix); dp_matrix = NULL; }
}

string Alignment::seqA(){
    return string(row_seq);
}

string Alignment::seqB(){
    return string(col_seq);
}

bool Alignment::canStep(){

    return (!aligned && !failure);

}

bool Alignment::isAligned(){
    return aligned;
}

short int Alignment::upperBound(){ return uBound; }
short int Alignment::lowerBound(){ return lBound; }

short int Alignment::setUpperBound(short int x, short int y, int score){
    int d = min(max_cols - x, max_rows - y); // get the maximum score along the diagonal
    int h = max(max_cols - x - d, max_rows - y - d);

    assert(h >= 0);
    uBound = score + d * match - h * gap;
    return uBound;
}

short int Alignment::setLowerBound(short int x, short int y, int score){
    int perimeter     = (max_cols - x) + (max_rows - y);
    int mismatchPath  = min(max_cols - x, max_rows - y);
    int mismatchPerimeter = max(max_cols - x - mismatchPath, max_rows - y - mismatchPath);

    assert(perimeter >= 0);
    assert(mismatchPath >= 0);
    assert(mismatchPerimeter >= 0);

    int perimeterScore = perimeter * gap;
    int mismatchScore = mismatchPath * mismatch + mismatchPerimeter * gap;

    lBound = score + max(perimeterScore, mismatchScore);
    return lBound;
}

short int Alignment::resize(short int l, short int r){
    left = l;
    right = r;
    short int d = max_rows - max_cols;    
    size = max_rows 
        * (left + right + 1) 
        - (left * (left + 1) / 2) 
        - (right * (right + 1) / 2)
        - (right * d)
        - (d * (d+1) / 2); // account for left hand overflow and the diagonal

    //printf("Left: %d\tRight: %d\td: %d\tSize: %d\n", left, right, d, size);
    if(dp_matrix){ free(dp_matrix); }
    dp_matrix = (short int*) malloc(sizeof(short int) * size);
    return size;
}
short int Alignment::step(){

    printf("[%d, %d]\n", lowerBound(), upperBound());
    assert( canStep() );
    //printf("Size at Step: %d\n", size);
    /*
       printf("Rows: %d\tCols: %d\tLeft: %d\tRight: %d\n", max_rows, max_cols, left, right);
       printf("%d - %d - %d - %d - %d = %d\n", (max_rows * (left + right + 1)), (left * (left + 1) / 2), (right * (right + 1) / 2), (right * d), (d * (d+1) / 2), size); 
     */
    bool boundary = false;
    int ii = 0;

    for(int y = 0; y < max_rows; ++y ){
        int start = (y >= left ? y - left : 0 );
        int stop  = (max_cols < y + right + 1 ? max_cols : y + right + 1);
        int lmax = MIN_INT;
        int rmax = MIN_INT;
        int line_max = MIN_INT;

        int max_x, max_y;
        for( int x = start ; x < stop; ++x){ 
            // add a matrix entry
            int t = getScore(x, y, ii);

            // we need to know if we've run into a boundary
            if( x == start && x > 0){ lmax = t; }
            if( x == stop - 1 && x < max_cols - 1){ rmax = t; }
            if( t > line_max ){
                line_max = t; 
                max_x = x; 
                max_y = y; 
            }

            ++ii; // advance the index
        }

        // if we've hit a boundary, short circuit our way out of the loop
        // and increase the boundary sizes
        if( max_cols > y + right + 1 ){

            if( line_max == lmax && left != max_cols / 2){ 
                boundary = true;
                left = min(left * 2, max_cols / 2);
            }

            if( line_max == rmax && right != max_cols / 2){
                boundary = true;
                right = min(right * 2, max_cols / 2);
            }
            
            if( boundary ){ 
                // printf("Hit the boundary...(%d < %d - %d)\n", y, max_rows, right);
                resize(left, right);
                setUpperBound(max_x, max_y, line_max);
                setLowerBound(max_x, max_y, line_max);
                return MIN_INT;
            }
        }
    }


    int s = dp_matrix[size - 1];
    // printf("%d\t%d\n", ii, size);
    //printf("Rows: %d\tCols: %d\tLeft: %d\tRight: %d\n", max_rows, max_cols, left, right);
    //printf("%d - %d - %d - %d - %d = %d\n", (max_rows * (left + right + 1)), (left * (left + 1) / 2), (right * (right + 1) / 2), (right * d), (d * (d+1) / 2), size); 
    //printf("Size: %d\tScore: %d\tLeft: %d\tRight: %d\n", size, s, left, right);

    assert(boundary || ii == size);

    if( !boundary ){ 
        aligned = true; 
        score = s;
        return s; 
    }

    return MIN_INT;
}

short int Alignment::align(){

    int score;
    while(canStep()){
        score = step();
    }

    return score;
}

Alignment* Alignment::round_robin( list<string> references, string target ){
    queue<Alignment*> robin;
    list<string>::iterator refItr;
    Alignment *bestAlign;
    Alignment *algn;
    int bestScore = MIN_INT;
    

    for( refItr = references.begin(); refItr != references.end(); refItr++ ){
        algn = new Alignment((*refItr), target);
        robin.push(algn);
    }

    int robinSize = robin.size();
    while( !robin.empty() ){
        // ensure that I'm not growing the robin accidentally.
        assert( robin.size() <= robinSize );
        robinSize = robin.size();

        

        algn = robin.front();
        robin.pop();

        // If this alignment might be better than the best lower bound, continue.
        if( algn->upperBound() >= bestScore ){

            assert(algn->canStep());
            algn->step();
            
            printf("[%d, %d]\n", algn->lowerBound(), algn->upperBound());

            if( algn->lowerBound() > bestScore ){
                bestScore = algn->lowerBound();
                bestAlign = algn;
            }

            if( algn->upperBound() >= bestScore && algn->canStep() ){
                robin.push(algn);
            }
        } else {
            delete algn;
        }                 
    }

    if( !bestAlign->isAligned() ){
        bestAlign->align();
    }

    return bestAlign;
}


int main( int argc, void** argv ){

    list<string> ref;
    ref.push_front(string("abcdefghijkl"));
    ref.push_front(string("abceefghijkl"));
    ref.push_front(string("abceefghijxl"));
    ref.push_front(string("abcdefghxxxijkl"));
    ref.push_front(string("xxxxxxxxxxxx"));

    list<string>::iterator itr;

    for( itr = ref.begin() ; itr != ref.end() ; itr++ ){
        Alignment n = Alignment(*itr, string("abcdefghijkl"));
        n.align();
        
        printf("%s: %d\n", (*itr).c_str(), n.globalScore());
    }
/*
    Alignment* a = Alignment::round_robin(ref, string("abcdefghijkl"));
    printf("[%d,%d]: %d %s\n", a->lowerBound(), a->upperBound(), a->globalScore(), a->seqA().c_str());
    */
    return 0;

}
