// k-band alignment
//   Sequence alignment using a diagonal band.
//   Where n,m are inputs and e is the edit distance:
//     O(e^2) + O(n) + O(m)

#include "kalign.h"

using namespace std;

inline int Alignment::getScore(
        int x,
        int y,
        int index){
    
    int gap = -1;
    int match = 1;
    int mismatch = -1;


    int diag, top, prv;
    int di, ti, li;
    int overflow  = 0;
    int underflow = 0;

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
        top = dp_matrix[ti] + gap;
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


    // Sanity check our indices
    assert( li >= -1 && li <= size );
    assert( di >= -1 && li <= size );
    assert( ti >= -1 && li <= size );
    assert( index >= 0 && index <= size );

    //Dump everything that we could plausibly want to know about the program
    // state at this point:
    // sprintf(stderr, "x:%3d y:%3d [%3d] %3d %3d (%3d %3d %3d [%c, %c] %s)\n", x, y, index, diag, dp_matrix[index], di, ti, li, row_seq[y], col_seq[x], col_seq);
    
    return dp_matrix[index];
}

int Alignment::globalScore(){
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

    
    row_seq = (char*) malloc(sizeof(char) * (max_rows + 1));
    col_seq = (char*) malloc(sizeof(char) * (max_cols + 1));

    strcpy(row_seq, seqA.c_str());
    strcpy(col_seq, seqB.c_str());

    right = 1;
    left  = (max_cols < max_rows ? max_rows - max_cols : 1); // ensure that we reach the bottom corner

    uBound = setUpperBound(0,0,0);
    lBound = setLowerBound(0,0,0);
    score = lBound;

    aligned = false;
    resize(1,1);
}

Alignment::~Alignment(){
    if( dp_matrix ){ free(dp_matrix); dp_matrix = NULL; }
    free(row_seq);
    free(col_seq);
}

string Alignment::seqA(){
    return string(row_seq);
}

string Alignment::seqB(){
    return string(col_seq);
}

bool Alignment::canStep(){
    return (!aligned);
}

bool Alignment::isAligned(){
    return aligned;
}

int Alignment::upperBound(){ return uBound; }
int Alignment::lowerBound(){ return lBound; }

int Alignment::setUpperBound(int x, int y, int score){
    int match = 1;
    int gap = -1;
    int d = min(max_cols - x, max_rows - y); // get the maximum score along the diagonal
    int h = max(max_cols - x - d, max_rows - y - d);

    assert(h >= 0);
    uBound = score + d * match - h * gap;
    return uBound;
}

int Alignment::setLowerBound(int x, int y, int score){
    int match = 1;
    int gap = -1;
    int mismatch = -1;
    
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

int Alignment::resize(int l, int r){
    left = l;
    right = r;
    int d = max_rows - max_cols;    
    size = max_rows 
        * (left + right + 1) 
        - (left * (left + 1) / 2) 
        - (right * (right + 1) / 2)
        - (right * d)
        - (d * (d+1) / 2); // account for left hand overflow and the diagonal
    
    if(dp_matrix){ free(dp_matrix); }
    // printf("Size: %d %d\n", size, sizeof(int)* size);
    dp_matrix = (int*) malloc(sizeof(int) * size);
    return size;
}
int Alignment::step(){

    assert( canStep() );
    
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

        // if we've hit a boundary, circuit our way out of the loop
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
                resize(left, right);
                setUpperBound(max_x, max_y, line_max);
                setLowerBound(max_x, max_y, line_max);
                return MIN_INT;
            }
        }
    }


    int s = dp_matrix[size - 1];

    // Ensure that my matrix diagonal size calculations were correct.
    assert(boundary || ii == size);

    if( !boundary ){ 
        aligned = true; 
        score = s;
        return s; 
    }

    return MIN_INT;
}

int Alignment::align(){

    int score;
    while(canStep()){
        score = step();
    }

    return score;
}

Alignment* Alignment::round_robin( list<string> references, string target ){
    queue<Alignment*> robin;
    list<Alignment*>  alignedList;
    list<string>::iterator refItr;
    Alignment *algn;
    int bestLowerBound = MIN_INT;
    

    for( refItr = references.begin(); refItr != references.end(); refItr++ ){
        algn = new Alignment((*refItr), target);
        robin.push(algn);
    }

    unsigned int robinSize = robin.size();
    while( !robin.empty() ){
        // ensure that I'm not growing the robin accidentally.
        assert( robin.size() <= robinSize );
        robinSize = robin.size();

        // initialize...
        algn = robin.front();
        robin.pop();

        // If this alignment might be better than the best lower bound, continue.
        if( algn->upperBound() >= bestLowerBound ){

            // ensure that we don't have any dead alignments
            assert(algn->canStep());
            algn->step();
            
            // Should we raise our lower bound?
            if( algn->lowerBound() > bestLowerBound ){
                bestLowerBound = algn->lowerBound();
            }

            // If we can beat the lower bound, chuck it in the queue
            if( algn->upperBound() >= bestLowerBound ){
                if( algn->isAligned() ){
                    alignedList.push_back(algn);
                } else {
                    robin.push(algn); // Insert incomplete alignments into the robin
                }
            } else {
                delete algn;
            }
        } else {
            delete algn;
        }
    }

    // pick out the best alignment
    // NOTE it's conceivable that we might want to examine suboptimal alignments
    list<Alignment*>::iterator algItr;
    Alignment* bestAlignment = NULL;
    for( algItr = alignedList.begin(); algItr != alignedList.end(); algItr++ ){
        if( bestAlignment ){
            if( bestAlignment->globalScore() < (*algItr)->globalScore() ){
                delete bestAlignment;
                bestAlignment = (*algItr);
            } else {
                delete (*algItr);
            }
        } else {
            bestAlignment = (*algItr);
        }
    }
    return bestAlignment;
}


int main( int argc, char** argv ){
/*
    list<string> ref;
    ref.push_front(string("abcdefghijkl"));
    ref.push_front(string("abceefghijkl"));
    ref.push_front(string("abceefghijxl"));
    ref.push_front(string("abcdefghxxxijkl"));
    ref.push_front(string("xxxxxxxxxxxx"));

    list<string>::iterator itr;
    
    printf("Performing all alignments: \n");
    for( itr = ref.begin() ; itr != ref.end() ; itr++ ){
        Alignment n = Alignment(*itr, string("abcdefghijkl"));
        n.align();
        
        printf("%s: %d\n", (*itr).c_str(), n.globalScore());
    }

    printf("-- -- -- -- -- -- -- --\n");
   
    */

    for( int jj = 0; jj < 10000; ++jj ){
    string x;
    string y;
    for( int ii = 0; ii < 500; ++ii ){
        x += "x";
        y += "y";
    }
    Alignment n = Alignment(x, y);
    n.align();

    }   
    /* Alignment* a = Alignment::round_robin(ref, string("abcdefghijkl"));
    printf("[%d,%d]: %d %s\n", a->lowerBound(), a->upperBound(), a->globalScore(), a->seqA().c_str());
    */

    return 0;

}
