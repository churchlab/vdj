// Multi-aligner

#include "malign.h"

using namespace std;
/*
 * MAligner
 */

MAligner::MAligner(int nseq, char** seqs, char** names){

   global_matrix = (dp_matrix*) malloc(sizeof(dp_matrix));
   global_matrix->matrix = NULL;
   global_matrix->size   = 0;

    for( int ii = 0; ii < nseq; ++ii ){
        MAlignerEntry *en = new MAlignerEntry(names[ii], seqs[ii], global_matrix);
        entries.push_back(en);
    }

}

const char* MAligner::align(char* input){
    list<MAlingerEntry*>::iterator entry_itr;

    for( entry_itr = entries.begin() ; entry_itr != entries.end() ; entry_itr++ ){
        (*entry_itr)->initialize(input);
    }


}

/*
 * MAlignerEntry
 */
MAlignerEntry::MAlignerEntry(char *nm, char *seq, dp_matrix *dp){
    this->name = nm;
    this->refSequence = seq;
    this->dpm = dp;
    this->max_rows = strlen(seq);
    this->max_cols = -1; // initialize this when we load it with a sequence

    this->gap      = -1;
    this->match    =  1;
    this->mismatch = -1;

    this->initialized = false;
}

int MAlignerEntry::upperBound(){
    return uBound;
}

int MAlignerEntry::lowerBound(){
    return lBound;
}

int MAlignerEntry::grow(bool growLeft, bool growRight){
                
    if( growLeft ){
        this->left = min(this->left * 2, this->max_cols / 2);
    }

    if( growRight ){
        this->right = min(this->right * 2, this->max_cols / 2);
    }
            
    int l = this->left;
    int r = this->right;
    int d = this->max_rows - this->max_cols;    
    this->size = 
        this->max_rows 
        * (l + r + 1) 
        - (l * (l + 1) / 2) 
        - (r * (r + 1) / 2)
        - (r * d)
        - (d * (d+1) / 2); // account for left hand overflow and the diagonal
    
    // Ask for a bigger scratch space if we need it
    if(this->dpm->size < this->size ){ 
        free(this->dpm->matrix);
        this->dpm->matrix = (int*) malloc(sizeof(int) * this->size);
        this->dpm->size = this->size;
    }
    
    return this->size;
}

bool MAlignerEntry::initialize(char* testSeq, int len=-1){
    this->testSequence = testSeq;
    this->max_cols = (len > 0 ? len : strlen(testSeq));
    this->initialized = true;
    this->left = 1;
    this->right = 1;
    grow(false, false);
    return true;
}

bool MAlignerEntry::step(){

    int max_rows = this->max_rows;
    int max_cols = this->max_cols;

    int *prevRow = NULL;
    int *thisRow = this->dpm->matrix;

    bool boundary = false;

    const int left  = this->left;
    const int right = this->right;

    for(int y = 0; y < max_rows; ++y ){
        int start = (y >= left ? y - left : 0 );
        int stop  = (max_cols < y + right + 1 ? max_cols : y + right + 1);
        int lmax = MINVAL;
        int rmax = MINVAL;
        int line_max = MINVAL;

        int max_x, max_y;
        int rowOffset = 0;
        for( int x = start ; x < stop; ++x){ 
            // add a matrix entry
            int t = this->scoreDP(x, y, rowOffset, prevRow, thisRow);

            // we need to know if we've run into a boundary
            if( x == start && x > 0){ lmax = t; }
            if( x == stop - 1 && x < max_cols - 1){ rmax = t; }
            if( t > line_max ){
                line_max = t; 
                max_x = x; 
                max_y = y; 
            }

            ++rowOffset;
        }

        prevRow = thisRow;
        thisRow = thisRow + rowOffset; //TODO off by one?

        // if we've hit a boundary, short circuit our way out of the loop
        // and increase the boundary sizes
        if( max_cols > y + right + 1 ){
            bool grow_left = false;
            bool grow_right = false;

            if( line_max == lmax
                && this->left != this->max_cols / 2 ){
                boundary = true;
                grow_left = true;
            }

            if( line_max == rmax
                && this->right != this->max_cols / 2){
                boundary = true;
                grow_right = true;
            }
            
            if( boundary ){ 
                grow(grow_left, grow_right);
             //   setUpperBound(max_x, max_y, line_max);
             //   setLowerBound(max_x, max_y, line_max);
            
                return MINVAL;
            }
        }
    }

    int s = this->dpm->matrix[this->size - 1];

    if( !boundary ){ 
        this->aligned = true; 
        this->score = s;
        return s; 
    }

    return MINVAL;
}

int MAlignerEntry::scoreDP(
        int x,
        int y,
        int rowOffset,
        int *prevRow,
        int *thisRow) {

    int leftVal;
    int upVal;
    int upleftVal;

    char* row_seq = this->refSequence;
    char* col_seq = this->testSequence;

    if(rowOffset > 0 ){
        leftVal = thisRow[rowOffset - 1] + this->gap;
    } else {
        leftVal = MINVAL;
    }

    if(y > 0 && x > 0){
        upleftVal = prevRow[rowOffset] + (row_seq[y] == col_seq[x] ? this->match : this->mismatch );
    }

    if((y > 0)
       && ( prevRow + rowOffset + 1 < thisRow )){
        upVal = prevRow[rowOffset + 1] + this->gap;
    } else {
        upVal = MINVAL;
    }

    thisRow[rowOffset] = max(upleftVal, max(leftVal, upVal));
    return thisRow[rowOffset];
}
