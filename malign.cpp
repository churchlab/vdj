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
        entries.insert(pair<string,MAlignerEntry*>(string(seqs[ii]),en));
    }

}

priority_queue<MAlignerEntry*> MAligner::align(char* input){
    queue<MAlignerEntry*> robin; 
    map<string, MAlignerEntry*>::iterator entry_itr;

    int len = strlen(input);

    for( entry_itr = entries.begin(); entry_itr != entries.end(); entry_itr++ ){
        MAlignerEntry *e = (*entry_itr).second;
        e->initialize(input, len);
        robin.push(e);
    }

    return this->roundRobin(robin);

}

priority_queue<MAlignerEntry*> MAligner::alignWith(char* input, int nnames, char** names){

    queue<MAlignerEntry*> robin;
    int len = strlen(input);

    map<string,MAlignerEntry*>::iterator search_itr;
    for( int ii = 0; ii < nnames ; ++ii ){
        search_itr = entries.find(string(names[ii]));
        if( search_itr != entries.end() ){
            (*search_itr).second->initialize(input, len);
            robin.push( (*search_itr).second );
        }
    }

    return this->roundRobin(robin);

}

priority_queue<MAlignerEntry*> MAligner::roundRobin(queue<MAlignerEntry*> robin){
    
    priority_queue<MAlignerEntry*> results;
    int bestLowerBound = MINVAL;
    
    while( !robin.empty() ){
        MAlignerEntry *mae = robin.front();
        robin.pop();

        if( mae->upperBound() > bestLowerBound ){
            mae->step();

            if( mae->lowerBound() > bestLowerBound ){
                bestLowerBound = mae->lowerBound();
            }

            if( mae->upperBound() >= bestLowerBound
                && !mae->isAligned() ){
                robin.push(mae);
            }
        }

        if( mae->isAligned() ){
            results.push(mae); 
        }
    }

    return results;
}

bool  MAligner::initialize(char* input){
    map<string, MAlignerEntry*>::iterator entry_itr;

    int len = strlen(input);

    for( entry_itr = entries.begin() ; entry_itr != entries.end() ; entry_itr++ ){
        (*entry_itr).second->initialize(input, len);
    }



    return true;
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
    this->setLowerBound(0,0,0);
    this->setUpperBound(0,0,0);
}

bool MAlignerEntry::isAligned() {
    return aligned;
}

const bool MAlignerEntry::operator< (MAlignerEntry& entry){
    return this->getScore() > entry.getScore(); 
}
const bool MAlignerEntry::operator> (MAlignerEntry& entry){
    return this->getScore() < entry.getScore();
}

int MAlignerEntry::getScore(){
    return score; 
}

int MAlignerEntry::upperBound(){
    return this->uBound;
}

int MAlignerEntry::lowerBound(){
    return this->lBound;
}

int MAlignerEntry::setUpperBound(int x, int y, int score){
    int d = min(this->max_cols - x, this->max_rows - y); // get the maximum score along the diagonal
    int h = max(this->max_cols - x - d, this->max_rows - y - d);

    assert(h >= 0);
    this->uBound = score + d * this->match - h * this->gap;
    return uBound;
}

int MAlignerEntry::setLowerBound(int x, int y, int score){
    int perimeter     = (this->max_cols - x) + (this->max_rows - y);
    int mismatchPath  = min(this->max_cols - x, this->max_rows - y);
    int mismatchPerimeter = 
          max(this->max_cols - x - mismatchPath, this->max_rows - y - mismatchPath);

    assert(perimeter >= 0);
    assert(mismatchPath >= 0);
    assert(mismatchPerimeter >= 0);

    int perimeterScore = perimeter * this->gap;
    int mismatchScore = mismatchPath * this->mismatch + mismatchPerimeter * this->gap;

    this->lBound = score + max(perimeterScore, mismatchScore);
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
                this->grow(grow_left, grow_right);
                this->setUpperBound(max_x, max_y, line_max);
                this->setLowerBound(max_x, max_y, line_max);
            
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
