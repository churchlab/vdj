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

MAligner::~MAligner(){
    if( global_matrix->matrix != NULL ){
        free(global_matrix->matrix);
    }

    free(global_matrix);
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
MAlignerEntry::MAlignerEntry(const char *nm, const char *seq, dp_matrix *dp){
    
    this->name = (char*) malloc(strlen(nm) + 1);
    this->refSequence = (char*) malloc(strlen(seq) + 1);
    this->testSequence = NULL;
    strcpy(this->name, nm);
    strcpy(this->refSequence, seq);

    this->ref_length = strlen(this->refSequence);

    this->dpm = dp;
    this->num_rows = NULL;
    this->num_cols = NULL; // initialize this when we load it with a sequence

    this->uBound = MINVAL;
    this->lBound = MINVAL;
    this->matrix_size = MINVAL;

    this->gap      = -1;
    this->match    =  1;
    this->mismatch = -1;

    this->initialized = false;
    this->setLowerBound(0,0,0);
    this->setUpperBound(0,0,0);
}

MAlignerEntry::~MAlignerEntry(){
    free(this->name);
    free(this->refSequence);
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
    
    if( !this->initialized ){ return MINVAL; }
    
    int d = min(this->num_cols - x, this->num_rows - y); // get the maximum score along the diagonal
    int h = max(this->num_cols - x - d, this->num_rows - y - d);

    assert(h >= 0);
    this->uBound = score + d * this->match - h * this->gap;
    return uBound;
}

int MAlignerEntry::setLowerBound(int x, int y, int score){

    if( !this->initialized ){ return MINVAL; }

    int perimeter     = (this->num_cols - x) + (this->num_rows - y);
    int mismatchPath  = min(this->num_cols - x, this->num_rows - y);
    int mismatchPerimeter = 
          max(this->num_cols - x - mismatchPath, this->num_rows - y - mismatchPath);
/*
    assert(perimeter >= 0);
    assert(mismatchPath >= 0);
    assert(mismatchPerimeter >= 0);
*/
    int perimeterScore = perimeter * this->gap;
    int mismatchScore = mismatchPath * this->mismatch + mismatchPerimeter * this->gap;

    this->lBound = score + max(perimeterScore, mismatchScore);
    return lBound;
}

int MAlignerEntry::grow(bool growLeft, bool growRight){
                
    
    if( growLeft ){
        this->left = min(this->left * 2, this->num_cols / 2);
    }

    if( growRight ){
        this->right = min(this->right * 2, this->num_cols / 2);
    }
            
    int l = this->left;
    int r = this->right;
    int d = this->num_rows - this->num_cols;    
    this->matrix_size = 
        this->num_rows 
        * (l + r + 1) 
        - (l * (l + 1) / 2) 
        - (r * (r + 1) / 2)
        - (r * d)
        - (d * (d+1) / 2); // account for left hand overflow and the diagonal
  
    printf("diff: %d Left: %d Right: %d Rows: %d Cols: %d Size: %d\n", d, l, r, this->num_rows, this->num_cols, this->matrix_size);
    // Ask for a bigger scratch space if we need it
    if(this->dpm->size < this->matrix_size ){ 
        free(this->dpm->matrix);
        this->dpm->matrix = (int*) malloc(sizeof(int) * this->matrix_size);
        this->dpm->size = this->matrix_size;
    }
    
    return this->matrix_size;
}

bool MAlignerEntry::initialize(const char* testSeq, int len=0){
   
    this->test_length = strlen(testSeq);
    
    //printf("Initialized entry with '%s' and length %d...\n", testSeq, len);
    this->testSequence = (char*) malloc(this->test_length + 1);
    strcpy(this->testSequence, testSeq);
    
    // test will be the row
    if( this->test_length > this->ref_length ) {
        this->row_seq  = this->testSequence;
        this->num_rows = this->test_length;
        this->col_seq  = this->refSequence;
        this->num_cols = this->ref_length;
    // reference will be the row
    } else {
        this->row_seq  = this->refSequence;
        this->num_rows = this->ref_length;
        this->col_seq  = this->testSequence;
        this->num_cols = this->test_length;
    }
    this->initialized = true;
    this->left  = (this->num_cols < this->num_rows ? this->num_rows - this->num_cols : 1);
    this->right = 1;
    this->aligned = false;
    grow(false, false);
    return true;
}

int MAlignerEntry::align(){

    if( !this->initialized ) { return MINVAL; }

    while( !this->isAligned() ){
        this->step();
    }

    return getScore();

}

bool MAlignerEntry::step(){

    int num_rows = this->num_rows;
    int num_cols = this->num_cols;

    int *prevRow = NULL;
    int *thisRow = this->dpm->matrix;

    bool boundary = false;

    const int left  = this->left;
    const int right = this->right;

    int ii = 0;

    for(int y = 0; y < num_rows; ++y ){

        int start = (y >= left ? y - left : 0 );
        int stop  = (num_cols < y + right + 1 ? num_cols : y + right + 1);
        int lmax = MINVAL;
        int rmax = MINVAL;
        int line_max = MINVAL;

        int max_x, max_y;
        int rowOffset = 0;
       
        // This offset is necessary for the cases where the lefthand side of 
        // the search window falls outside of the DP matrix.
        // By doing this we can avoid any sort of zany offset garbage.
        if( start == 0 ){
            prevRow = prevRow - 1;
        }

        for( int x = start ; x < stop; ++x){ 
            assert( (ii > 0 && !(x == 0 && y == 0)) || ii == 0 );
            
            // add a matrix entry
            int t = this->scoreDP(x, y, rowOffset, prevRow, thisRow);

            // we need to know if we've run into a boundary
            if( x == start && x > 0){ lmax = t; }
            if( x == stop - 1 && x < num_cols - 1){ rmax = t; }
            if( t > line_max ){
                line_max = t; 
                max_x = x; 
                max_y = y; 
            }

            ++rowOffset;
            ++ii;
        }

        prevRow = thisRow;
        thisRow = thisRow + rowOffset; //TODO off by one?

        // if we've hit a boundary, short circuit our way out of the loop
        // and increase the boundary sizes
        if( num_cols > y + right + 1 ){
            bool grow_left = false;
            bool grow_right = false;

            if( line_max == lmax
                && this->left != this->num_cols / 2 ){
                boundary = true;
                grow_left = true;
            }

            if( line_max == rmax
                && this->right != this->num_cols / 2){
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

    int s = this->dpm->matrix[ii - 1];
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

    int leftVal = MINVAL;
    int upVal = MINVAL;
    int upleftVal = this->gap * x + this->gap * y;

    const char* row_seq = this->row_seq;
    const char* col_seq = this->col_seq;

    if(rowOffset > 0 ){
        leftVal = thisRow[rowOffset - 1] + this->gap;
    } else {
        leftVal = MINVAL;
    }

    if(y > 0 && x > 0){
        upleftVal = prevRow[rowOffset] + (row_seq[y] == col_seq[x] ? this->match : this->mismatch );
    } else {
        upleftVal = upleftVal + (row_seq[y] == col_seq[x] ? this->match : this->mismatch);
    }

    if((y > 0)
       && ( prevRow + rowOffset + 1 < thisRow )){
        upVal = prevRow[rowOffset + 1] + this->gap;
    } else {
        upVal = MINVAL;
    }

    thisRow[rowOffset] = max(upleftVal, max(leftVal, upVal));
    
    //printf("<%p, %p, %p> [x:%d y:%d] %d (%p)\n", upVal, leftVal, upleftVal, x, y, *(thisRow + rowOffset), (thisRow + rowOffset));
    return thisRow[rowOffset];
}

int main(int argc, char** argv){

    dp_matrix *dpm = (dp_matrix*) malloc(sizeof(dp_matrix));
    dpm->matrix = NULL;
    dpm->size = 0;

    MAlignerEntry *aligner = new MAlignerEntry(string("Sample").c_str(), string("ABCDEFGHIJKL").c_str(), dpm);

    aligner->initialize(string("ABCDEFGHIJKLMNOP").c_str());
    aligner->align();
    printf("%d\n", aligner->getScore());
    //aligner->initialize(string("XXXXXXXXXXXX").c_str());
    //aligner->align();
    //aligner->initialize(string("XXXXXXX").c_str());
    //aligner->align();
    //aligner->initialize(string("ABCDEFG").c_str());
    //aligner->align();

    delete aligner;

    return 0;
}
