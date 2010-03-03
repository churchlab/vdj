// Multi-aligner

#include "malign.h"

using namespace std;
/*
 * MAlignerCore
 */

MAlignerCore::MAlignerCore(){

    global_matrix = (dp_matrix*) malloc(sizeof(dp_matrix));
    global_matrix->matrix = NULL;
    global_matrix->size   = 0;
    global_matrix->_owner  = NULL;
    _id = 0;
}

void MAlignerCore::addEntry(string name, string sequence){

    MAlignerEntry *en = new MAlignerEntry(name, sequence, global_matrix, _id++);
    entries.insert(pair<string,MAlignerEntry*>(name,en));
}

MAlignerCore::~MAlignerCore(){
  
    map<string, MAlignerEntry*>::iterator m_itr;

    for( m_itr = entries.begin(); m_itr != entries.end() ; m_itr++){
        delete (*m_itr).second;
    }
   
    entries.clear();

    if( global_matrix->matrix != NULL ){
        free(global_matrix->matrix);
    }

    free(global_matrix);
}

MAE_queue MAlignerCore::align(string input){
    queue<MAlignerEntry*> robin; 
    map<string, MAlignerEntry*>::iterator entry_itr;

    int len = (int) input.length();

    for( entry_itr = entries.begin(); entry_itr != entries.end(); entry_itr++ ){
        MAlignerEntry *e = (*entry_itr).second;
        e->initialize(input, len);
        robin.push(e);
    }

    return this->roundRobin(robin);

}
/*
MAE_queue MAlignerCore::alignWith(string input, list<string> refs){

    queue<MAlignerEntry*> robin; 
    int len = (int) input.length();
    list<string>::iterator ref_itr;
    map<string,MAlignerEntry*>::iterator entry_itr;
    for( ref_itr = refs.begin() ; ref_itr != refs.end(); ref_itr++ ){
        entry_itr = entries.find(*ref_itr);
        if( entry_itr != entries.end() ){
            MAlignerEntry* mae = (*entry_itr).second;
            mae->initialize(input, len);
            robin.push(mae);
        }
    }

    MAE_queue res = this->roundRobin(robin);
    if( (int) res.size() > 0 ){
        return res.top().getName()
    }
}
*/

pair<string, string> MAlignerCore::bestAlign(string input){
    MAE_queue pq = align(input);
    if( !pq.empty() ){
        printf("Score from alignment: %d\n", pq.top()->getScore());
        return pq.top()->getAlignment();
    }

    return pair<string,string>(string(""), string(""));
}

MAE_queue MAlignerCore::roundRobin(queue<MAlignerEntry*> robin){

    MAE_queue results;
    int bestLowerBound = MINVAL;

    while( !robin.empty() ){
        MAlignerEntry *mae = robin.front();
        robin.pop();

        if( mae->upperBound() >= bestLowerBound ){
            mae->step();

            if( mae->lowerBound() > bestLowerBound ){
                bestLowerBound = mae->lowerBound();
            }

            if( mae->upperBound() >= bestLowerBound
                    && !mae->isAligned() ){
                robin.push(mae);
                mae->getAlignment();
            } 
        }

        if( mae->isAligned() ){
            results.push(mae); 
        }
    }
    
    return results;
}

/*
 * MAlignerEntry
 */
MAlignerEntry::MAlignerEntry(string name, string seq, dp_matrix *dp, int id){

    // copy the reference into place
    this->_id = id;
    this->_name = name;
    this->refSequence = (char*) malloc(seq.length() + 1);
    this->testSequence = NULL;
    strcpy(this->refSequence, seq.c_str());

    this->ref_length = strlen(this->refSequence);

    this->dpm = dp;
    this->num_rows = NULL;
    this->num_cols = NULL; // initialize this when we load it with a sequence

    this->uBound = MAXVAL;
    this->lBound = MINVAL;
    this->matrix_size = MINVAL;

    this->doGrow = false;

    this->gap      = -1;
    this->match    =  1;
    this->mismatch = -1;

    this->initialized = false;
    this->setLowerBound(0,0,0);
    this->setUpperBound(0,0,0);
}

MAlignerEntry::~MAlignerEntry(){
    free(this->refSequence);
    if( this->testSequence ){ free(this->testSequence); } 
}

bool operator> (MAlignerEntry& lhs, MAlignerEntry& rhs){
    return lhs.getScore() < rhs.getScore(); 
}
bool operator< (MAlignerEntry& lhs, MAlignerEntry& rhs){
    return lhs.getScore() < rhs.getScore();
}

bool CompareMEntry::operator()(MAlignerEntry* x, MAlignerEntry* y) const
{
    return *x < *y;
}


string MAlignerEntry::getName(){
    return this->_name;

}

bool MAlignerEntry::isAligned() {
    return aligned;
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
        this->_left = min(this->_left * 2, this->num_cols / 2);
    }

    if( growRight ){
        this->_right = min(this->_right * 2, this->num_cols / 2);
    }

    int l = this->_left;
    int r = this->_right;
    int d = this->num_rows - this->num_cols;    
    this->matrix_size = 
        this->num_rows 
        * (l + r + 1) 
        - (l * (l + 1) / 2) 
        - (r * (r + 1) / 2)
        - (r * d)
        - (d * (d+1) / 2); // account for _left hand overflow and the diagonal

    // Ask for a bigger scratch space if we need it
    if(this->doGrow && this->dpm->size < this->matrix_size ){ 
        free(this->dpm->matrix);
        this->dpm->matrix = (int*) malloc(sizeof(int) * this->matrix_size);
        this->dpm->size = this->matrix_size;
    }

    return this->matrix_size;
}

bool MAlignerEntry::initialize(string testSeq, int len=0){
    this->test_length = (int) testSeq.length();

    if( this->testSequence ){ free(this->testSequence); }
    //printf("Initialized entry with '%s' and length %d...\n", testSeq, len);
    this->testSequence = (char*) malloc(this->test_length + 1);
    strcpy(this->testSequence, testSeq.c_str());

    // test will be the row
    if( this->test_length > this->ref_length ) {
        this->ref_is_row  = false;
        this->row_seq  = this->testSequence;
        this->num_rows = this->test_length;
        this->col_seq  = this->refSequence;
        this->num_cols = this->ref_length;
        // reference will be the row
    } else {
        this->ref_is_row  = true;
        this->row_seq  = this->refSequence;
        this->num_rows = this->ref_length;
        this->col_seq  = this->testSequence;
        this->num_cols = this->test_length;
    }

    this->initialized = true;
    this->_left  = (this->num_cols < this->num_rows ? this->num_rows - this->num_cols : 1);
    this->_right = 1;
    this->aligned = false;
    this->doGrow = true;
    
    this->uBound = MAXVAL;
    this->lBound = MINVAL;

    grow(false, false);
    return true;
}

int MAlignerEntry::align(){

    if( !this->initialized ) { return MINVAL; }

    do{
        this->step();
    }while( !this->isAligned() );

    return getScore();

}

pair<string, string> MAlignerEntry::getAlignment(){

    if( !this->isAligned() || !this->myMemory() ){
        this->align();
    }

    int y = this->num_rows - 1;
    int x = this->num_cols - 1;
    const int _left  = this->_left;
    const int _right = this->_right;
    
    string rowString = "";
    string colString = "";

    assert((int) this->offsetData.size() == this->num_rows);
    int *prevRow = this->offsetData.front().first;
    int *currRow = this->offsetData.front().second;
    
    int *index = currRow + _left - ( this->num_rows - this->num_cols);
    int offset = index - currRow + 1;
    printf("Offset: %d Score: %d\n", offset, *index);
    while( x > 0 && y > 0 ){
        int* up   = prevRow + offset;
        int* upD  = prevRow + offset - 1;
        int* left = index - 1;

        assert(upD < currRow);

        int upScore   = MINVAL;
        int upDScore  = MINVAL;
        int leftScore = MINVAL;
        
        if( y > 0 && up < currRow ){
            upScore = *up;
        }

        if( x > 0 && left >= currRow ){
            leftScore = *left;
        }

        if( x > 0 && y > 0 ){
            upDScore = *upD;
        }
        //printf("left: %d diag: %d up: %d\n", left - currRow, upD - currRow, up - currRow);
        //printf("(%d,%d)[%ld]: left: %d diag: %d up: %d offset: %d here: %d\n", x, y,(index - this->dpm->matrix),leftScore, upDScore, upScore, offset, *index);

        if( upDScore >= upScore &&
            upDScore >= leftScore ){
            // either a match or mismatch
            rowString += row_seq[y];
            colString += col_seq[x];
            x--;
            y--;
            index = upD;
            offsetData.pop_front();
        } else if( leftScore > upScore ){
            // insert a row gap
            rowString += "-";
            colString += col_seq[x];
            x--;
            offset--;
            index = left;
        } else {
            // insert a column gap
            rowString += row_seq[y];
            colString += "-";
            y--;
            index = up;
            offset++;
            offsetData.pop_front();
        }

        prevRow = this->offsetData.front().first;
        currRow = this->offsetData.front().second;
    }
   
    if( x > 0 ){

        rowString += "-";
        colString += col_seq[x--];

    }

    if( y > 0 ){
        colString += "-";
        rowString += row_seq[y--];
    }

    printf("Score: %d\n", this->getScore());

    if( ref_is_row ){
        return pair<string,string>(rowString, colString);
    } else {
        return pair<string,string>(colString, rowString);
    }
}



bool MAlignerEntry::myMemory(){

    return ( this->dpm->_owner == this );

}

bool MAlignerEntry::step(){

    if( this->doGrow ){
        this->grow(false, false);
        this->doGrow = false;
    }

    int num_rows = this->num_rows;
    int num_cols = this->num_cols;

    int *prevRow = NULL;
    int *thisRow = this->dpm->matrix;

    // Take control of the memory so that we can output the backtrace.
    this->dpm->_owner = this;

    bool boundary = false;

    const int _left  = this->_left;
    const int _right = this->_right;

    int ii = 0;

    offsetData.clear();

    for(int y = 0; y < num_rows; ++y ){

        int start = (y >= _left ? y - _left : 0 );
        int stop  = (num_cols < y + _right + 1 ? num_cols : y + _right + 1);
        int lmax = MINVAL;
        int rmax = MINVAL;
        int line_max = MINVAL;

        int max_x = MINVAL;
        int max_y = MINVAL;
        int rowOffset = 0;

        // This offset is necessary for the cases where the lefthand side of 
        // the search window falls outside of the DP matrix.
        // By doing this we can avoid any sort of zany offset garbage.
        // This offset shift the previous row over by one, putting the up and
        // diagonal indices in the _right place.
        if( prevRow && start == 0 ){
            prevRow = prevRow - 1;
        }

        offsetData.push_front(pair<int*, int*>(prevRow, thisRow));

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
        if( num_cols > y + _right + 1 ){
            bool grow_left = false;
            bool grow_right = false;

            if( line_max == lmax
                    && this->_left != this->num_cols / 2 ){
                boundary = true;
                grow_left = true;
            }

            if( line_max == rmax
                    && this->_right != this->num_cols / 2){
                boundary = true;
                grow_right = true;
            }

            if( boundary ){ 
                this->doGrow = false;
                this->grow(grow_left, grow_right);
                this->doGrow = true;

                this->setUpperBound(max_x, max_y, line_max);
                this->setLowerBound(max_x, max_y, line_max);

                return MINVAL;
            }
        }
    }

    int s = this->dpm->matrix[ii - 1];
    _index = this->dpm->matrix + ii - 1;
    if( !boundary ){ 
        this->aligned = true; 
        this->score = s;
        this->uBound = s;
        this->lBound = s;
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

    return thisRow[rowOffset];
}
