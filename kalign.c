// k-band alignment
//   Sequence alignment using a diagonal band.
//   Where n,m are inputs and e is the edit distance:
//     O(e^2) + O(n) + O(m)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#ifndef GAP
#define GAP (-1)
#endif

#ifndef MATCH
#define MATCH (1)
#endif

#ifndef MISMATCH
#define MISMATCH (-1)
#endif

#ifndef MIN_INT
#define MIN_INT (-16384)
#endif

inline int score(
    short int *dp_matrix,
    unsigned short int x,
    unsigned short int y,
    unsigned short int index,
    unsigned short int max_col,
    unsigned char L,
    unsigned char R,
    bool match){

    short int diag, top, left;
    short int di, ti, li;
    unsigned short int overflow  = 0;
    unsigned short int underflow = 0;

    // Get the overflow for the previous row
    if( R + y - max_col > 0 ){
        overflow = R + y - max_col;
    }

    // How much did the previous row underflow by?
    if( y - L < 1 ){
        underflow = L - y;
    }

    // If the cell to our left exists, get the value
    if( y - x < L ) {
        li = index - 1;
        left = dp_matrix[li] + GAP;
    } else {
        li = -1;
        left = MIN_INT;
    }
    
    // Do the same with the cell above us
    if( x - y < R && y != 0 ){
        ti = index - (L - underflow) - (R - overflow);
        top = dp_matrix[ti] + GAP;
    } else {
        ti = -1;
        top = MIN_INT;
    }

    // and finally the diagonal
    if( y > 0 && x > 0 ){
        di = index - (L - underflow) - (R - overflow) - 1;
        diag = dp_matrix[di];
        diag += ( match ? MATCH : MISMATCH );
    } else {
        diag = GAP * x + GAP * y + ( match ? MATCH : MISMATCH );
        di = -1;
    }
   // take the max of the three neighbors
    dp_matrix[index] = (diag > top ? (diag > left ? diag : left) : (top > left ? top : left));

    // printf("x:%3d y:%3d [%3d] %3d %3d (%3d %3d %3d %2d)\n", x, y, index, diag, dp_matrix[index], di, ti, li, underflow);
    return dp_matrix[index];
}

int align( char* seqA, char* seqB ){

    int max_rows = strlen(seqA);
    int max_cols = strlen(seqB);
    
    // Let's avoid dealing with twice as many edge cases:
    if( max_cols > max_rows ){
        return align(seqB, seqA);
    }

    int right = 1;
    int left  = 1;
    //ensure that we reach the bottom corner.
    //TODO should I guarantee that this is a power of 2?
    if( max_cols < max_rows ){
        left = max_rows - max_cols; 
    }
    int d = max_rows - max_cols;    
   int ii; // index variable
   while( left < max_cols / 2 && right < max_cols / 2 ){
        ii = 0; // reset the index
        int size = max_rows 
                   * (left + right + 1) 
                   - (left * (left + 1) / 2) 
                   - (right * (right + 1) / 2)
                   - (right * d)
                   - (d * (d+1) / 2); // account for left hand overflow and the diagonal
    /*
        printf("Rows: %d\tCols: %d\tLeft: %d\tRight: %d\n", max_rows, max_cols, left, right);
        printf("%d - %d - %d - %d - %d = %d\n", (max_rows * (left + right + 1)), (left * (left + 1) / 2), (right * (right + 1) / 2), (right * d), (d * (d+1) / 2), size); 
      */
        bool boundary = false;
        short int *dp_matrix = (short int*) malloc(sizeof(short int) * size);


        for(int y = 0; y < max_rows; ++y ){
            int start = (y >= left ? y - left : 0 );
            int stop  = (max_cols < y + right + 1 ? max_cols : y + right + 1);
            int lmax = MIN_INT;
            int rmax = MIN_INT;
            int line_max = MIN_INT;
            
            for( int x = start ; x < stop; ++x){ 
                // add a matrix entry
                int t = score(dp_matrix, x, y, ii, max_cols, left, right, seqB[x] == seqA[y]);
                    
                // we need to know if we've run into a boundary
                if( x == start && x > 0){ lmax = t; }
                if( x == stop - 1 && x < max_cols - 1){ rmax = t; }
                if( t > line_max ){ line_max = t; }

                ++ii; // advance the index
            }
            
            // if we've hit a boundary, short circuit our way out of the loop
            // and increase the boundary sizes
            // unless we're in the end game
            if( max_cols > y + right + 1 ){
                if( line_max == lmax){ left  *= 2; boundary = true; }
                if( line_max == rmax ){ right *= 2; boundary = true; }
                if( boundary ){ 
                    printf("Hit the boundary...(%d < %d - %d)\n", y, max_rows, right); 
                    break;
                }
            }
        }

    
        int s = dp_matrix[size - 1];
        printf("%d\t%d\n", ii, size);
        printf("Rows: %d\tCols: %d\tLeft: %d\tRight: %d\n", max_rows, max_cols, left, right);
        printf("%d - %d - %d - %d - %d = %d\n", (max_rows * (left + right + 1)), (left * (left + 1) / 2), (right * (right + 1) / 2), (right * d), (d * (d+1) / 2), size); 
        printf("Size: %d\tScore: %d\tLeft: %d\tRight: %d\n", size, s, left, right);
        assert(boundary || ii == size);
        free(dp_matrix);
        
        if( !boundary ){ return s; }
    }
    
    return MIN_INT;
}

int main( int argc, void** argv ){

    printf("%d\n", align("aaxaaxxaaaaaa", "aaaaaaaxaaa"));
    return 0;

}
