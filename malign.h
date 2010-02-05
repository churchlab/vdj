#include <list>

#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <utility>

#ifndef MINVAL
#define MINVAL (-16384)
#endif

typedef struct{
    int *matrix;
    int size;
} dp_matrix;

class MAlignerEntry {
            public:
                MAlignerEntry(char *nm, char *seq, dp_matrix *dp);
                bool initialize(char*, int);
                int align();
                bool step();
                int upperBound();
                int lowerBound();
                int gap, match, mismatch;
            private:

                int scoreDP(int, int, int, int*, int*);
                
                char *name;
                char *refSequence;
                char *testSequence;
                dp_matrix *dpm;
                int uBound;
                int lBound;
                int left;
                int right;
                int grow(bool, bool);
                int size;

                int max_rows;
                int max_cols;

                int score;
                bool aligned;
                bool initialized;
   

};
 
class MAligner {
    public:
        MAligner(int nseq, char** seqs, char** names);
        const char* align(char* input);
        const char* alignWith(char* input, int nnames, char** names);
    private:
        std::list<MAlignerEntry*> entries;
        dp_matrix *global_matrix;
};
   
