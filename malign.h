#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <assert.h>

#include <algorithm>
#include <list>
#include <map>
#include <stack>
#include <string>
#include <queue>
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
                MAlignerEntry(const char *nm, const char *seq, dp_matrix *dp);
                ~MAlignerEntry();
                bool initialize(const char*, int);
                bool isAligned();
                bool isInitialized();
                int align();
                bool step();
                int upperBound();
                int lowerBound();
                int getScore();
                int gap, match, mismatch;

                const bool operator< (MAlignerEntry&);
                const bool operator> (MAlignerEntry&);
            private:

                int scoreDP(int, int, int, int*, int*);
                char *name;
                char *refSequence;
                int ref_length;

                char *testSequence;
                int test_length;

                // #rows >= #cols;
                char *row_seq;
                char *col_seq;
                int num_rows;
                int num_cols;

                dp_matrix *dpm;
                int matrix_size;
                
                int uBound;
                int lBound;
                int left;
                int right;

                int grow(bool, bool);
                int setLowerBound(int, int, int);
                int setUpperBound(int, int, int);

                int score;
                bool aligned;
                bool initialized;
};
 
class MAligner {
    public:
        MAligner(int nseq, char** seqs, char** names);
        ~MAligner();
        bool initialize(char* input);
        std::priority_queue<MAlignerEntry*> align(char* input);
        std::priority_queue<MAlignerEntry*> alignWith(char* input, int nnames, char** names);
        std::priority_queue<MAlignerEntry*> roundRobin(std::queue<MAlignerEntry*>);

    private:
        std::map<std::string,MAlignerEntry*> entries;
        dp_matrix *global_matrix;
};
   
