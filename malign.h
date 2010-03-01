#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <iostream>
#include <fstream>
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

#ifndef MAXVAL
#define MAXVAL (16384)
#endif

typedef struct{
    int *matrix;
    int size;
} dp_matrix;

class MAlignerEntry {
            public:
                MAlignerEntry(std::string, std::string, dp_matrix *dp, int);
                ~MAlignerEntry();
                bool initialize(std::string, int);
                bool isAligned();
                bool isInitialized();
                int align();
                bool step();
                int upperBound();
                int lowerBound();
                int getScore();
                int getId(){ return _id; }
                int gap, match, mismatch;
                std::string getName();
            private:
                int _id;
                int scoreDP(int, int, int, int*, int*);
                std::string _name;
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
                bool doGrow;
                int setLowerBound(int, int, int);
                int setUpperBound(int, int, int);

                int score;
                bool aligned;
                bool initialized;
};

bool operator>(MAlignerEntry&, MAlignerEntry&);
bool operator<(MAlignerEntry&, MAlignerEntry&);

struct CompareMEntry : public std::binary_function<MAlignerEntry*, MAlignerEntry*, bool>{
    bool operator()(MAlignerEntry*, MAlignerEntry*) const;
};

typedef std::priority_queue<MAlignerEntry*, std::vector<MAlignerEntry*>, CompareMEntry> MAE_queue;

class MAlignerCore {
    public:
      
        MAlignerCore();
        ~MAlignerCore();
        void addEntry(std::string name, std::string sequence);
        std::string bestAlign(std::string input);
        std::string alignWith(std::string input, std::list<std::string> refs);
        MAE_queue align(std::string input);
        MAE_queue roundRobin(std::queue<MAlignerEntry*>);

    private:
        int _id;
        std::map<std::string,MAlignerEntry*> entries;
        dp_matrix *global_matrix;
};
   
