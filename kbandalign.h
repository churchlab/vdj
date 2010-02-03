#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <string>
#include <algorithm>
#include <utility>
#include <list>
#include <queue>

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

#ifndef MAX_INT
#define MAX_INT (16384)
#endif

class Alignment {
    public:
        Alignment(std::string seqA, std::string seqB);
        ~Alignment();
        std::string seqA();
        std::string seqB();
        int upperBound();
        int lowerBound();
        int globalScore();
        int align();
        int step();
        bool canStep();
        bool isAligned();
        static Alignment* round_robin(std::list<std::string>, std::string);
    private:
        int getScore(int x,  int y,  int index);
        int resize(int, int);
        int setUpperBound(int x, int y, int score = 0);
        int setLowerBound(int x, int y, int score = 0);
        int size;
        int *dp_matrix;
        int left;
        int right;
        char* row_seq;
        char* col_seq;
        int max_rows;
        int max_cols;
        
        int uBound;
        int lBound;
        int score;

        bool aligned;
};

inline int score(int*, int,  int,  int,  int,  char,  char, bool);
