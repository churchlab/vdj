#include <stdlib.h>
#include <queue>
using namespace std;

int binomialTable[] = { 0, 1, 2, 3, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 10, 11, 11, 12, 12, 13, 13, 13, 14, 14, 15, 15, 15, 16, 16, 16, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20, 21, 21, 22, 22, 22, 23, 23, 23, 24, 24, 24, 25, 25, 25, 26, 26, 26, 27, 27, 27, 28, 28, 28, 29, 29, 29, 30, 30, 30, 31, 31, 31 };

enum nucleotide { A, C, T, G };
enum quality { nil, Q10, Q20, Q30, Q40, Q50 };

nucleotide complement(nucleotide n){
    switch(n){
        case A: return T;
        case C: return G;
        case T: return A;
        case G: return G;
    }
}

typedef struct {
    nucleotide *_seq;
    quality    *_qual;
    int        length;
} qread;

qread align(qread readA, qread readB){
    int minLen = min(readA.length, readB.length);
    priority_queue<pair<double, int> > res;

    nucleotide *seqA = readA._seq;
    nucleotide *seqB = readB._seq;

    for(int ii = 10 ; ii < minLen; ii++){
        int hits = 0;
        for(int jj = 0; jj < ii; jj++){
            nucleotide baseA = seqA[readA.length - 1 - jj];
            nucleotide baseB = seqB[readB.length - 1 - jj];
            hits += (baseA == complement(baseB) ? 1 : 0);
        }

        if( hits > binomialTable[ii] ){
            double score = (double) hits / (double) ii;
            res.push(pair<double, int>(score, ii));
        }
    }

        qread result;

        if( res.empty() ){
            result.length = readA.length + readB.length + 1;
            result._seq = (nucleotide*) malloc(sizeof(nucleotide) * result.length);
            result._qual = (quality*) malloc(sizeof(quality) * result.length);
        } else {
            pair<double, int> so = res.top();

        }

        return result;

}
