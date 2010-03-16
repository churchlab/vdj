#include <assert.h>

#include <algorithm>
#include <queue>
#include <vector>
#include <string>

#ifndef __MAXIMUM_LIKELIHOOD_ALIGNER__
#define __MAXIMUM_LIKELIHOOD_ALIGNER__

double phredTable[] = {-0.6931471805599453,-0.10536051565782628,-1.005033585350145e-2,-1.0005003335835344e-3,-1.0000500033334732e-4,-1.0000050000287824e-5,-1.000000500029089e-6,-1.0000000494736474e-7,-1.0000000100247594e-8,-9.999999722180686e-10,-1.000000082790371e-10};

int binomialTable[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12, 13, 14, 15, 16, 17, 18, 18, 19, 20, 20, 21, 21, 22, 23, 23, 24, 24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 32, 32, 33, 33, 34, 34, 35, 35, 36, 36, 36, 37, 37, 38, 38, 39, 39, 40, 40, 41, 41, 41, 42, 42, 43, 43, 44, 44, 44, 45, 45, 46, 46, 47, 47, 47, 48, 48, 49, 49, 50, 50, 50, 51, 51, 52, 52, 52, 53, 53, 54, 54, 55, 55, 55, 56, 56, 57 };

class qread {
    public:
        std::string      _seq;
        std::vector<int> _qual;
};

class LLCell {
    public:
        LLCell();
        LLCell(double l, char b, int q);
        LLCell(char chrA, int qualA, char chrB, int qualB);
        void adjustLikelihood(double adj);
        char getBase();
        int  getQuality();
        double getLikelihood();
        bool isMatch();
    private:
        char _base;
        int _quality;
        double _likelihood;
        bool _match;
};
 
#endif
