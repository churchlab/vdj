#include <assert.h>

#include <algorithm>
#include <queue>
#include <vector>
#include <string>

#ifndef __MAXIMUM_LIKELIHOOD_ALIGNER__
#define __MAXIMUM_LIKELIHOOD_ALIGNER__

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

class MLA {
    public:
        MLA();
        MLA(std::vector<int>);
        ~MLA();
        static double phred2log(int);
        bool significantMatch(int, int);
    private:
        static double phredTable[];
        static int defaultBinomialTable[];
        std::vector<int> binomialTable;

};

class LikelihoodAligner {
    public:
        LikelihoodAligner();
        LikelihoodAligner(std::vector<int>);
        ~LikelihoodAligner();
        qread align(qread, qread);
    private:
        MLA _table;
};
#endif
