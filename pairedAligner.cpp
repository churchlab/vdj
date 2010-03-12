#include <assert.h>

#include <algorithm>
#include <queue>
#include <vector>
#include <string>

using namespace std;

double phredTable[] = {-0.6931471805599453,-0.10536051565782628,-1.005033585350145e-2,-1.0005003335835344e-3,-1.0000500033334732e-4,-1.0000050000287824e-5,-1.000000500029089e-6,-1.0000000494736474e-7,-1.0000000100247594e-8,-9.999999722180686e-10,-1.000000082790371e-10};

class LLCell {
    public:
        LLCell(){
            _base = '-';
            _quality = 0;
            _likelihood = 0.0;
        }
        LLCell(double l, char b, int q){
            _base = b;
            _quality = q;
            _likelihood = l;
            _match = false;
        }
        LLCell(char chrA, int qualA, char chrB, int qualB){
            _match = false;
            if( chrA == chrB ){
                _quality = qualA + qualB;
                _base = chrA;
                _match = true;
            } else if( qualA >= qualB ){
                _quality = qualA - qualB;
                _base = chrA;
            } else {
                _quality = qualB - qualA;
                _base = chrB;
            }
            _likelihood = phredTable[_quality/10];
        }
        void adjustLikelihood(double adj){ _quality += adj; }
        char getBase(){ return _base; }
        int  getQuality(){ return _quality; }
        double getLikelihood(){ return _likelihood; }
        bool isMatch(){ return _match; }
    private:
        char _base;
        int _quality;
        double _likelihood;
        bool _match;

};

class qread {
    public:
        std::string      _seq;
        std::vector<int> _qual;
};


    
    int binomialTable[] = { 0, 1, 2, 3, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 10, 11, 11, 12, 12, 13, 13, 13, 14, 14, 15, 15, 15, 16, 16, 16, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20, 21, 21, 22, 22, 22, 23, 23, 23, 24, 24, 24, 25, 25, 25, 26, 26, 26, 27, 27, 27, 28, 28, 28, 29, 29, 29, 30, 30, 30, 31, 31, 31 };

bool significantMatch(int matches, int overlap){
    return matches > binomialTable[overlap];
}

bool goodBase(char n){
    switch(n){
        case 'a':
        case 'A': 
        case 'c':
        case 'C':
        case 't':
        case 'T':
        case 'g':
        case 'G': return true;
        default:  return false;
    }
}

char normalizeBase(char n){
    switch(n){
        case 'a':
        case 'A': return 'A';
        case 'c':
        case 'C': return 'C';
        case 't':
        case 'T': return 'T';
        case 'g':
        case 'G': return 'G';
        default:  return n;
    }
}

char complement(char n){
    switch(n){
        case 'a': return 't';
        case 'A': return 'T';
        case 'c': return 'g';
        case 'C': return 'G';
        case 't': return 'a';
        case 'T': return 'A';
        case 'g': return 'c';
case 'G': return 'C';
        default:  return n;
    }
}

string reverse_complement(string s){
    string res;
    res.resize(s.size());
    reverse(s.begin(), s.end());
    transform(s.begin(),s.end(), res.begin(), complement);
    return res;
}

double phred2log(int q){
    return phredTable[q/10];

}

qread align(qread readA, qread readB){

    assert(readA._qual.size() == readA._seq.size());
    assert(readB._qual.size() == readB._seq.size());

    string seqA = readA._seq;
    string seqB = reverse_complement(readB._seq);
    vector<int> qualA = readA._qual;
    vector<int> qualB = readB._qual;
    reverse(qualB.begin(), qualB.end());

    int lenA = seqA.size();
    int lenB = seqB.size();

    LLCell dpm[lenA + 1][lenB + 1];

    dpm[0][0] = LLCell(
            phred2log(qualA[0]),
            seqA[0],
            qualA[0]);

    int ii, jj;
    for(ii = 1; ii < lenA; ii++){
       dpm[ii][0] =
           LLCell(
               dpm[ii - 1][0].getLikelihood() + phred2log(qualA[ii]),
                seqA[ii],
                qualA[ii]);
    }

    dpm[ii][0] =
        LLCell(
            dpm[ii - 1][0].getLikelihood() + phred2log(0),
            '-',
            0);

    for(jj = 1; jj < lenB; jj++){
        for(ii = 0; ii < lenA - 1; ii++){
            double likelihood;
            char chrA = seqA[ii];
            char chrB = seqB[jj];
            char consensusBase;
            int qual;
            dpm[ii][jj] = LLCell(
                    chrA,
                    qualA[ii],
                    chrB,
                    qualB[jj]);
            dpm[ii][jj].adjustLikelihood(dpm[ii - 1][jj - 1].getLikelihood());
        }
    }

    for(ii = 1; ii < lenB; ii++){
        dpm[lenA][ii] =
            LLCell(
                dpm[lenA][ii - 1].getLikelihood() + phred2log(qualB[ii]),
                seqB[ii],
                qualB[ii]);
    }


    ii = lenB;
    jj = lenA;

    string consensus;
 
    int matches = 0;
    int overlap = 0;

    printf("Log-Likelihood: %f\n", dpm[ii][jj].getLikelihood());
    while( ii > 0 && jj > 0 ){

        double upLikelihood   = dpm[ii][jj - 1].getLikelihood();
        double leftLikelihood = dpm[ii - 1][jj].getLikelihood();
        double diagLikelihood = dpm[ii - 1][jj - 1].getLikelihood();

        consensus += dpm[ii][jj].getBase();
        if( diagLikelihood >= upLikelihood &&
            diagLikelihood >= leftLikelihood ){ 
            
            if( dpm[ii][jj].isMatch() ){
                matches++;
            }
            overlap++;

            ii--;
            jj--;
        } else if( upLikelihood >= leftLikelihood ){
            jj--;
        } else {
            ii--;
        }
    }

    while(ii > 0){
        consensus += dpm[ii][jj].getBase();
        ii--;
    }

    while(jj > 0){
        consensus += dpm[ii][jj].getBase();
        jj--;
    }

    printf("%s\n", (significantMatch(overlap, matches) ? "SIGNIFICANT" : "NOT SIGNIFICANT"));
    reverse(consensus.begin(), consensus.end());
    return readA;

}

int main(int argc, char **argv){
    qread polyA, polyT;

    
    polyA._seq = string("AAACCCTTTT");
    polyA._qual = vector<int>(polyA._seq.size(),50);
    polyT._seq = string("GGGAAAAGG");
    polyT._qual = vector<int>(polyT._seq.size(),50);
    qread r = align(polyA, polyT);
    printf("%s\n", r._seq.c_str());
    return 0;
}
