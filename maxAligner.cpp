#include <assert.h>
#include <climits>
#include <cfloat>

#include <algorithm>
#include <queue>
#include <vector>
#include <string>

#include "maxAligner.h"

using namespace std;

LLCell::LLCell(){
    _base = 'x';
    _quality = 0;
    _likelihood = -5.0;
}
LLCell::LLCell(double l, char b, int q){
    _base = b;
    _quality = q;
    _likelihood = l;
    _match = false;
}
LLCell::LLCell(char chrA, int qualA, char chrB, int qualB){

    /* When faced with contradictory observations, prefer the more
     * likely base.
     * What score should we assign in the face of contradictory information?
     * Let B be some base, 
     * BQn = obs B with quality n
     * ~BQm = obs not-B wih quality m
     * P(B | BQn ^ ~BQm) = P(BQn | B) P(~BQm | B) / P(BQn ^ ~BQm)
     * P(BQn ^ ~BQm) = P(BQn) * P(~BQm) + (1 - P(BQn)) * (1 - P(~BQm))
     *
     * Unsurprisingly, this becomes vanishingly close to subtracting
     * Phred scores.
     */
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
    _likelihood = MLA::phred2log(_quality);
}
void LLCell::adjustLikelihood(double adj){ _likelihood += adj; }
char LLCell::getBase(){ return _base; }
int  LLCell::getQuality(){ return _quality; }
double LLCell::getLikelihood(){ return _likelihood; }
bool LLCell::isMatch(){ return _match; }

bool MLA::significantMatch(int matches, int overlap){
    //printf("Found %d matches with an overlap of %d (hurdle: %d)\n", matches, overlap, MLA::binomialTable[overlap]);
    return matches >= MLA::binomialTable[overlap];
}

double MLA::phredTable[] = {-0.6931471805599453,-0.10536051565782628,-1.005033585350145e-2,-1.0005003335835344e-3,-1.0000500033334732e-4,-1.0000050000287824e-5,-1.000000500029089e-6,-1.0000000494736474e-7,-1.0000000100247594e-8,-9.999999722180686e-10,-1.000000082790371e-10};
int MLA::binomialTable[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12, 13, 14, 15, 16, 17, 18, 18, 19, 20, 20, 21, 21, 22, 23, 23, 24, 24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 32, 32, 33, 33, 34, 34, 35, 35, 36, 36, 36, 37, 37, 38, 38, 39, 39, 40, 40, 41, 41, 41, 42, 42, 43, 43, 44, 44, 44, 45, 45, 46, 46, 47, 47, 47, 48, 48, 49, 49, 50, 50, 50, 51, 51, 52, 52, 52, 53, 53, 54, 54, 55, 55, 55, 56, 56, 57 };


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

double MLA::phred2log(int q){
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

    LLCell dpm[lenA+1][lenB+1];

    int ii;
    int jj;

    dpm[0][0] = LLCell(0, '*', 0);

    for( ii = 1 ; ii < lenA; ii++ ){
        dpm[ii][0] = LLCell(
                dpm[ii-1][0].getLikelihood() + MLA::phred2log(qualA[ii - 1]),
                seqA[ii - 1],
                qualA[ii - 1]);
    }

    dpm[lenA][0] = LLCell(
            seqA[lenA-1], qualA[lenA-1],
            seqB[0],qualB[0]
            );

    for( jj = 1 ; jj < lenB; jj++ ){
        dpm[0][jj] = LLCell(
                dpm[0][jj-1].getLikelihood() + MLA::phred2log(qualB[jj - 1]),
                seqB[jj - 1],
                qualB[jj - 1]);
    }

    dpm[0][lenB] = LLCell(
            seqA[0], qualA[0],
            seqB[lenB-1],qualB[lenB-1]);


    dpm[0][0] = LLCell(seqA[0],qualA[0],seqB[0],qualB[0]);
    
    for( jj = 1; jj < lenB; jj++ ){
        for( ii = 1; ii < lenA; ii++ ){
            dpm[ii][jj] = LLCell(seqA[ii - 1],qualA[ii - 1],seqB[jj - 1],qualB[jj - 1]);
            dpm[ii][jj].adjustLikelihood(dpm[ii-1][jj-1].getLikelihood());
        }
    }

    for( ii = 1 ; ii < lenA ; ii++ ){
        dpm[ii][lenB] = 
            LLCell( MLA::phred2log(qualA[ii]),
                    seqA[ii],
                    qualA[ii]);

        double bt = max(
                dpm[ii - 1][lenB].getLikelihood(),
                dpm[ii - 1][lenB -1].getLikelihood());

        dpm[lenA][jj].adjustLikelihood(bt);



        dpm[ii][lenB].adjustLikelihood( dpm[ii - 1][lenB].getLikelihood());

    }

    for( jj = 1; jj < lenB; jj++ ){
        dpm[lenA][jj] = 
            LLCell( MLA::phred2log(qualB[jj]),
                    seqB[jj],
                    qualB[jj]);

        double bt = max(
                dpm[lenA][jj - 1].getLikelihood(),
                dpm[lenA-1][jj - 1].getLikelihood());

        dpm[lenA][jj].adjustLikelihood(bt);

    }

    dpm[lenA][lenB] = LLCell(-DBL_MAX, '*', 0);

    ii = lenA;
    jj = lenB;

    string consensus;

    int matches = 0;
    int overlap = 0;

    vector<int> qualities;

    double upLikelihood, leftLikelihood, diagLikelihood;

    // Backtrace
    while( ii > 0 && jj > 0 ){

        upLikelihood   = -DBL_MAX;
        leftLikelihood = -DBL_MAX; 
        diagLikelihood = -DBL_MAX;

        if( 0 == ii || lenA == ii ){
            upLikelihood = dpm[ii][jj - 1].getLikelihood();
        }

        if( 0 == jj || lenB == jj ){
            leftLikelihood = dpm[ii - 1][jj].getLikelihood();
        }

        diagLikelihood = dpm[ii - 1][jj - 1].getLikelihood();

        if( diagLikelihood >= upLikelihood &&
                diagLikelihood >= leftLikelihood ){ 
            ii--;
            jj--;
            // Move diagonally through the overlap
            if( dpm[ii][jj].isMatch() ){
                matches++;
            }
            overlap++;


        } else if( upLikelihood >= leftLikelihood ){
            // Move up
            jj--;
        } else {
            // Move left
            ii--;
        }

        consensus += dpm[ii][jj].getBase();
        qualities.push_back(dpm[ii][jj].getQuality());
    }

    while(ii > 0){
        consensus += dpm[ii][jj].getBase();
        qualities.push_back(dpm[ii][jj].getQuality());
        ii--;
    }

    while(jj > 0){
        consensus += dpm[ii][jj].getBase();
        qualities.push_back(dpm[ii][jj].getQuality());
        jj--;
    }

    bool significant = MLA::significantMatch(matches, overlap);

    qread consensusRead;

    if( significant ){
        reverse(consensus.begin(), consensus.end());
        reverse(qualities.begin(), qualities.end());

        consensusRead._qual = qualities;
        consensusRead._seq  = consensus;
    } else {
        consensusRead._seq = seqA;
        consensusRead._seq += "-" + seqB;
        consensusRead._qual = qualA;
        consensusRead._qual.push_back(0);
        consensusRead._qual.insert(consensusRead._qual.end(), qualB.begin(), qualB.end());
    }

    return consensusRead;
}
/*
int main(int argc, char **argv){
    qread polyA, polyT;

    polyA._seq = string("AAAAAAAAAAA");
    polyA._qual = vector<int>(polyA._seq.size(),40);
    polyT._seq = string("GGGGGGGGTTTTG");
    polyT._qual = vector<int>(polyT._seq.size(),50);
    qread r = align(polyA, polyT);

    printf("%s\n", r._seq.c_str());
    vector<int>::iterator itr;

    for( itr = r._qual.begin() ; itr != r._qual.end(); itr++ ){
        printf("%d\t", *itr);
    }
    printf("\n");
    return 0;
} */
