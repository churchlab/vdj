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

double MLA::phredTable[] = {-2.218712235567753, -1.5814737534084538,-0.9968430440078468,-0.6955244713323139,-0.5076758736967449,-0.3801304080661717,-0.2892681872016077,-0.22255151597283288,-0.17255657291383555,-0.13455196028877356,-0.10536051565782628,-8.276530266916031e-2,-6.517417319943086e-2,-5.1418274157962233e-2,-4.062484422164181e-2,-3.2133574023085816e-2,-2.5439727534223144e-2,-2.0154364761232325e-2,-1.5975869246699935e-2,-1.2669170208635413e-2,-1.005033585350145e-2,-7.974998278512694e-3,-6.329562931111249e-3,-5.024473890985146e-3,-3.989017266406581e-3,-3.167288226157358e-3,-2.5150465111819995e-3,-1.9972555025509443e-3,-1.586150464280187e-3,-1.2597185241064498e-3,-1.0005003335835344e-3,-7.946438805584922e-4,-6.311564818346209e-4,-5.013128699287788e-4,-3.9818643625136844e-4,-3.1627777656023393e-4,-2.51220196302187e-4,-1.9954613950356738e-4,-1.5850188000545264e-4,-1.2590046631056348e-4,-1.0000500033334732e-4,-7.943597842624206e-5,-6.309772506765721e-5,-5.011997934787248e-5,-3.981150952301729e-5,-3.162327661227747e-5,-2.5119179799065442e-5,-1.9952822205921343e-5,-1.5849057520282164e-5,-1.2589333363280912e-5,-1.0000050000287824e-5,-7.94331389524571e-6,-6.309593350271555e-6,-5.011884895732116e-6,-3.981079630003642e-6,-3.16228266018417e-6,-2.51188958631412e-6,-1.995264305462918e-6,-1.5848944484078974e-6,-1.258926204253983e-6,-1.000000500029089e-6,-7.94328550222186e-7,-6.309575435365988e-7,-5.011873591979585e-7,-3.981072498074477e-7,-3.162278159729645e-7,-2.5118867465567504e-7,-1.9952625144221186e-7,-1.584893317835645e-7,-1.2589254909403657e-7,-1.0000000494736474e-7,-7.943282658471264e-8,-6.309573640552222e-8,-5.0118724649609444e-8,-3.981071784872332e-8,-3.162277711949987e-8,-2.5118864637361162e-8,-1.9952623308329733e-8,-1.5848932070211884e-8,-1.2589254162894963e-8,-1.0000000100247594e-8,-7.943282396744903e-9,-6.309573439199526e-9,-5.0118723741305115e-9,-3.981071702449905e-9,-3.162277666949986e-9,-2.5118864797519564e-9,-1.99526229071369e-9,-1.58489321792216e-9,-1.258925409157477e-9,-9.999999722180686e-10,-7.94328270141873e-10,-6.309573976396216e-10,-5.011872029760115e-10,-3.9810721394070965e-10,-3.1622771073384733e-10,-2.511886254868043e-10,-1.9952628440337202e-10,-1.5848933278141135e-10,-1.25892518639967e-10,-1.000000082790371e-10}; 
    
    -0.6931471805599453,-0.10536051565782628,-1.005033585350145e-2,-1.0005003335835344e-3,-1.0000500033334732e-4,-1.0000050000287824e-5,-1.000000500029089e-6,-1.0000000494736474e-7,-1.0000000100247594e-8,-9.999999722180686e-10,-1.000000082790371e-10};
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
    return phredTable[q];

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
