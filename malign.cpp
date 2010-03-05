// Multi-aligner

#include "malign.h"

using namespace std;
/*
 * MAlignerCore
 */

MAlignerCore::MAlignerCore(){

}

void MAlignerCore::addEntry(string name, string sequence){
    entries.insert(pair<string, BandedAligner*>(name, new BandedAligner(name, sequence)));

}

MAlignerCore::~MAlignerCore(){
  
    entries.clear();

}

list<pair<int, BandedAligner*> > MAlignerCore::align(string input){
    queue<BandedAligner*> robin; 
    map<string, BandedAligner*>::iterator entry_itr;

    list<pair<int, BandedAligner*> > res;

    for( entry_itr = entries.begin(); entry_itr != entries.end(); entry_itr++ ){
        BandedAligner *e = (*entry_itr).second;
        e->initialize(input);
        e->align();
        e->getBacktrace();
        res.push_back(pair<int, BandedAligner*>(e->getScore(), e));
    }
    return res;

}
/*
MAE_queue MAlignerCore::alignWith(string input, list<string> refs){

    queue<MAlignerEntry*> robin; 
    int len = (int) input.length();
    list<string>::iterator ref_itr;
    map<string,MAlignerEntry*>::iterator entry_itr;
    for( ref_itr = refs.begin() ; ref_itr != refs.end(); ref_itr++ ){
        entry_itr = entries.find(*ref_itr);
        if( entry_itr != entries.end() ){
            MAlignerEntry* algn = (*entry_itr).second;
            algn->initialize(input, len);
            robin.push(algn);
        }
    }

    MAE_queue res = this->roundRobin(robin);
    if( (int) res.size() > 0 ){
        return res.top().getName()
    }
}
*/

pair<string, pair<string, string> > MAlignerCore::bestAlign(string input){
    list<pair<int, BandedAligner*> > res = align(input);

    int bestScore = -8388608;
    BandedAligner* bestAlignment = NULL;

    list<pair<int, BandedAligner*> >::iterator itr;

    for( itr = res.begin(); itr != res.end(); itr++ ){
        if((*itr).first > bestScore){
            bestScore = (*itr).first;
            bestAlignment = (*itr).second;
        }
    }

    if( bestAlignment ){
        pair<string, string> bt = bestAlignment->getBacktrace();
        return pair<string, pair<string, string> >(bestAlignment->getName(), bt);
    }

    pair<string,string> dummy(string(""),string("")); 
    return pair<string, pair<string, string> >(string(""), dummy);
}

list<pair<int, BandedAligner*> > MAlignerCore::roundRobin(queue<BandedAligner*> robin){

    list<pair<int, BandedAligner*> > results;
    int bestLowerBound = -8388608;
    int bestScore = bestLowerBound;
    BandedAligner* bestAlignment = NULL;

    while( !robin.empty() ){
        BandedAligner *algn = robin.front();
        robin.pop();

        if( algn->upperBound() >= bestLowerBound ){
            algn->step();

            if( algn->lowerBound() > bestLowerBound ){
                bestLowerBound = algn->lowerBound();
            }

            if( algn->upperBound() >= bestLowerBound
                    && !algn->isAligned() ){
                robin.push(algn);
            } 
        }

        if( algn->isAligned() ){
            if( algn->getScore() > bestScore ){
                bestScore = algn->getScore();
                bestAlignment = algn;
            }
            results.push_back(pair<int, BandedAligner*>(algn->getScore(), algn)); 
        }
    }
    printf("Got here..\n");
    if( bestAlignment ){
        pair<string,string> bt = bestAlignment->getBacktrace();
    }
    return results;
}
