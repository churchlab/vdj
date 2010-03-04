// Multi-aligner

#include "malign.h"

using namespace std;
/*
 * MAlignerCore
 */

MAlignerCore::MAlignerCore(){

}

void MAlignerCore::addEntry(string name, string sequence){

}

MAlignerCore::~MAlignerCore(){
  
    entries.clear();

}

list<BandedAligner*> MAlignerCore::align(string input){
    queue<BandedAligner*> robin; 
    map<string, BandedAligner*>::iterator entry_itr;

    int len = (int) input.length();

    for( entry_itr = entries.begin(); entry_itr != entries.end(); entry_itr++ ){
        BandedAligner *e = (*entry_itr).second;
        e->initialize(input);
        robin.push(e);
    }

    return this->roundRobin(robin);

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

pair<string, string> MAlignerCore::bestAlign(string input){
    return pair<string,string>(string(""), string(""));
}

list<BandedAligner*> MAlignerCore::roundRobin(queue<BandedAligner*> robin){

    list<BandedAligner*> results;
    int bestLowerBound = MINVAL;

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
                algn->getAlignment();
            } 
        }

        if( algn->isAligned() ){
            results.push_back(algn); 
        }
    }
    
    return results;
}
