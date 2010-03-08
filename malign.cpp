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

    map<string, BandedAligner*>::iterator entry_itr;
    for( entry_itr = entries.begin(); entry_itr != entries.end() ; entry_itr++ ){
        delete (*entry_itr).second;
    }
    entries.clear();

}

pair<string, pair<string, string> > MAlignerCore::align(string input, list<string> refs){
    queue<BandedAligner*> robin; 
    
    if( refs.empty() ){
        
        map<string, BandedAligner*>::iterator entry_itr;
        
        for( entry_itr = entries.begin(); entry_itr != entries.end(); entry_itr++ ){
            BandedAligner *e = (*entry_itr).second;
            e->initialize(input);
            robin.push(e);
        }

    } else {

        list<string>::iterator sel_itr;
        map<string, BandedAligner*>::iterator ref_itr;

        // extract all our keys
        for( sel_itr = refs.begin() ; sel_itr != refs.end() ; sel_itr++ ){
            ref_itr = entries.find(*sel_itr);
            if( ref_itr != entries.end() ){
                BandedAligner *e = (*ref_itr).second;
                e->initialize(input);
                robin.push(e);
            }
        }
    }

    multimap<int, BandedAligner*> algnResult = roundRobin(robin);


    if( !algnResult.empty() ){
        multimap<int, BandedAligner*>::reverse_iterator res_itr;
        res_itr = algnResult.rbegin(); // get the last element
        BandedAligner *bestAlignment = (*res_itr).second;
        pair<string, string> bt = bestAlignment->getBacktrace();
        return pair<string, pair<string, string> >(bestAlignment->getName(), bt);
    }

    pair<string,string> dummy(string(""),string("")); 
    return pair<string, pair<string, string> >(string(""), dummy);
}

multimap<int, BandedAligner*> MAlignerCore::roundRobin(queue<BandedAligner*> robin){

    multimap<int, BandedAligner* > results;
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
            results.insert(pair<int, BandedAligner*>(algn->getScore(), algn)); 
        }
    }
    
    return results;
}
/*
int main() {

    MAlignerCore a;
    a.addEntry("polyATCG", "ATCGATCG");
    pair<string, pair<string, string> > alignment = a.align("GATCGAT");
    string name = alignment.first;
    pair<string, string> bt = alignment.second;

    printf("%s\n\t%s\n\t%s\n", name.c_str(), bt.first.c_str(), bt.second.c_str());
    return 0;

}
*/
