#include "hasher.h"

using namespace std;
//

SequenceHasher::SequenceHasher(){

    _ref = NULL;
    _initialized = false;

}

SequenceHasher::~SequenceHasher(){
    if( _ref ){
        delete _ref;
    }

    list<ObservationSet*>::iterator obs_itr;

    for( obs_itr = _obs.begin(); obs_itr != _obs.end(); obs_itr++ ){
        delete (*obs_itr);
    }
    _obs.clear();
}

priority_queue<pair<double, string> > SequenceHasher::hash(string sequence){

    if( !_initialized ){
        initialize();
    }

    ObservationSet os(sequence);
    list<ObservationSet*>::iterator l_itr;

    priority_queue<pair<double, string> > rpq;
    for(l_itr = _obs.begin(); l_itr != _obs.end(); l_itr++){

        double likelihood = (*l_itr)->likelihood(&os, _ref);
        string nm = (*l_itr)->getName();
        rpq.push(pair<double, string>(likelihood, nm));
    }
    
    while(!rpq.empty()){
        double likelihood = rpq.top().first;
        string nm = rpq.top().second;
        
        printf("%s : %f\t", nm.c_str(), likelihood);
        rpq.pop();
    }

    printf("\n");
    return rpq;
}

void SequenceHasher::addReference(string name, string sequence){
    
    //TODO Ensure that arguments are in a consistent order
    ObservationSet* ob = new ObservationSet(sequence, name);
    _obs.push_back(ob);

    _initialized = false;
    return;
}

void SequenceHasher::initialize(){
    _initialized = true;
    if( !_ref ){ delete _ref; }

    _ref = new ReferenceSet(_obs);
    list<ObservationSet*>::iterator obs_itr;
    return;
}

int main(int argc, char** argv){
    
    SequenceHasher sq;
    sq.addReference(string("Poly-A"),string("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    sq.addReference(string("Poly-C"),string("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"));
    sq.addReference(string("Poly-G"),string("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"));
    /*sq.addReference(string("Poly-T"),string("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    sq.addReference(string("AthenC"),string("AAAAAAAAAAAAGAAAAAAAAAAAAAACCCCCCCCCCCCCCCCC"));
    sq.addReference(string("ATCG"),string("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"));
    */
    sq.initialize();
    sq.hash(string("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    sq.hash(string("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    sq.hash(string("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"));
    sq.hash(string("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"));
    sq.hash(string("CACACACACACACACACACACACACACACACACACACACACACACACA"));
    sq.hash(string("AAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCC"));

    return 0;
}

