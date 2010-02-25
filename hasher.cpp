#include "hasher.h"

using namespace std;

SequenceHasher::SequenceHasher(){

    _ref = NULL;
    _initialized = false;

}

SequenceHasher::~SequenceHasher(){
    if( _ref ){
        delete _ref;
    }

    list<ObservationSet*>::iterator obs_itr;
    list<LikelihoodSet*>::iterator like_itr;

    for( obs_itr = _obs.begin(); obs_itr != _obs.end(); obs_itr++ ){
        delete (*obs_itr);
    }
    _obs.clear();

    for( like_itr = _likelihoods.begin(); like_itr != _likelihoods.end(); like_itr++ ){
        delete (*like_itr);
    }
    _likelihoods.clear();

}

void SequenceHasher::hash(string sequence){

    if( !_initialized ){
        initialize();
    }

    ObservationSet os(sequence);
    list<LikelihoodSet*>::iterator l_itr;

    priority_queue<pair<double, string> > rpq;
    printf("%s: ", sequence.c_str());
    printf("[ ");
    for(l_itr = _likelihoods.begin(); l_itr != _likelihoods.end(); l_itr++){

        double likelihood = (*l_itr)->likelihood(&os);
        string nm = (*l_itr)->getName();
        rpq.push(pair<double, string>(likelihood, nm));
    }
    
    while(!rpq.empty()){
        double likelihood = rpq.top().first;
        string nm = rpq.top().second;
        printf("'%s' : %f\t", nm.c_str(), likelihood);
        rpq.pop();
    }
    printf("]\n");

    return;
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

    list<LikelihoodSet*>::iterator like_itr;

    for( like_itr = _likelihoods.begin(); like_itr != _likelihoods.end(); like_itr++ ){
        delete (*like_itr);
    }
    _likelihoods.clear();
    
    if( !_ref ){ delete _ref; }

    _ref = new ReferenceSet(_obs);
    list<ObservationSet*>::iterator obs_itr;
    for( obs_itr = _obs.begin(); obs_itr != _obs.end(); obs_itr++ ){
        LikelihoodSet *ls = new LikelihoodSet(*obs_itr, _ref);
        _likelihoods.push_back(ls);
    }
    return;
}
/*
int main(int argc, char** argv){
    
    SequenceHasher sq;
    sq.addReference(string("Poly-A"),string("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    sq.addReference(string("Poly-C"),string("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"));
    sq.addReference(string("Poly-G"),string("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"));
    sq.addReference(string("Poly-T"),string("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    sq.addReference(string("AthenC"),string("AAAAAAAAAAAAGAAAAAAAAAAAAAACCCCCCCCCCCCCCCCC"));
    sq.addReference(string("ATCG"),string("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"));
    sq.initialize();
    sq.hash(string("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    sq.hash(string("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"));
    sq.hash(string("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"));
    sq.hash(string("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"));
    sq.hash(string("CACACACACACACACACACACACACACACACACACACACACACACACA"));
    sq.hash(string("AAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCC"));

    return 0;
}
*/
