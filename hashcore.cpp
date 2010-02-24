#include "hashcore.h"

using namespace std;

unsigned short getNucleotide(char ch){
       switch(ch) {
        case 'a':
        case 'A': return A;
        case 't':
        case 'T': return T;
        case 'c':
        case 'C': return C;
        case 'g': 
        case 'G': return G;
        default:  return 5;
    }
}

OddsTable::OddsTable(int n){
    _size = n;
    _contains = (double*) malloc(sizeof(double) * (n - 1));
    _doesnot  = (double*) malloc(sizeof(double) * (n - 1));
    for( int c = 1; c < n ; c++ ){
        double cGivenF   = -1 * log(c + 1);
        double cGivenNF  = -1 * log(c * (n - c + 1));

        _contains[c-1] = cGivenF;
        _doesnot[c-1]  = cGivenNF;
        //TODO This is wrong, but it doesn't matter unless we actually want to
        // extract probabilities. This ought to be fixed.
        //printf("p(c|~f) = %f\n", cGivenNF);
        //printf("p(c|f) = %f\tp(c|~f)= %f (%d, %d)\n", cGivenF, cGivenNF, c, n);       
        //ncGivenF  = -1 * log((c + 1) * (n - c));
        //ncGivenNF = -1 * log(n - c + 1);
    }

    //printf("null score: %f\n", _featurelessOdds);
}

OddsTable::~OddsTable(){
    free(_contains);
    free(_doesnot);
}

double OddsTable::oddsGivenFeature(int c){

    assert( c > 0 && c <= _size );
    //printf("%d: %f\n", (c - 1), _contains[c-1]);
    return _contains[c-1];
}

double OddsTable::oddsWithoutFeature(int c){
    assert( c > 0 && c <= _size );
    return _doesnot[c - 1];
}

/*
 * FeatureSet
 */

string FeatureSet::getName() {
    return _name;
}

/*
string FeatureSet::getSequence() {
    return string(_sequence);
}
*/
/*
 * ReferenceSet
 */

ReferenceSet::ReferenceSet(list<ObservationSet*> observations ){

    _refs = new map<unsigned long, int>();
    _withOdds = new map<unsigned long, double>();
    //_withoutOdds = new map<unsigned long, double>();

    list<ObservationSet*>::iterator obs_itr;
    set<unsigned long>::iterator feat_itr;
    
    _size = observations.size();

    //TODO First get it working, then get it right. I can do this in linear time.
    //     I'm taking O(n lg n)
    for( obs_itr = observations.begin() ; obs_itr != observations.end() ; obs_itr++ ){
        set<unsigned long>* features = (*obs_itr)->getFeatureSet();
        
        for( feat_itr = features->begin() ; feat_itr != features->end() ; feat_itr++ ){
           insertBump(_refs, (*feat_itr)); 
        }
    }

    _odds = new OddsTable(_size);
    
    map<unsigned long, int>::iterator ref_itr, tmp_itr;
  
    /*
    // Throw out all universal features
    for( ref_itr = _refs->begin(); ref_itr != _refs->end() ; ){
        if( (*ref_itr).second == _size ){
            tmp_itr = ref_itr;
            ref_itr++;
            _refs->erase(tmp_itr);
        } else {
            ref_itr++;
        }
    }
*/
    _nullOdds = 0.0;

    for( ref_itr = _refs->begin(); ref_itr != _refs->end() ; ref_itr++ ){
        double with    = _odds->oddsGivenFeature((*ref_itr).second); 
        double without = _odds->oddsWithoutFeature((*ref_itr).second); 
        
        _withOdds->insert(pair<unsigned long, double>((*ref_itr).first, with + without));
        //_withoutOdds->insert(pair<unsigned long, double>((*ref_itr).first, without));
        _nullOdds += without;    
    }
}

double ReferenceSet::getNull(){ return _nullOdds; }

LikelihoodSet* ReferenceSet::makeLikelihood(ObservationSet* obs){
    //TODO This is a stopgap. This code needs to be cleaned up.
    return new LikelihoodSet(obs, this);
}

/*
 * ObservationSet: Raw observation
 */

//TODO I really don't like this type. Improve it.
ObservationSet::ObservationSet(string sequence, string name){

    _name = name;
    _sequence = sequence;
    
    _obsset = extractFeatures(sequence, new map<unsigned long, int>());
   // printf("Observation Set size: %d\n", _obsset->size());
}

//TODO call the supertype deconstructor
ObservationSet::~ObservationSet(){
    delete _obsset;
}

set<unsigned long>* ObservationSet::getFeatureSet() {
    map<unsigned long, int>::iterator f_itr;
    set<unsigned long> *features = new set<unsigned long>();

    for( f_itr = _obsset->begin(); f_itr != _obsset->end(); f_itr++ ){
        features->insert((*f_itr).first);    
    }

    return features;
}

/*
 * LikelihoodSet: Observations in log-likelihood space
 */
LikelihoodSet::LikelihoodSet(ObservationSet* obs, ReferenceSet* ref){
   
    map<unsigned long, double>::iterator ref_itr;
    set<unsigned long>::iterator obs_itr; 
    set<unsigned long> *obs_features = obs->getFeatureSet();

    _lset = new map<unsigned long, double>();

    //printf("obs length: %d\n", obs->_obsset->size());
    //TODO First get it working, then get it right.
    //TODO Make this robust -- handle keys not found
    for(obs_itr = obs_features->begin() ; obs_itr != obs_features->end(); obs_itr++){
       double val = (*ref->_withOdds->find(*obs_itr)).second;

       printf("%f\n", val);
       _lset->insert(pair<unsigned long, double>(*obs_itr, val));
    }
    _name = obs->getName();
    _nullOdds = ref->getNull();
    delete obs_features;
}

double LikelihoodSet::likelihood(ObservationSet* os){

    map<unsigned long, int>::iterator os_itr;
    map<unsigned long, double>::iterator lset_itr;
    double like = _nullOdds;
    for( os_itr = os->_obsset->begin() ; os_itr != os->_obsset->end(); os_itr++ ){
        lset_itr = _lset->find(os_itr->first);

        if( lset_itr != _lset->end() ){
            like += lset_itr->second;
        }
    }

    return like;
    
}

// This code was not generated by hand. Modifying it is ill-advised.
inline void runCombs( map<unsigned long, int>* hm, unsigned long xx, unsigned long mask, int ii ){

if(!((0x3B2D7 & mask)) && (ii >  18)) { insertBump(hm, ((0x38000 & xx) >> 6) | ((0x3000 & xx) >> 5) | ((0x200 & xx) >> 3) | ((0xc0 & xx) >> 2) | ((0x10 & xx) >> 1) | ((0x7 & xx) >> 0)); }
if(!((0x3C44D7 & mask)) && (ii >  22)) { insertBump(hm, ((0x3c0000 & xx) >> 10) | ((0x4000 & xx) >> 7) | ((0x400 & xx) >> 4) | ((0xc0 & xx) >> 2) | ((0x10 & xx) >> 1) | ((0x7 & xx) >> 0)); }
if(!((0xFFF & mask)) && (ii >  12)) { insertBump(hm, ((0xfff & xx) >> 0)); }
if(!((0x1A18AF & mask)) && (ii >  21)) { insertBump(hm, ((0x180000 & xx) >> 10) | ((0x20000 & xx) >> 9) | ((0x1800 & xx) >> 5) | ((0x80 & xx) >> 2) | ((0x20 & xx) >> 1) | ((0xf & xx) >> 0)); }
if(!((0x7747 & mask)) && (ii >  15)) { insertBump(hm, ((0x7000 & xx) >> 5) | ((0x700 & xx) >> 4) | ((0x40 & xx) >> 3) | ((0x7 & xx) >> 0)); }

    return;
}

inline void insertBump(map<unsigned long, int>* mp, unsigned long key){

    map<unsigned long, int>::iterator itr = mp->lower_bound(key);

    // Are we in the right place?
    if( (*itr).first == key ){
        int val = (*itr).second + 1;
        mp->erase(key);
        mp->insert(pair<unsigned long, int>(key, val));
    } else {
        // otherwise add it to the set
        mp->insert(itr, pair<unsigned long, int>(key, 1));
    }
    return;
}

map<unsigned long, int>* extractFeatures(string seq, map<unsigned long, int> *res = NULL){
    int len = seq.length();

    if( !res ){
        res = new map<unsigned long, int>();
    }

    unsigned long acc = 0;
    unsigned long nmask = 0;
    unsigned short nuc;
    
    int ii;
    for( ii = 0; ii < len; ++ii) {
        nuc = getNucleotide(seq[ii]);
        acc <<= 2;
        acc += nuc % 4;
        
        nmask <<= 2;
        if( nuc > G ){
            nmask += 3;
        }

        runCombs(res, acc, nmask, ii);    
    }
   
    return res;
}

map<unsigned long, int>* makeFeatureSet(list<char*> sequences){

    list<char*>::iterator seq_itr;
    map<unsigned long, int> *featureSet = new map<unsigned long, int>();

    for( seq_itr = sequences.begin(); seq_itr != sequences.end(); seq_itr++ ){
        extractFeatures(*seq_itr, featureSet);
    }

    return featureSet;
}

int main(int argc, char** argv){
    
    char* seq = "AAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAA";

    list<ObservationSet*> obsList;
    for(int ii = 0; ii < 10 ; ++ii ){
        obsList.push_back(new ObservationSet("Test", seq));    
    }

    ObservationSet *observations = new ObservationSet("Obs", seq);
    ReferenceSet *ref = new ReferenceSet(obsList);
    LikelihoodSet *ls = new LikelihoodSet(observations, ref);



    printf("Null: %f\n", ref->getNull());
    printf("Likelihood: %f\n", ls->likelihood(observations));
    
    list<char*> sequences;
    sequences.push_back(seq);

    map<unsigned long, int>* f = extractFeatures(seq);

    map<unsigned long, int>::iterator itr;

    for( itr = f->begin(); itr != f->end() ; itr++){
        printf("0x%X: %d\n", (unsigned int) (*itr).first, (*itr).second);
    }

    delete f;
    

    delete observations;
    return 0;

}
