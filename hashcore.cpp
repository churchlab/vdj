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
        double cd = c;
        double nd = n;

        double p = 1 / (cd+1);
        double np = 1 / (cd * (nd - cd + 1));

        double cGivenF   = log(p  / (1 - p));
        double cGivenNF  = log(np / (1 - np));

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
    return _contains[c-1];
}

double OddsTable::oddsWithoutFeature(int c){
    assert( c > 0 && c <= _size );
    return _doesnot[c - 1];
}

/*
 * ReferenceSet
 */

ReferenceSet::ReferenceSet(list<ObservationSet*> observations ){

    _refs = new map<unsigned long, int>();
    _withOdds = new map<unsigned long, double>();
    _withoutOdds = new map<unsigned long, double>();

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

    for( ref_itr = _refs->begin(); ref_itr != _refs->end() ; ref_itr++ ){
        double with    = _odds->oddsGivenFeature((*ref_itr).second); 
        double without = _odds->oddsWithoutFeature((*ref_itr).second); 
        
        _withOdds->insert(pair<unsigned long, double>((*ref_itr).first, with));
        _withoutOdds->insert(pair<unsigned long, double>((*ref_itr).first, without));
    }
}

double ReferenceSet::getScore(unsigned long key, bool present){
    
    map<unsigned long, double>::iterator odds_itr;
    double score = 0.0;

    

    if( present ){
        odds_itr = _withOdds->find(key);
        if( odds_itr != _withOdds->end() ){
            score = (*odds_itr).second;
        }
    } else {
         odds_itr = _withoutOdds->find(key);
        if( odds_itr != _withoutOdds->end() ){
            score = (*odds_itr).second;
        }
    }

    return score;
}

/*
 * ObservationSet: Raw observation
 */
ObservationSet::ObservationSet(string sequence, string name){

    _name = name;
    _sequence = sequence;
    _obsset = extractFeatures(sequence);
   // printf("Observation Set size: %d\n", _obsset->size());
}

ObservationSet::~ObservationSet(){
    delete _obsset;
}

string ObservationSet::getName() {
    return _name;
}

set<unsigned long>* ObservationSet::getFeatureSet() {
    return new set<unsigned long>(*_obsset);
}

double ObservationSet::likelihood(ObservationSet* os, ReferenceSet* ref){

    set<unsigned long>::iterator t_itr = os->_obsset->begin();
    set<unsigned long>::iterator t_itr_end = os->_obsset->end();

    set<unsigned long>::iterator r_itr     = this->_obsset->begin();
    set<unsigned long>::iterator r_itr_end = this->_obsset->end();

    set<unsigned long> intersection;
    set<unsigned long> difference;

    set_intersection(
            t_itr, t_itr_end,
            r_itr, r_itr_end,
            inserter(intersection, intersection.begin()));
    
    set_symmetric_difference(
            t_itr, t_itr_end,
            r_itr, r_itr_end,
            inserter(difference, difference.begin()));
   
    double score = 0.0;
    
    set<unsigned long>::iterator v_itr;

        for( v_itr = intersection.begin(); v_itr != intersection.end() ; v_itr++ ){
            score += ref->getScore( *v_itr, true);
        }

        for( v_itr = difference.begin(); v_itr != difference.end() ; v_itr++ ){
            score += ref->getScore( *v_itr, false);
        }
    
    return score;
}

// This code was not generated by hand. Modifying it is ill-advised.
inline void runCombs( set<unsigned long>* hm, unsigned long xx, unsigned long mask, int ii ){

if(!((0x3B2D7 & mask)) && (ii >  18)) { hm->insert(((0x38000 & xx) >> 6) | ((0x3000 & xx) >> 5) | ((0x200 & xx) >> 3) | ((0xc0 & xx) >> 2) | ((0x10 & xx) >> 1) | ((0x7 & xx) >> 0)); }
if(!((0x3C44D7 & mask)) && (ii >  22)) { hm->insert(((0x3c0000 & xx) >> 10) | ((0x4000 & xx) >> 7) | ((0x400 & xx) >> 4) | ((0xc0 & xx) >> 2) | ((0x10 & xx) >> 1) | ((0x7 & xx) >> 0)); }
if(!((0xFFF & mask)) && (ii >  12)) { hm->insert(((0xfff & xx) >> 0)); }
if(!((0x1A18AF & mask)) && (ii >  21)) { hm->insert(((0x180000 & xx) >> 10) | ((0x20000 & xx) >> 9) | ((0x1800 & xx) >> 5) | ((0x80 & xx) >> 2) | ((0x20 & xx) >> 1) | ((0xf & xx) >> 0)); }
if(!((0x7747 & mask)) && (ii >  15)) { hm->insert(((0x7000 & xx) >> 5) | ((0x700 & xx) >> 4) | ((0x40 & xx) >> 3) | ((0x7 & xx) >> 0)); }

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

set<unsigned long>* extractFeatures(string seq, set<unsigned long> *res){
    
    int len = seq.length();

    if( !res ){
        res = new set<unsigned long>();
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
