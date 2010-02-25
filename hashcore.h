#include <utility>
#include <list>
#include <set>
#include <map>
#include <string>
#include <iterator>

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// TODO Use templates so that I can abstract out the feature sets and inherit from the supertype


class OddsTable {
    public:
        OddsTable(int);
        ~OddsTable();
        double oddsGivenFeature(int);
        double oddsWithoutFeature(int);
    private:
        int     _size;
        double *_contains;
        double *_doesnot;
};

class FeatureSet {
    public:
        std::string getName();
    protected:
        std::string _sequence;
        std::string _name;
};

class LikelihoodSet;
class ObservationSet;
class ReferenceSet;

class ReferenceSet: public FeatureSet {
    public:
        ReferenceSet(std::list<ObservationSet*>);
        bool addObservation(ObservationSet*);
        int size();
        double getNull();
        LikelihoodSet *makeLikelihood(ObservationSet*);
        friend class LikelihoodSet;
    protected:
        std::map<unsigned long, int> *_refs;
        std::map<unsigned long, double> *_withOdds;
        //std::map<unsigned long, double> *_withoutOdds;
 
    private:
        int _size;
        OddsTable *_odds;
        double _nullOdds;
};

class ObservationSet: public FeatureSet {
    public:
        ObservationSet(std::string sequence, std::string name = std::string(""));
        ~ObservationSet();
        std::set<unsigned long>* getFeatureSet();
    protected:
        std::map<unsigned long, int>* _obsset;
    private:
        friend class LikelihoodSet;
};

class LikelihoodSet: public FeatureSet {
    public:
        LikelihoodSet(ObservationSet* obs, ReferenceSet* corpus);
        ~LikelihoodSet();
        double likelihood(ObservationSet*);
    private:
        std::map<unsigned long, double>* _lset; 
        double _nullOdds;
};


enum nucleotide {A,T,C,G};
unsigned short getNucleotide(char);
void runCombs(std::map<unsigned long, int>*, unsigned long, unsigned long, int);
void insertBump(std::map<unsigned long, int>*, unsigned long);
std::map<unsigned long, int>* extractFeatures(std::string, std::map<unsigned long, int>* = NULL);

