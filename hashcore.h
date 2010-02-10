#include <utility>
#include <list>
#include <set>
#include <map>
#include <string>
#include <iterator>

#include <assert.h>
#include <math.h>
#include <stdlib.h>

// #include <Python.h>

// TODO Use templates so that I can abstract out the feature sets and inherit from the supertype

enum nucleotide {A,T,C,G};

class OddsTable {
    public:
        OddsTable(int);
        ~OddsTable();

        double getFeaturelessOdds();
        double oddsGivenFeature(int);
        double oddsWithoutFeature(int);
    private:
        int     _size;
        double _featurelessOdds;
        double *_contains;
        double *_doesnot;
};

class FeatureSet {
    public:
        std::string getName();
        std::set<unsigned long> *getFeatureSet();
    protected:
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
        std::map<unsigned long, double> *_withoutOdds;
 
    private:
        int _size;
        OddsTable *_odds;
        double _nullOdds;
        bool addPseudocounts();
        bool _pseudocounts;
};

class ObservationSet: public FeatureSet {
    public:
        ObservationSet(char *name, char *sequence);
        ObservationSet(std::string name, std::string sequence);
        ~ObservationSet();
        std::set<unsigned long> *getFeatureSet();
    protected:
        std::map<unsigned long, int>* _obsset;
    private:
        char *_sequence;
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



unsigned short getNucleotide(char);
void runCombs(std::map<unsigned long, int>*, unsigned long, unsigned long, int);
void insertBump(std::map<unsigned long, int>*, unsigned long);
std::map<unsigned long, int>* extractFeatures(char*, std::map<unsigned long, int>*);
void pseudocount(std::map<unsigned long, int>*);

