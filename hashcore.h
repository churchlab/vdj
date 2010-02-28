#include <utility>
#include <list>
#include <set>
#include <map>
#include <string>
#include <iterator>
#include <vector>

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
 };

class ObservationSet;
class ReferenceSet;

class ReferenceSet: public FeatureSet {
    public:
        ReferenceSet(std::list<ObservationSet*>);
        bool addObservation(ObservationSet*);
        int size();
        double getScore(unsigned long, bool);
    protected:
        std::map<unsigned long, int> *_refs;
        std::map<unsigned long, double> *_withOdds;
        std::map<unsigned long, double> *_withoutOdds;
 
    private:
        int _size;
        OddsTable *_odds;
};

class ObservationSet {
    public:
        ObservationSet(std::string sequence, std::string name = std::string(""));
        ~ObservationSet();
        std::string getName();
        double likelihood(ObservationSet*, ReferenceSet*);
        std::set<unsigned long>* getFeatureSet();
    protected:
        std::string _sequence;
        std::string _name;
        std::set<unsigned long>* _obsset;
};

enum nucleotide {A,T,C,G};
unsigned short getNucleotide(char);
void runCombs(std::map<unsigned long, int>*, unsigned long, unsigned long, int);
void insertBump(std::map<unsigned long, int>*, unsigned long);
std::set<unsigned long>* extractFeatures(std::string, std::set<unsigned long>* = NULL);

