#include <list>

#include <string>
#include <queue>

#include "hashcore.h"

class SequenceHasher {
	public:
        SequenceHasher();
        ~SequenceHasher();

        //python methods
        void addReference(std::string, std::string);
        void initialize();
        void hash(std::string);
    private:
        bool _initialized;
        std::list<ObservationSet*> _obs;
        ReferenceSet *_ref;
        std::list<LikelihoodSet*> _likelihoods;
}; 
