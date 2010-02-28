#include <list>

#include <string>
#include <queue>
#include <utility>

#include "hashcore.h"

class SequenceHasher {
	public:
        SequenceHasher();
        ~SequenceHasher();

        //python methods
        void addReference(std::string, std::string);
        void initialize();
        std::priority_queue<std::pair<double, std::string> > hash(std::string);
    private:
        bool _initialized;
        std::list<ObservationSet*> _obs;
        ReferenceSet *_ref;
}; 
