#ifndef __MALIGNER_CORE__
#define __MALIGNER_CORE__

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <iostream>
#include <fstream>
#include <assert.h>

#include <algorithm>
#include <list>
#include <map>
#include <stack>
#include <string>
#include <queue>
#include <utility>
#include <vector>

#include "bandedalignment.h"
class MAlignerCore {
    public:
      
        MAlignerCore();
        ~MAlignerCore();
        void addEntry(std::string name, std::string sequence);
        std::pair<std::string, std::string> bestAlign(std::string input);
        std::string alignWith(std::string input, std::list<std::string> refs);
        std::list<std::pair<int, BandedAligner*> > align(std::string input);
        std::list<std::pair<int, BandedAligner*> > roundRobin(std::queue<BandedAligner*>);

    private:
        std::map<std::string,BandedAligner*> entries;
};

#endif
