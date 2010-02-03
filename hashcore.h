#include <utility>
#include <list>
#include <set>
#include <map>
#include <string>
#include <iterator>

#include <stdlib.h>

#include <Python.h>

enum nucleotide {A,T,C,G};
unsigned short getNucleotide(char);
void runCombs(std::map<unsigned long, int>*, unsigned long, unsigned long, int);
void insertBump(std::map<unsigned long, int>*, unsigned long);
std::map<unsigned long, int>* loadSequence(char*);
