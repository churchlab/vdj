#include "pairedAlign.h"

using namespace std;

bool goodBase(char n){
    switch(n){
        case 'a':
        case 'A': 
        case 'c':
        case 'C':
        case 't':
        case 'T':
        case 'g':
        case 'G': return true;
        default:  return false;
    }
}

char normalizeBase(char n){
    switch(n){
        case 'a':
        case 'A': return 'A';
        case 'c':
        case 'C': return 'C';
        case 't':
        case 'T': return 'T';
        case 'g':
        case 'G': return 'G';
        default:  return n;
    }
}

char complement(char n){
    switch(n){
        case 'a': return 't';
        case 'A': return 'T';
        case 'c': return 'g';
        case 'C': return 'G';
        case 't': return 'a';
        case 'T': return 'A';
        case 'g': return 'c';
        case 'G': return 'C';
        default:  return n;
    }
}

qread align(qread readA, qread readB){

    assert(readA._qual.size() == readA._seq.size());
    assert(readB._qual.size() == readB._seq.size());

    string seqA = readA._seq;
    string seqB = readB._seq;
    vector<int> qualA = readA._qual;
    vector<int> qualB = readB._qual;


    int minLen = min(seqA.size(), seqB.size());
    priority_queue<pair<double, int> > res;

     for(int ii = 10 ; ii < minLen; ii++){
        int hits = 0;
        for(int jj = 0; jj < ii; jj++){
            char baseA = seqA[seqA.size() - 1 - jj];
            char baseB = seqB[seqB.size() - 1 - jj];
            if( goodBase(baseA) && goodBase(baseB) ){
                hits += (normalizeBase(baseA) == complement(baseB) ? 1 : 0);
            }
        }

        if( hits > binomialTable[ii] ){
            double score = (double) hits / (double) ii;
            res.push(pair<double, int>(score, ii));
        }
    }

        qread result;
        int overlap;

        // Reverse complement one of the read pair ends.
        reverse(seqB.begin(), seqB.end());
        string::iterator pos;
        for( pos = seqB.begin(); pos != seqB.end() ; pos++ ){
            *pos = complement(*pos);
        }
        reverse(qualB.begin(), qualB.end()); 

        int res_length = 0;
        // We have no evidence that indicates that our paired reads overlap.
        if( res.empty() ){
            overlap = -1;

            res_length = seqA.size() + seqB.size() + 1;
            result._qual.reserve(res_length);
            result._qual = qualA;
            result._seq = seqA;
            // insert a Q00 gap between the sequences
            result._seq += '-';
            result._qual.push_back(0);

            result._seq += readB._seq; // tack on the other read
            result._qual.insert(result._qual.end(), qualB.begin(), qualB.end());
        } else {
            overlap = res.top().second;
            res_length = seqA.size() + seqB.size() - overlap;
            result._qual.reserve(res_length);
            result._seq = seqA.substr(0, seqA.size() - overlap);    
        
            string::iterator pos;
            for( pos = seqB.begin(); pos != seqB.end() ; pos++ ){
                *pos = complement(*pos);
            }

            int al = seqA.size();
            int bl = seqB.size();

            for(int ii = 0; ii < overlap; ii++ ){
                int aind = al - overlap + ii - 1;
                if( seqA[aind] == seqB[ii] ){
                    result._seq += seqB[ii];
                    result._qual.push_back(qualA[aind] + qualB[ii]);
                } else {
                    if( qualA[aind] >= qualB[ii] ){
                        result._seq += seqA[aind];
                    } else {
                        result._seq += seqB[ii];
                    }
                        result._qual.push_back(0);
                }
            }

            if( overlap < bl ){
                result._seq += seqB.substr(overlap, bl - overlap);
                result._qual.insert(result._qual.end(), result._qual.begin() + overlap, result._qual.end());
            }
        }

        assert(result._qual.size() == result._seq.size());
        return result;
}
