#include <assert.h>

#include <queue>
#include <vector>
#include <string>

using namespace std;

class qread {
    public:
        std::string      _seq;
        std::vector<int> _qual;
};


int binomialTable[] = { 0, 1, 2, 3, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 10, 11, 11, 12, 12, 13, 13, 13, 14, 14, 15, 15, 15, 16, 16, 16, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20, 21, 21, 22, 22, 22, 23, 23, 23, 24, 24, 24, 25, 25, 25, 26, 26, 26, 27, 27, 27, 28, 28, 28, 29, 29, 29, 30, 30, 30, 31, 31, 31 };

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

string reverse_complement(string s){
    string res;
    res.resize(s.size());
    reverse(s.begin(), s.end());
    transform(s.begin(),s.end(), res.begin(), complement);
    return res;
}


qread align(qread readA, qread readB){

    assert(readA._qual.size() == readA._seq.size());
    assert(readB._qual.size() == readB._seq.size());

    string seqA = readA._seq;
    string seqB = reverse_complement(readB._seq);
    vector<int> qualA = readA._qual;
    vector<int> qualB = readB._qual;
    reverse(qualB.begin(), qualB.end());


    //printf("SeqA: %s (%d)\nSeqB: %s (%d)\n", seqA.c_str(), seqA.size(), seqB.c_str(), seqB.size());
    int minLen = min(seqA.size(), seqB.size());
    priority_queue<pair<double, int> > res;

     for(int ii = 5 ; ii <= minLen; ii++){
        int hits = 0;
        //printf("%d\t%d\n", seqA.size() - ii, ii);
        //printf("%s\n%s\n\n", seqA.substr(seqA.size() - ii, ii).c_str(), seqB.substr(0,ii).c_str());
        for(int jj = 0; jj < ii; jj++){
            char baseA = seqA[seqA.size() - 1 - jj];
            char baseB = seqB[ii - jj - 1];
            //printf("%c %c (%d %d)\n", baseA, baseB, seqA.size() - 1 - jj, ii - jj);
            if( goodBase(baseA) && goodBase(baseB) ){
                hits += (normalizeBase(baseA) == normalizeBase(baseB) ? 1 : 0);
            }
        }

        if( hits > binomialTable[ii] ){
            double score = (double) hits / (double) ii;
            //printf("Pushing <%f,%d>\n", score, ii);
            res.push(pair<double, int>(score, ii));
        }
    }

        qread result;
        int overlap;

        int res_length = 0;
        // If we have no evidence that indicates that our paired reads overlap:
        //  i.e. no one crossed the binomial threshhold.
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
            printf("%s\n", result._seq.c_str());
            string bPart = seqB.substr(overlap, seqB.size() - overlap);

            printf("%s\n", bPart.c_str());

            result._qual.insert(result._qual.end(), qualA.begin(), qualA.end() - overlap);

            int al = seqA.size();
            int bl = seqB.size();

            for(int ii = 0; ii < overlap; ii++ ){
                int aind = al - overlap + ii;
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
                result._seq += bPart;
                result._qual.insert(result._qual.end(), qualB.begin() + overlap, qualB.end());
            }
        }

        printf("Assembly: %s (%d)\n", result._seq.c_str(), (int) result._seq.size());
        printf("Qual length: %d\n", (int) result._qual.size());
        vector<int>::iterator q_itr;

        for( q_itr = result._qual.begin(); q_itr != result._qual.end(); q_itr++ ){
            printf("%d ", *q_itr);
        }
        printf("\n");
        //assert(result._qual.size() == result._seq.size());
        return result;
}

int main(int argc, char **argv){
    qread polyA, polyT;

    
    polyA._seq = string("AAACCCTTTT");
    polyA._qual = vector<int>(polyA._seq.size(),10);
    polyT._seq = string("GGGAAAAGG");
    polyT._qual = vector<int>(polyT._seq.size(),10);
    qread r = align(polyA, polyT);
    printf("%s\n", r._seq.c_str());
    return 0;
}
