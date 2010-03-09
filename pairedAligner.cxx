#include "pairedAligner.hxx"

using namespace std;
using namespace Py;

bool PAlign::goodBase(char n){
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

char PAlign::normalizeBase(char n){
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

char PAlign::complement(char n){
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

qread PAlign::align(qread readA, qread readB){

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
       
            result._qual.insert(result._qual.end(), qualA.begin(), qualA.end() - overlap);

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

PairedAligner::PairedAligner(){
}

PairedAligner::~PairedAligner(){
}

Object PairedAligner::py_align( const Tuple &args ){
    args.verify_length(2,4);

    String seqA;
    String seqB;
    vector<int> qualA;
    vector<int> qualB;
    
    if( args.length() == 2 ){
        seqA = args[0];
        seqB = args[1];
        qualA = vector<int>(seqA.size(), 0);
        qualB = vector<int>(seqB.size(), 0);
    }

    if( args.length() == 3 ){
        args.verify_length(2);
    }

    if( args.length() == 4 ){
        seqA = String(args[0]);
        seqB = String(args[2]);
        List qA = List(args[1]);
        List qB = List(args[3]);
        
        qualA.reserve(qA.size());
        qualB.reserve(qB.size());

        List::iterator itr;

        for( itr = qA.begin() ; itr != qA.end() ; itr ++ ){
            qualA.push_back(Int(*itr));
        }

        for( itr = qB.begin() ; itr != qB.end() ; itr ++ ){
            qualB.push_back(Int(*itr));
        }
    }

    qread readA;
        readA._seq  = seqA;
        readA._qual = qualA;
    qread readB;
        readB._seq  = seqB;
        readB._qual = qualB;
    qread res = PAlign::align(readA, readB);

    List quals;
    Tuple r(2);
    r[0] = String(res._seq);
    
    vector<int>::iterator q_itr;
    for( q_itr = res._qual.begin(); q_itr != res._qual.end(); q_itr++ ){
        quals.append( Int(*q_itr) );
    }
    //FIXME I'm discarding qualities for now.
    //r[1] = quals;
    return r;

}

Object PairedAligner::getattr( const char *name ){
    return getattr_methods(name);    
}

Object PairedAligner::repr(){
    return Py::String("Paired End Aligner Object");
}

void PairedAligner::init_type(){
    behaviors().name("PairedAligner");
    behaviors().doc("PairedAligner objects: nil");
    behaviors().supportGetattr();
    behaviors().supportRepr();

    add_varargs_method("align",     &PairedAligner::py_align, "align(seq,[qual],seq,[qual]): align a sequence against the references");
    add_varargs_method("reference_count", &PairedAligner::reference_count);
}

class paligner_module : public Py::ExtensionModule<paligner_module>
{
public:
    paligner_module()
    : Py::ExtensionModule<paligner_module>( "pairedAligner" ) // this must be name of the file on disk e.g. simple.so or simple.pyd
    {
        PairedAligner::init_type();
        add_varargs_method("PairedAligner",&paligner_module::new_paligner,"PairedAligner()");
        initialize( "documentation for the simple module" );
    }

    virtual ~paligner_module()
    {}

private:
    Object new_paligner(const Py::Tuple& args){
        return asObject(new PairedAligner());
    }

};

extern "C" void initpairedAligner()
{
#if defined(PY_WIN32_DELAYLOAD_PYTHON_DLL)
    Py::InitialisePythonIndirectPy::Interface();
#endif
    static paligner_module* paligner = new paligner_module;
}
