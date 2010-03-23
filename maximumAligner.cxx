#include "maximumAligner.hxx"

using namespace std;
using namespace Py;

MaximumAligner::MaximumAligner( const Tuple &args ){
    args.verify_length(0,1);

    if( args.length() == 1 ){
        vector<int> binomial;
        List pyBin = args[0];
        List::iterator itr;
        
        for(itr = pyBin.begin(); itr != pyBin.end(); itr++ ){
            binomial.push_back(Int(*itr));
        }
        _aligner = LikelihoodAligner(binomial); 
    } else {
        _aligner = LikelihoodAligner(); 
    }
}

MaximumAligner::~MaximumAligner(){

}

Object MaximumAligner::py_align( const Tuple &args ){
    args.verify_length(2,4);

    String seqA;
    String seqB;
    vector<int> qualA;
    vector<int> qualB;
    
    if( args.length() == 2 ){
        seqA = args[0];
        seqB = args[1];
        qualA = vector<int>(seqA.size(), 10);
        qualB = vector<int>(seqB.size(), 10);
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
    qread res = _aligner.align(readA, readB);

    List quals;
    Tuple r(2);
    r[0] = String(res._seq);
    
    vector<int>::iterator q_itr;
    for( q_itr = res._qual.begin(); q_itr != res._qual.end(); q_itr++ ){
        quals.append( Int(*q_itr) );
    }
    
    r[1] = quals;
    return r;
}

Object MaximumAligner::getattr( const char *name ){
    return getattr_methods(name);    
}

Object MaximumAligner::repr(){
    return Py::String("Paired End Aligner Object");
}

void MaximumAligner::init_type(){
    behaviors().name("MaximumAligner");
    behaviors().doc("MaximumAligner objects: nil");
    behaviors().supportGetattr();
    behaviors().supportRepr();

    add_varargs_method("align",     &MaximumAligner::py_align, "align(seq,[qual],seq,[qual]): align a sequence against the references");
    add_varargs_method("reference_count", &MaximumAligner::reference_count);
}

class maxaligner_module : public Py::ExtensionModule<maxaligner_module>
{
public:
    maxaligner_module()
    : Py::ExtensionModule<maxaligner_module>( "maximumAligner" ) // this must be name of the file on disk e.g. simple.so or simple.pyd
    {
        MaximumAligner::init_type();
        add_varargs_method("MaximumAligner",&maxaligner_module::new_maxaligner,"MaximumAligner()");
        initialize( "documentation for the simple module" );
    }

    virtual ~maxaligner_module()
    {}

private:
    Object new_maxaligner(const Py::Tuple& args){
        return asObject(new MaximumAligner(args));
    }

};

extern "C" void initmaximumAligner()
{
#if defined(PY_WIN32_DELAYLOAD_PYTHON_DLL)
    Py::InitialisePythonIndirectPy::Interface();
#endif
    static maxaligner_module* maxaligner = new maxaligner_module;
}
