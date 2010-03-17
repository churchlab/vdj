#include "pairedAligner.hxx"

using namespace std;
using namespace Py;

Object PairedAligner::py_align( const Tuple &args ){
    return args;
    /*
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
    qread res = align(readA, readB);

    List quals;
    Tuple r(2);
    r[0] = String(res._seq);
    
    vector<int>::iterator q_itr;
    for( q_itr = res._qual.begin(); q_itr != res._qual.end(); q_itr++ ){
        quals.append( Int(*q_itr) );
    }
    
    r[1] = quals;
    return r; */
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
