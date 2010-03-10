#include "maligner.hxx"

using namespace std;
using namespace Py;

MAligner::MAligner(){
    _mcore = new MAlignerCore();
}

MAligner::~MAligner(){
    delete _mcore;
}

Object MAligner::addEntry( const Tuple &args ){
    args.verify_length(2);
    String name = args[0];
    String sequence = args[1];
    
    _mcore->addEntry(name, sequence);
    return args[0];
}

Object MAligner::align( const Tuple &args ){
    args.verify_length(1, 2);

    String sequence = args[0];
    pair<string, pair<string, string> > res;
    string name;
    pair<string, string> trace;
    if( args.length() == 1 ){
        res = _mcore->align(sequence);
    } else {
        List py_refs = args[1];
        SeqBase<Object>::iterator itr;

        list<string> refs;
        for( itr = py_refs.begin() ; itr != py_refs.end(); itr++ ){
            string s = String(*itr);
            refs.push_back(s);    
        }

        res = _mcore->align(sequence,refs);
    }
        
    name = res.first;
    trace = res.second;
    Tuple t(3);
    t[0] = String(name);
    t[1] = String(trace.first);
    t[2] = String(trace.second);
    return t;
}

/*
Object MAligner::alignWith( const Tuple &args ){
    args.verify_length(2);
    String sequence = args[0];
    List   refs = args[1];

    // dummy return
    return String("");
}
*/

Object MAligner::getattr( const char *name ){
    return getattr_methods(name);    
}

Object MAligner::repr(){
    return Py::String("Multiple Alignment Object");
}

void MAligner::init_type(){
    behaviors().name("MAligner");
    behaviors().doc("MAligner objects: nil");
    behaviors().supportGetattr();
    behaviors().supportRepr();

    add_varargs_method("addReference",  &MAligner::addEntry, "addReference(name, sequence): add an aligner reference");
    add_varargs_method("align",     &MAligner::align, "align(sequence,[refs]): align a sequence against the references");
    add_varargs_method("reference_count", &MAligner::reference_count);
}

class maligner_module : public Py::ExtensionModule<maligner_module>
{
public:
    maligner_module()
    : Py::ExtensionModule<maligner_module>( "maligner" ) // this must be name of the file on disk e.g. simple.so or simple.pyd
    {
        MAligner::init_type();
        add_varargs_method("MAligner",&maligner_module::new_maligner,"MAligner()");
        initialize( "documentation for the simple module" );
    }

    virtual ~maligner_module()
    {}

private:
    Object new_maligner(const Py::Tuple& args){
        return asObject(new MAligner());
    }

};

extern "C" void initmaligner()
{
#if defined(PY_WIN32_DELAYLOAD_PYTHON_DLL)
    Py::InitialisePythonIndirectPy::Interface();
#endif
    static maligner_module* maligner = new maligner_module;
}
