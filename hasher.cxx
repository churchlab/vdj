#include "hasher.hxx"

using namespace std;
using namespace Py;

SequenceHasher::SequenceHasher(){

    _ref = NULL;
    _initialized = false;

}

SequenceHasher::~SequenceHasher(){
    if( _ref ){
        delete _ref;
    }

    list<ObservationSet*>::iterator obs_itr;

    for( obs_itr = _obs.begin(); obs_itr != _obs.end(); obs_itr++ ){
        delete (*obs_itr);
    }
    _obs.clear();
}

priority_queue<pair<double, string> > SequenceHasher::hash(string sequence){

    if( !_initialized ){
        initialize();
    }

    ObservationSet os(sequence);
    list<ObservationSet*>::iterator l_itr;

    priority_queue<pair<double, string> > rpq;
    for(l_itr = _obs.begin(); l_itr != _obs.end(); l_itr++){

        double likelihood = (*l_itr)->likelihood(&os, _ref);
        string nm = (*l_itr)->getName();
        rpq.push(pair<double, string>(likelihood, nm));
    }
    
    return rpq;
}

void SequenceHasher::addReference(string name, string sequence){

    //TODO Ensure that arguments are in a consistent order
    ObservationSet* ob = new ObservationSet(sequence, name);
    _obs.push_back(ob);

    _initialized = false;
    return;
}

void SequenceHasher::initialize(){
    _initialized = true;
    if( !_ref ){ delete _ref; }

    _ref = new ReferenceSet(_obs);
    list<ObservationSet*>::iterator obs_itr;

    return;
}

Object SequenceHasher::py_addReference( const Tuple &args ){
    args.verify_length(2);
    String name = args[0];
    String sequence = args[1];

    addReference(name, sequence);

    return args[0];
}

Object SequenceHasher::py_hash( const Tuple &args ){
    args.verify_length(1);
    String sequence = args[0];
    priority_queue<pair<double, string> > results = hash(sequence);

    List d;
    
    while( !results.empty() ){
        string name = results.top().second;
        d.append(String(name));
        results.pop();
    }

    return d;
}

Object SequenceHasher::getattr( const char *name ){
    return getattr_methods(name);
}

Object SequenceHasher::repr(){
    return Py::String("A Hashed based Naive Bayes Classifier for sequence classification");
}

void SequenceHasher::init_type(){
    behaviors().name("SequenceHasher");
    behaviors().doc("SequenceHasher objects: nil");
    behaviors().supportGetattr();
    behaviors().supportRepr();

    add_varargs_method("addReference",  &SequenceHasher::py_addReference, "addReference(name, sequence)");
    add_varargs_method("hash",          &SequenceHasher::py_hash,          "hash(sequence)");
    add_varargs_method("reference_count", &SequenceHasher::reference_count);
}

class hasher_module : public Py::ExtensionModule<hasher_module>
{
public:
    hasher_module()
    : Py::ExtensionModule<hasher_module>("hasher")
    {
        SequenceHasher::init_type();
        add_varargs_method("SequenceHasher", &hasher_module::new_hasher);
        initialize( "docs" );
    }

    virtual ~hasher_module()
    {}

private:
    Object new_hasher(const Py::Tuple& args){
        return asObject(new SequenceHasher());
    }
};

extern "C" void inithasher()
{
#if defined(PY_WIN32_DELAYLOAD_PYTHON_DLL)
    Py::InitialisePythonIndirectPy::Interface();
#endif
    static hasher_module* hasher = new hasher_module;
}
