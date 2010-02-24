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

    _obs.clear();
    _likelihoods.clear();

}

priority_queue<pair<double, string> > SequenceHasher::hash(string sequence){

    if( !_initialized ){
        initialize();
    }

    ObservationSet os(sequence);
    list<LikelihoodSet*>::iterator l_itr;

    priority_queue<pair<double, string> > rpq;

    for(l_itr = _likelihoods.begin(); l_itr != _likelihoods.end(); l_itr++){

        double likelihood = (*l_itr)->likelihood(&os);
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
    _likelihoods.clear();
    
    if( !_ref ){ delete _ref; }

    _ref = new ReferenceSet(_obs);
    list<ObservationSet*>::iterator obs_itr;
    for( obs_itr = _obs.begin(); obs_itr != _obs.end(); obs_itr++ ){
        LikelihoodSet *ls = new LikelihoodSet(*obs_itr, _ref);
        _likelihoods.push_back(ls);
    }
    return;
}

Object SequenceHasher::py_hash( const Tuple &args ){


    args.verify_length(1);
    String sequence = args[0];

    priority_queue<pair<double, string> > rpq = hash(sequence);

    Dict result;
    for(int ii = 0; ii < 5 && ii < (int) rpq.size(); ++ii){
        pair<double, string> entry = rpq.top();
        result[String(entry.second)] = Float(entry.first);
        rpq.pop();
    }

    return result;
}

Object SequenceHasher::py_addReference( const Tuple &args ){
    args.verify_length(2);
    String name = args[0];
    String sequence = args[1];
    
    addReference(name, sequence);
    return args[0];
}

Object SequenceHasher::py_initialize( const Tuple &args ){

    initialize();
    return Int((int) _obs.size());
}

Object SequenceHasher::getattr( const char *name ){
    return getattr_methods(name);    
}

Object SequenceHasher::repr(){
    return Py::String("Multiple Alignment Object");
}

void SequenceHasher::init_type(){
    behaviors().name("SequenceHasher");
    behaviors().doc("SequenceHasher objects: nil");
    behaviors().supportGetattr();
    behaviors().supportRepr();

    add_varargs_method("addReference",  &SequenceHasher::py_addReference, "addReference(name, sequence): add a reference entry");
    add_varargs_method("initialize",  &SequenceHasher::py_initialize, "initialize(): initialize the hasher");
    add_varargs_method("hash",  &SequenceHasher::py_hash, "hash(sequence): hash a sequence, returning the best matches in order");
    add_varargs_method("reference_count", &SequenceHasher::reference_count);
}

class hasher_module : public Py::ExtensionModule<hasher_module>
{
public:
    hasher_module()
    : Py::ExtensionModule<hasher_module>( "sequenceHasher" ) // this must be name of the file on disk e.g. simple.so or simple.pyd
    {
        SequenceHasher::init_type();
        add_varargs_method("SequenceHasher",&hasher_module::new_seqHasher,"SequenceHasher()");
        initialize( "documentation for the simple module" );
    }

    virtual ~hasher_module()
    {}

private:
    Object new_seqHasher(const Py::Tuple& args){
        return asObject(new SequenceHasher());
    }

};

extern "C" void initsequenceHasher()
{
#if defined(PY_WIN32_DELAYLOAD_PYTHON_DLL)
    Py::InitialisePythonIndirectPy::Interface();
#endif
    static hasher_module* seqHasher = new hasher_module();
}
