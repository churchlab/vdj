#include <stdio.h>

#include "hello.hxx"

using namespace Py;

Counter::Counter(){
	_counter = 0;
}

Counter::~Counter(){}

Object Counter::getattr( const char *name ){
  /*
    if(std::string(name) == "count")
            return Int(_counter);
    if(std::string(name) == "name")
            return String("Counter Object");
    */
    return getattr_methods(name);    
}

Object Counter::repr(){
    return Py::String("Hello World!");
}

Object Counter::py_count( const Py::Tuple &args ){
    return Int(++_counter);
}

void Counter::init_type(){
    behaviors().name("Counter");
    behaviors().doc("Counter objects: nil");
    behaviors().supportGetattr();
    behaviors().supportRepr();

    add_varargs_method("countUp", &Counter::py_count, "print out a count");
}

class hello_module : public Py::ExtensionModule<hello_module>
{
public:
    hello_module()
    : Py::ExtensionModule<hello_module>( "hello" ) // this must be name of the file on disk e.g. simple.so or simple.pyd
    {
        Counter::init_type();

        add_varargs_method("hello", &hello_module::func, "documentation for func()");
        add_varargs_method("Counter",&hello_module::new_counter,"Counter()");


        // after initialize the moduleDictionary with exist
        initialize( "documentation for the simple module" );
    }

    virtual ~hello_module()
    {}

private:
    Object new_counter(const Py::Tuple& args){
        return asObject(new Counter() );
    }

    Object func( const Py::Tuple &args )
    {
        return String("HELLO WORLD!");
    }

};

extern "C" void inithello()
{
#if defined(PY_WIN32_DELAYLOAD_PYTHON_DLL)
    Py::InitialisePythonIndirectPy::Interface();
#endif

    static hello_module* simple = new hello_module;
}
