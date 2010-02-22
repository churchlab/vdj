#include <CXX/Extensions.hxx>
#include <CXX/Objects.hxx>
#include <Python.h>

#include <string>

class Counter : public Py::PythonExtension<Counter> {
	public:
        Counter();
        ~Counter();
        int hello();

        //python methods
        Py::Object py_count(const Py::Tuple &args);

        //pycxx methods
        static void init_type();
        virtual Py::Object repr();
        virtual Py::Object getattr( const char *name);
    private:
        int _counter;
}; 
