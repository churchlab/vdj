#include <CXX/Extensions.hxx>
#include <CXX/Objects.hxx>
#include <Python.h>

#include <string>

#include "malign.h"

class MAligner : public Py::PythonExtension<MAligner> {
	public:
        MAligner();
        ~MAligner();

        //python methods
        Py::Object addEntry(const Py::Tuple &args);
        Py::Object align(const Py::Tuple &args);

        //pycxx methods
        static void init_type();
        virtual Py::Object repr();
        virtual Py::Object getattr( const char *name);
        Py::Object reference_count(const Py::Tuple& args){
            return Py::Int(this->ob_refcnt);
        }
    private:
        MAlignerCore *_mcore;
}; 
