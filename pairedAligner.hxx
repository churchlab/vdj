#include <CXX/Extensions.hxx>
#include <CXX/Objects.hxx>
#include <Python.h>

#include <assert.h>
#include <queue>
#include <string>
#include <vector>

//#include "pairedAligner.h"

class PairedAligner : public Py::PythonExtension<PairedAligner> {
	public:
        PairedAligner();
        ~PairedAligner();

        //python methods
        Py::Object py_align(const Py::Tuple &args);

        //pycxx methods
        static void init_type();
        virtual Py::Object repr();
        virtual Py::Object getattr( const char *name);
        Py::Object reference_count(const Py::Tuple& args){
            return Py::Int(this->ob_refcnt);
        }
    private:
};

