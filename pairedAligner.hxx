#include <CXX/Extensions.hxx>
#include <CXX/Objects.hxx>
#include <Python.h>

#include <assert.h>

#include <queue>
#include <string>
#include <vector>

typedef struct {
    std::string      _seq;
    std::vector<int> _qual;
} qread;

class PAlign {
    public:
        static qread align(qread, qread);
    private:
        static std::string reverse_complement(std::string s);
        static bool goodBase(char);
        static char normalizeBase(char);
        static char complement(char);
        static int binomialTable[];
};

int PAlign::binomialTable[] = { 0, 1, 2, 3, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 10, 11, 11, 12, 12, 13, 13, 13, 14, 14, 15, 15, 15, 16, 16, 16, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20, 21, 21, 22, 22, 22, 23, 23, 23, 24, 24, 24, 25, 25, 25, 26, 26, 26, 27, 27, 27, 28, 28, 28, 29, 29, 29, 30, 30, 30, 31, 31, 31 };


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
}; 
