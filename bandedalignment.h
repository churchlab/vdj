#include "bandedmatrix.h"
#include <string>

class BandedAlignment {
    public:
        BandedAlignment(std::string);
        void initialize(std::string);
        int align();
    private:
        std::string _ref;
        std::string _test;
        bool _initialized;
        BandedMatrix<int> _matrix;
};
