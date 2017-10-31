#include <pybind11/pybind11.h>
namespace py = pybind11;

void init_high_level(py::module&);
void init_low_level(py::module&);

PYBIND11_MODULE(pyEM2, m) {
    init_high_level(m);
    init_low_level(m);
}
