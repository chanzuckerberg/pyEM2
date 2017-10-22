#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ExpressionMatrix.hpp"

namespace py = pybind11;
namespace czi_em2 = ChanZuckerberg::ExpressionMatrix2;

// The constructor signature for ExpressionMatrix is a directory name string
// and an ExpressionMatrixCreationParameters POD with four uint64_ts. In python
// we'd rather not require separate parameters object; instead we just want some
// optional keyword arguments. This function wraps the creation of the parameters
// object and the ExpressionMatrix itself in a way that can be exposed to
// py::init
std::unique_ptr<czi_em2::ExpressionMatrix> intialize_expression_matrix(
    std::string directory_name, uint64_t cell_capacity, uint64_t gene_capacity,
    uint64_t cell_metadata_name_capacity, uint64_t cell_metadata_value_capacity)
{
    auto creation_params = czi_em2::ExpressionMatrixCreationParameters();
    creation_params.geneCapacity = gene_capacity;
    creation_params.cellCapacity = cell_capacity;
    creation_params.cellMetaDataNameCapacity = cell_metadata_name_capacity;
    creation_params.cellMetaDataValueCapacity = cell_metadata_value_capacity;

    return std::unique_ptr<czi_em2::ExpressionMatrix>(
        new czi_em2::ExpressionMatrix(directory_name, creation_params));
}

PYBIND11_MODULE(pyEM2, m) {
    
    // There are two ways to construct an ExpressionMatrix in the C++ code,
    // 1. With a directory_name that doesn't exist and a creation parameters object
    // 2. With an existing directory_name that was created with constructor 1
    //
    // We're going to use method 1 in the __init__ of the python ExpressionMatrix
    // class, and we'll expose method 2 via a static method called
    // "from_existing_directory".
    py::class_<czi_em2::ExpressionMatrix>(m, "ExpressionMatrix")
        .def(py::init(&intialize_expression_matrix),
                py::arg("directory_name"),
                py::arg("cell_capacity")=1<<18,
                py::arg("gene_capacity")=1<<24,
                py::arg("cell_metadata_name_capacity")=1<<16,
                py::arg("cell_metadata_value_capacity")=1<<28)

        .def_static("from_existing_directory",
            [](std::string directory_name) {
                return std::unique_ptr<czi_em2::ExpressionMatrix>(
                    new czi_em2::ExpressionMatrix(directory_name));
            },
            py::arg("existing_em2_directory")
        )

        // Binding to EM.addGene, but takes a list of strings
        .def("add_genes", [](czi_em2::ExpressionMatrix& e, const std::vector<std::string>& gene_names) {
            for(auto gene_name: gene_names) {
                e.addGene(gene_name);
            };
        }, py::arg("gene_name_list"))



    ;
}
