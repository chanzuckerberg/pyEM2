#include <Eigen/Sparse>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "ClusterGraph.hpp"
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

void init_high_level(py::module& m) {
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
                    new czi_em2::ExpressionMatrix(directory_name, false));
            },
            py::arg("existing_em2_directory")
        )

        .def_static("from_csc_matrix",
            [](const Eigen::SparseMatrix<float>& csc_mtx, std::string directory_name) {
                auto em_ptr = intialize_expression_matrix(
                    directory_name,
                    csc_mtx.rows(),
                    csc_mtx.cols(),
                    csc_mtx.cols() * 100,
                    csc_mtx.cols() * 100);

                // Create some fake gene names for now
                for(int i=0; i<csc_mtx.rows(); ++i) {
                    std::string gene_name = "gene" + std::to_string(i);
                }

                std::vector<std::pair<std::string, float> > cell_expression_counts;
                for(int k=0; k < csc_mtx.outerSize(); ++k) {
                    cell_expression_counts.clear();
                    for(Eigen::SparseMatrix<float>::InnerIterator it(csc_mtx, k); it; ++it) {
                       std::string gene_name = "gene" + std::to_string(it.row());
                       cell_expression_counts.push_back(std::make_pair(gene_name, it.value()));
                    }
                    std::string cell_name = "cell" + std::to_string(k);
                    std::vector<std::pair<std::string, std::string> > cell_metadata;
                    cell_metadata.push_back(make_pair("CellName", cell_name));
                    em_ptr->addCell(cell_metadata, cell_expression_counts);
                }
            return em_ptr;
        })

        // Binding to EM.addGene, but takes a list of strings
        .def("add_genes", [](czi_em2::ExpressionMatrix& e, const std::vector<std::string>& gene_names) {
            for(auto gene_name: gene_names) {
                e.addGene(gene_name);
            };
        }, py::arg("gene_name_list"))

        .def_property_readonly("genes", [](czi_em2::ExpressionMatrix& e) {
            std::vector<czi_em2::GeneId> gene_ids(e.geneSets["AllGenes"].begin(),
                                                  e.geneSets["AllGenes"].end());
            std::vector<std::string> gene_names;
            for (auto gene_id : gene_ids) {
                gene_names.push_back(e.geneName(gene_id));
            }
            return gene_names;
        })

        // Binding to EM.addCell(metadata, counts). Args are two dicts, so on the
        // python side you don't have to worry about pairs.
        .def("add_cell",
            [](czi_em2::ExpressionMatrix& e,
               const std::map<std::string, std::string>& metadata,
               const std::map<std::string, float>& counts) {

                std::vector<std::pair<std::string, std::string> > cxx_metadata;
                std::vector<std::pair<std::string, float> > cxx_counts;
                cxx_metadata.assign(metadata.begin(), metadata.end());
                cxx_counts.assign(counts.begin(), counts.end());
                e.addCell(cxx_metadata, cxx_counts);

        }, py::arg("metadata"), py::arg("counts"))
    ;
}
