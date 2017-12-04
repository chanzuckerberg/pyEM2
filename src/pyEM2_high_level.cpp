// Bindings for the "high-level" API that aims to transform the C++ interface into
// a pythonic one.

#include <Eigen/Sparse>            // For working with scipy.sparse matrices
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include "ClusterGraph.hpp"
#include "ExpressionMatrix.hpp"
#include "SimilarPairs.hpp"

namespace py = pybind11;
namespace czi_em2 = ChanZuckerberg::ExpressionMatrix2;


//////////////////////////////////////////////////////////////////////////////////////
// C++ Wrappers
//////////////////////////////////////////////////////////////////////////////////////

// The C++ version of EM2 keeps everything within a single ExpressionMatrix class. Ideally,
// in python, we'd like to be able to segregate the interface a little more and let users
// work with objects that have more tightly scoped responsibilities. We do that by creating
// little C++ wrapper classes that hold a shared_ptr to the big ExpressionMatrix object,
// some identifier like a name or id, and methods that delegate calls into the big EM object.


// The main wrapper class around ExpressionMatrix.
// This class has a couple responsiblities:
//   1. Enable the creation of shared_ptrs that other objects can hold.
//     We want to have potentially lots of objects that hold on to a pointer to the big EM
//     object, so it makes sense to use shared_ptrs for that. But, within a pybind11 method
//     we always get the object itself and never a pointer to it. So this class uses
//     enable_shared_from_this to let us create shared_ptrs that we can give to new objects
//     that we create.
//   2. Track counts of different objects created from the big EM object. The C++ interface
//     uses internal names and ids to refer to different entities like cells and graphs. In
//     python we're going to create new objects, so we don't want the user to have to worry
//     about the internal references. So we just keep counters, use them as names, and
//     increment them when we create something new.
class ExpressionMatrixWrapper : public std::enable_shared_from_this<ExpressionMatrixWrapper> {
public:

    // Keep a unique_ptr to the one and only EM object, as well as the directory name which
    // is needed to refer to memory mapped objects.
    std::unique_ptr<czi_em2::ExpressionMatrix> em_ptr;
    std::string directory_name;

    // Counters for different entity types.
    size_t cell_set_counter = 0;
    size_t gene_set_counter = 0;
    size_t similar_pairs_counter = 0;
    size_t cell_graph_counter = 0;
    size_t cluster_graph_counter = 0;

    // Returns a shared_ptr to this object.
    std::shared_ptr<ExpressionMatrixWrapper> getptr() {
        return shared_from_this();
    };
};

// The notion of "cell" in EM2. Cells have metadata and cells have expression counts.
// The expression counts are read only, but you can change the metadata.
class CellWrapper {
public:

    // Cells are fully defined by their CellId, so that's all we need to store.
    std::shared_ptr<ExpressionMatrixWrapper> wrapper_ptr;
    czi_em2::CellId cell_id;

    CellWrapper(std::shared_ptr<ExpressionMatrixWrapper> wrapper, czi_em2::CellId id) :
        wrapper_ptr(wrapper), cell_id(id) {};

    // Update the cell's metadata. This is unfortunately a little subtle:
    // cell.metadata['key'] = 'value' does nothing. But
    // cell.metadata = {'key': 'value'} updates the metadata.
    // There's got to be a way...
    void _set_metadata(py::dict& metadata) {
        for(auto& kv : metadata) {
            wrapper_ptr->em_ptr->setCellMetaData(cell_id, kv.first.cast<std::string>(), kv.second.cast<std::string>());
        }
    };

    // Get the cell's metadata. Note that this is a copy, so modifying doesn't affect
    // the underlying CellId in the EM.
    py::dict _get_metadata() const {
        py::dict metadata_dict;
        for(auto it : wrapper_ptr->em_ptr->getCellMetaData(cell_id)) {
            metadata_dict[py::cast(it.first)] = py::cast(it.second);
        }
        return metadata_dict;
    };

    // Handle CellName separately since it's required and has to be unique.
    void _set_name(std::string name) {
        wrapper_ptr->em_ptr->setCellMetaData(cell_id, "CellName", name);
    }

    std::string _get_name() const {
        return wrapper_ptr->em_ptr->getCellMetaData(cell_id, "CellName");
    }

    // Get a dict of {"gene_name": count, ...}
    std::map<std::string, float> _get_expression_counts() const {
        std::vector<czi_em2::CellId> id_vector;
        id_vector.push_back(cell_id);
        auto exp_counts_pairs = wrapper_ptr->em_ptr->getCellsExpressionCounts(id_vector);
        std::map<std::string, float> exp_counts;

        for(auto& pair : exp_counts_pairs[0]) {
            exp_counts[wrapper_ptr->em_ptr->geneName(pair.first)] = pair.second;
        }
        return exp_counts;
    }
};

// An EM2 "gene". Genes are pretty simple, they just have a name.
class GeneWrapper {
public:
    GeneWrapper(std::shared_ptr<ExpressionMatrixWrapper> wrapper, czi_em2::GeneId id) :
        wrapper_ptr(wrapper), gene_id(id) {};

    std::shared_ptr<ExpressionMatrixWrapper> wrapper_ptr;
    czi_em2::GeneId gene_id;

    std::string name() const {
        return wrapper_ptr->em_ptr->geneName(gene_id);
    }
};

// Expose the concept of similar pairs as an independent object. All a user really needs
// is the ability to create this object and the ability to query similar cells.
class SimilarPairsWrapper {
public:
    std::shared_ptr<ExpressionMatrixWrapper> wrapper_ptr;
    std::string full_name; // The entire name, including the "/SimilarPairs-"
    std::string name; // The name in the EM2 interface, so missing the directoryName/SimilarPairs...
    std::string cell_set_name;
    std::string gene_set_name;

    // Get a list of cell indexes that are similar to the cell at cell_index.
    std::vector<size_t> similar_cells(size_t cell_index) {
        czi_em2::SimilarPairs sps(full_name, false);
        auto cell_id = wrapper_ptr->em_ptr->getCellSet("AllCells")[cell_index];
        auto local_id = sps.getLocalCellId(cell_id);
        auto similar_cells = sps[local_id];
        std::vector<size_t> similar_indices;
        for(auto cell_pair : similar_cells) {
            auto global_id = sps.getGlobalCellId(cell_pair.first);
            similar_indices.push_back(wrapper_ptr->em_ptr->getCellSet("AllCells")[global_id]);
        }
        return similar_indices;
    }

};

// Wrap a CellGraphVertexInfo. This is a component of a CellGraph. The main thing here is to expose
// the CellWrapper rather than just the CellId. We're not telling people about CellIds.
class CellGraphVertexWrapper {
public:
    std::shared_ptr<ExpressionMatrixWrapper> wrapper_ptr;
    czi_em2::CellGraphVertexInfo cell_vertex;

    CellGraphVertexWrapper(std::shared_ptr<ExpressionMatrixWrapper> wrapper, czi_em2::CellGraphVertexInfo vertex) :
        wrapper_ptr(wrapper), cell_vertex(vertex) {};

    CellWrapper cell() {
        return CellWrapper(wrapper_ptr, cell_vertex.cellId);
    }

    float x() {return cell_vertex.x();}
    float y() {return cell_vertex.y();}
};

// Wrap a CellGraph.
class CellGraphWrapper {
public:
    std::shared_ptr<ExpressionMatrixWrapper> wrapper_ptr;
    std::string name;

    // Compute the layout on initialization. This is maybe not a great choice for huge datasets.
    CellGraphWrapper(std::shared_ptr<ExpressionMatrixWrapper> wrapper, std::string graph_name) :
            wrapper_ptr(wrapper), name(graph_name) {
        wrapper_ptr->em_ptr->computeCellGraphLayout(name);
    }

    // Return a list of CellGraphVertices
    std::vector<CellGraphVertexWrapper> vertices() {
        auto vinfos = wrapper_ptr->em_ptr->getCellGraphVertices(name);
        std::vector<CellGraphVertexWrapper> vwrappers;
        for(auto& vinfo : vinfos) {
            vwrappers.push_back(CellGraphVertexWrapper(wrapper_ptr, vinfo));
        }
        return vwrappers;
    }

    // Return a list of edges, which are pairs of CellWrappers
    std::vector<std::pair<CellWrapper, CellWrapper> > edges() {
        auto raw_edges = wrapper_ptr->em_ptr->getCellGraphEdges(name);
        std::vector<std::pair<CellWrapper, CellWrapper> > wrapped_edges;
        for(auto& raw_edge : raw_edges) {
            wrapped_edges.push_back(std::make_pair(CellWrapper(wrapper_ptr, raw_edge.first),
                                                   CellWrapper(wrapper_ptr, raw_edge.second)));

        }
        return wrapped_edges;
    }
};



//////////////////////////////////////////////////////////////////////////////////////
// Pybind code
//////////////////////////////////////////////////////////////////////////////////////


// The constructor signature for ExpressionMatrix is a directory name string
// and an ExpressionMatrixCreationParameters POD with four uint64_ts. In python
// we'd rather not require separate parameters object; instead we just want some
// optional keyword arguments. This function wraps the creation of the parameters
// object and the ExpressionMatrix itself in a way that can be exposed to
// py::init
std::shared_ptr<ExpressionMatrixWrapper> intialize_expression_matrix(
    std::string directory_name, uint64_t cell_capacity, uint64_t gene_capacity,
    uint64_t cell_metadata_name_capacity, uint64_t cell_metadata_value_capacity)
{
    czi_em2::ExpressionMatrixCreationParameters creation_params(
        gene_capacity,
        cell_capacity,
        cell_metadata_name_capacity,
        cell_metadata_value_capacity);

    std::unique_ptr<czi_em2::ExpressionMatrix> em_ptr(
        new czi_em2::ExpressionMatrix(directory_name, creation_params));

    auto wrapper_ptr = std::make_shared<ExpressionMatrixWrapper>(ExpressionMatrixWrapper());
    wrapper_ptr->em_ptr = std::move(em_ptr);
    wrapper_ptr->directory_name = directory_name;
    return wrapper_ptr;
}

void init_high_level(py::module& m) {
  // There are two ways to construct an ExpressionMatrix in the C++ code,
  // 1. With a directory_name that doesn't exist and a creation parameters object
  // 2. With an existing directory_name that was created with constructor 1
  //
  // We're going to use method 1 in the __init__ of the python ExpressionMatrix
  // class, and we'll expose method 2 via a static method called
  // "from_existing_directory".
    py::class_<ExpressionMatrixWrapper, std::shared_ptr<ExpressionMatrixWrapper> >(
      m,
      "ExpressionMatrix",
      "Matrix of gene expression values for different cells.")
      .def(py::init(
        [](std::string directory_name, uint64_t cell_capacity, uint64_t gene_capacity,
           uint64_t cell_metadata_name_capacity, uint64_t cell_metadata_value_capacity) {

           auto wrapper_ptr = intialize_expression_matrix(
              directory_name, cell_capacity, gene_capacity, cell_metadata_name_capacity,
              cell_metadata_value_capacity);
          return wrapper_ptr;
          }),
        "Create a new, empty ExpressionMatrix.",
        py::arg("directory_name"),
        py::arg("cell_capacity")=1<<18,
        py::arg("gene_capacity")=1<<24,
        py::arg("cell_metadata_name_capacity")=1<<16,
        py::arg("cell_metadata_value_capacity")=1<<28)

        .def_static("from_existing_directory",
            [](std::string directory_name) {
                std::unique_ptr<czi_em2::ExpressionMatrix> em_ptr(new czi_em2::ExpressionMatrix(directory_name, false));
                auto wrapper_ptr = std::make_shared<ExpressionMatrixWrapper>(ExpressionMatrixWrapper());
                wrapper_ptr->em_ptr = std::move(em_ptr);
                wrapper_ptr->directory_name = directory_name;
                return wrapper_ptr;
            },
            "Create an ExpressionMatrix from an previously-created directory.",
            py::arg("existing_em2_directory")
        )

        .def_static("from_csc_matrix",
            [](const Eigen::SparseMatrix<float, Eigen::ColMajor>& csc_mtx, std::string directory_name) {
                // Note that rows are genes, columns are cells, and the outer dimension must be columns.
                // Also note that this necessarily makes a copy of the matrix, vvv sad about that but there
                // doesn't seem to be a way around that with a sparse matrix.
                auto wrapper_ptr = intialize_expression_matrix(
                    directory_name,
                    csc_mtx.cols() * 2, // cell_capacity
                    csc_mtx.rows() * 2, // gene_capacity
                    csc_mtx.cols() * 30,
                    csc_mtx.cols() * 30);

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

                    wrapper_ptr->em_ptr->addCell(cell_metadata, cell_expression_counts);
                }
            return wrapper_ptr;
           },
           "Create an ExpressionMatrix from a scipy.sparse.csc_matrix. Note that dummy "
           "cell and gene names are generated automatically.",
           py::arg("column_major_sparse_matrix"),
           py::arg("directory_name")
        )

        // Binding to EM.addGene, but takes a list of strings
        .def("add_genes", [](ExpressionMatrixWrapper& self, const std::vector<std::string>& gene_names) {
            for(auto gene_name: gene_names) {
                self.em_ptr->addGene(gene_name);
            };
        },
        "Add a list of genes to the ExpressionMatrix.",
        py::arg("gene_name_list"))

        .def_property_readonly("genes",
            [](ExpressionMatrixWrapper& self) {
                py::list genes;
                for (auto gene_id : self.em_ptr->geneSets["AllGenes"]) {
                    GeneWrapper gene(self.getptr(), gene_id);
                    genes.append(py::cast(gene));
                }
                return genes;
        })

        .def_property_readonly("cells",
            [](ExpressionMatrixWrapper& self) {
                auto cell_ids(self.em_ptr->getCellSet("AllCells"));
                py::list cells;
                for(auto id : cell_ids) {
                    CellWrapper cell(self.getptr(), id);
                    cells.append(py::cast(cell));
                }
                return cells;
            })

        // Binding to EM.addCell(metadata, counts). Args are two dicts, so on the
        // python side you don't have to worry about pairs.
        .def("add_cell",
            [](ExpressionMatrixWrapper& self,
               const std::string name,
               const std::map<std::string, std::string>& metadata,
               const std::map<std::string, float>& counts) {

                std::vector<std::pair<std::string, std::string> > cxx_metadata;
                std::vector<std::pair<std::string, float> > cxx_counts;
                cxx_metadata.assign(metadata.begin(), metadata.end());
                cxx_metadata.push_back(std::make_pair("CellName", name));
                cxx_counts.assign(counts.begin(), counts.end());
                self.em_ptr->addCell(cxx_metadata, cxx_counts);

        },
        "Add a cell to the ExpressionMatrix.",
        py::arg("name"),
        py::arg("metadata"),
        py::arg("counts"))

    .def("create_similar_pairs",
        [](ExpressionMatrixWrapper& self,
           double similarity_threshold,
           size_t max_pairs_per_cell,
           size_t lsh_count,
           unsigned int seed,
           py::object& cell_set,
           py::object& gene_set
           ) {

        std::string cell_set_name;
        std::string gene_set_name;

        if(!cell_set.is_none()) {
           std::vector<czi_em2::CellId> cell_ids;
           for(auto cell : cell_set) {
               cell_ids.push_back(cell.cast<CellWrapper>().cell_id);
           }
           cell_set_name = std::to_string(self.cell_set_counter);
           self.cell_set_counter++;
           self.em_ptr->createCellSet(cell_set_name, cell_ids);
        } else {
            cell_set_name = "AllCells";
        }

        if(!gene_set.is_none()) {
           std::vector<std::string> gene_names;
           for(auto gene : gene_set) {
               gene_names.push_back(gene.cast<GeneWrapper>().name());
           }
           gene_set_name = std::to_string(self.gene_set_counter);
           self.gene_set_counter++;
           int ignored = 0;
           int empty = 0;
           self.em_ptr->createGeneSetFromGeneNames(gene_set_name, gene_names, ignored, empty);
        } else {
            gene_set_name = "AllGenes";
        }

        std::string name = std::to_string(self.similar_pairs_counter);
        self.similar_pairs_counter++;

        self.em_ptr->findSimilarPairs4(gene_set_name, cell_set_name, name, max_pairs_per_cell,
                                       similarity_threshold, lsh_count, seed);

        auto sp_ptr = std::unique_ptr<SimilarPairsWrapper>(new SimilarPairsWrapper());
        sp_ptr->wrapper_ptr = self.getptr();

        // findSimilarPairs4 prepends this to the name, and we'll need this to look it up
        // later.
        std::string full_name = self.directory_name + "/SimilarPairs-" + name;
        sp_ptr->full_name = full_name;
        sp_ptr->name = name;
        sp_ptr->cell_set_name = cell_set_name;
        sp_ptr->gene_set_name = gene_set_name;
        return sp_ptr;
        },
        "Create a new SimilarPairs object. The SimilarPairs object can then be used to "
        "identify approximate nearest neighbors for cells in the ExpressionMatrix.",
        py::arg("similarity_threshold")=0,
        py::arg("max_pairs_per_cell")=10,
        py::arg("lsh_count")=1024,
        py::arg("seed")=42,
        py::arg("cell_set")=py::none(),
        py::arg("gene_set")=py::none())

    .def("create_cell_graph",
        [](ExpressionMatrixWrapper& self,
           double similarity_threshold,
           size_t max_pairs_per_cell,
           const SimilarPairsWrapper& similar_pairs,
           py::object cell_set) {

        std::string cell_set_name;
        if(!cell_set.is_none()) {
           std::vector<czi_em2::CellId> cell_ids;
           for(auto& cell : cell_set) {
               cell_ids.push_back(cell.cast<CellWrapper>().cell_id);
           }
           cell_set_name = std::to_string(self.cell_set_counter);
           self.cell_set_counter++;
           self.em_ptr->createCellSet(cell_set_name, cell_ids);
        } else {
            cell_set_name = similar_pairs.cell_set_name;
        }

        std::string cell_graph_name = std::to_string(self.cell_graph_counter);
        self.cell_graph_counter++;

        self.em_ptr->createCellGraph(
            cell_graph_name, cell_set_name,
            similar_pairs.name, similarity_threshold,
            max_pairs_per_cell);

        std::unique_ptr<CellGraphWrapper> cg_ptr(new CellGraphWrapper(self.getptr(), cell_graph_name));

        return cg_ptr;
        },
        "Create a CellGraph.",
        py::arg("similarity_threshold")=0,
        py::arg("max_pairs_per_cell")=10,
        py::arg("similar_pairs"),
        py::arg("cell_set")=py::none())
    ;

    py::class_<SimilarPairsWrapper>(
      m,
      "SimilarPairs",
      "SimilarPairs objects allow approximate nearest neighbory queries for cells "
      "in an ExpressionMatrix.")
      .def("get_similar_cells",
         &SimilarPairsWrapper::similar_cells,
         "Get the indices of cells similar to the cells at cell_index.",
         py::arg("cell_index"))
    ;

    py::class_<CellWrapper>(
      m,
      "Cell",
      "A single cell from the ExpressionMatrix. A Cell has a name, metadata, and expression counts")
      .def(py::init<std::shared_ptr<ExpressionMatrixWrapper>, czi_em2::CellId>())
      .def_property("metadata", &CellWrapper::_get_metadata, &CellWrapper::_set_metadata)
      .def_property_readonly("expression_counts", &CellWrapper::_get_expression_counts)
      .def_property("name", &CellWrapper::_get_name, &CellWrapper::_set_name);
    ;

    py::class_<GeneWrapper>(
      m,
      "Gene",
      "A single gene from an ExpressionMatrix. Genes are simple, and just have a name.")
      .def(py::init<std::shared_ptr<ExpressionMatrixWrapper>, czi_em2::GeneId>())
      .def_property_readonly("name", &GeneWrapper::name)
      .def("__str__", &GeneWrapper::name)
    ;

    py::class_<CellGraphVertexWrapper>(
      m,
      "CellGraphVertex",
      "A vertex in a CellGraph.")
      .def(py::init<std::shared_ptr<ExpressionMatrixWrapper>, czi_em2::CellGraphVertexInfo>())
      .def_property_readonly("cell", &CellGraphVertexWrapper::cell)
      .def_property_readonly("x", &CellGraphVertexWrapper::x)
      .def_property_readonly("y", &CellGraphVertexWrapper::y)
    ;

    py::class_<CellGraphWrapper>(
      m,
      "CellGraph",
      "A graph of cells in an ExpressionMatrix.")
      .def(py::init<std::shared_ptr<ExpressionMatrixWrapper>, std::string>())
      .def_property_readonly("vertices", &CellGraphWrapper::vertices)
      .def_property_readonly("edges", &CellGraphWrapper::edges)
    ;
}
