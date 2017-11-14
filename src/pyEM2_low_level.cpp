// The "low-level" interface that directly binds C++ methods and classes

#include "ClusterGraph.hpp"
#include "ExpressionMatrix.hpp"
#include "MemoryMappedVector.hpp"
#include "MemoryMappedVectorOfLists.hpp"
using namespace ChanZuckerberg;
using namespace ExpressionMatrix2;

// Pybind11.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
using namespace pybind11;

class LowLevelExpressionMatrix : public ExpressionMatrix {
    // Use the same constructors as ExpressionMatrix
    using ExpressionMatrix::ExpressionMatrix;
};

void init_low_level(module& m) {
    
    module low_level_module = m.def_submodule("low_level", "Low-level access to C++ implementation of ExpressionMatrix2");

    low_level_module.doc() =
        "Software for analysis, visualization, and clustering "
        "of gene expression data from single-cell RNA sequencing.";

    // Class ExpressionMatrix.
    class_<LowLevelExpressionMatrix>(
        low_level_module,
        "ExpressionMatrix",
        "Top level object used to store and manipulate an expression matrix "
        "containing information on its genes and cells "
        "and many related data structures.")
       .def(init<string, ExpressionMatrixCreationParameters>(),
           "Construct a new, empty ExpressionMatrix.",
           arg("directoryName"),
           arg("parameters")
       )
       .def(init<string, bool>(),
           "Access an existing ExpressionMatrix.",
           arg("directoryName"),
           arg("allowReadOnly")=false
       )

       // Get the total number of genes or cells currently in the system.
       .def("geneCount",
           [](LowLevelExpressionMatrix& self) {return self.geneCount();},
           "Return the total number of genes."
       )
       .def("cellCount",
           [](LowLevelExpressionMatrix& self) {return self.cellCount();},
           "Return the total number of cells."
       )

       // Genes.
       .def("addGene",
           [](LowLevelExpressionMatrix& self, string geneName) {return self.addGene(geneName);},
           "Add a new gene with a given name.",
           arg("geneName")
           )
       .def("geneName",
           [](LowLevelExpressionMatrix& self, GeneId geneId) {return self.geneName(geneId);},
           "Return the name of a gene given its numeric id.",
           arg("geneId")
       )
       .def("geneIdFromName",
           [](LowLevelExpressionMatrix& self, string geneName) {return self.geneIdFromName(geneName);},
           "Return the numeric id of a gene given its name.",
           arg("geneName")
       )


       // Various ways to add cells.
       .def
       (
           "addCell",
          [](LowLevelExpressionMatrix& self, std::vector<std::pair<std::string, std::string> >& metadata,
             std::vector<std::pair<std::string, float> >& expressionCounts) {
                return self.addCell(metadata, expressionCounts);},
           "Adds a cell to the system. The cell expression counts "
           "and meta data are given in Python lists. "
           "metadata is a list of tuples with meta data "
           "(name, value) pairs. "
           "expressionCounts is a list of tuples (geneName, count). "
           "See `here <../../../PythonApi.html#addCell>`__ "
           "for an example. "
           "Returns the cell id of the cell that was just added. "
           "Cell ids begin at zero and increment by one each time a cell is added. ",
           arg("metaData"),
           arg("expressionCounts")
       )
       .def("addCells",
           [](LowLevelExpressionMatrix& self,
              string expressionCountsFileName,
              string expressionCountsFileSeparators,
              string cellMetaDataFileName,
              string cellMetaDataFileSeparators) {
             return self.addCells(expressionCountsFileName, expressionCountsFileSeparators,
                                  cellMetaDataFileName, cellMetaDataFileSeparators);},
           "Add cells and their meta data from input files in delimited format.",
           arg("expressionCountsFileName"),
           arg("expressionCountsFileSeparators") = ",",
           arg("cellMetaDataFileName"),
           arg("cellMetaDataFileSeparators") = ","
       )
       .def("addCellsFromHdf5",
           [](LowLevelExpressionMatrix& self,
              string fileName,
              string cellNamePrefix,
              vector<pair<string, string> > cellMetaDataArgument,
              double totalExpressionCountThreshold) {
             return self.addCellsFromHdf5(fileName, cellNamePrefix, cellMetaDataArgument,
                                          totalExpressionCountThreshold);},
           "Add cells from an input file in hdf5 format.",
           arg("fileName"),
           arg("cellNamePrefix"),
           arg("cellMetaDataArgument"),
           arg("totalExpressionCountThreshold")
       )
       .def("addCellsFromBioHub1",
           [](LowLevelExpressionMatrix& self,
              string expressionCountsFileName,
              size_t initialMetaDataCount,
              size_t finalMetaDataCount,
              string plateMetaDataFileName) {
             return self.addCellsFromBioHub1(expressionCountsFileName, initialMetaDataCount,
                                             finalMetaDataCount, plateMetaDataFileName);},
           "Add cells in the format used by the BioHub pipeline, July 2017.",
           arg("expressionCountsFileName"),
           arg("initialMetaDataCount"),
           arg("finalMetaDataCount"),
           arg("plateMetaDataFileName")
       )
       .def("addCellsFromBioHub2",
           [](LowLevelExpressionMatrix& self,
              string plateFileName,
              double totalExpressionCountThreshold) {
             return self.addCellsFromBioHub2(plateFileName, totalExpressionCountThreshold);},
           "Add cells in the format used by the BioHub pipeline, September 2017.",
           arg("plateFileName"),
           arg("totalExpressionCountThreshold")
       )
       .def("addCellMetaData",
           [](LowLevelExpressionMatrix& self, string cellMetaDataFileName) {
            return self.addCellMetaData(cellMetaDataFileName);},
           "Add meta data for existing cells from an input csv file.",
           arg("cellMetaDataFileName")
       )

       ;



    // Class ExpressionMatrixCreationParameters.
    class_<ExpressionMatrixCreationParameters>(
        low_level_module,
        "ExpressionMatrixCreationParameters",
        "Class used to store creation parameters for a new expression matrix.")
        .def(init<uint64_t, uint64_t, uint64_t, uint64_t>())
        .def_readwrite("geneCapacity", &ExpressionMatrixCreationParameters::geneCapacity)
        .def_readwrite("cellCapacity", &ExpressionMatrixCreationParameters::cellCapacity)
        .def_readwrite("cellMetaDataNameCapacity", &ExpressionMatrixCreationParameters::cellMetaDataNameCapacity)
        .def_readwrite("cellMetaDataValueCapacity", &ExpressionMatrixCreationParameters::cellMetaDataValueCapacity)
        ;



    // Class ClusterGraphCreationParameters.
    class_<ClusterGraphCreationParameters>(low_level_module, "ClusterGraphCreationParameters")
        .def(init<>())
        .def_readwrite("stableIterationCount", &ClusterGraphCreationParameters::stableIterationCount)
        .def_readwrite("maxIterationCount", &ClusterGraphCreationParameters::maxIterationCount)
        .def_readwrite("seed", &ClusterGraphCreationParameters::seed)
        .def_readwrite("minClusterSize", &ClusterGraphCreationParameters::minClusterSize)
        .def_readwrite("maxConnectivity", &ClusterGraphCreationParameters::maxConnectivity)
        .def_readwrite("similarityThreshold", &ClusterGraphCreationParameters::similarityThreshold)
        ;

    // Enum class NormalizationMethod.
    enum_<NormalizationMethod>(
        low_level_module,
        "NormalizationMethod",
        "Various ways to normalize gene expressions.")
        .value(normalizationMethodToShortString(NormalizationMethod::None).c_str(),       NormalizationMethod::None)
        .value(normalizationMethodToShortString(NormalizationMethod::L1).c_str(),         NormalizationMethod::L1)
        .value(normalizationMethodToShortString(NormalizationMethod::L2).c_str(),         NormalizationMethod::L2)
        .value(normalizationMethodToShortString(NormalizationMethod::Invalid).c_str(),    NormalizationMethod::Invalid)
        .export_values()
        ;
}
