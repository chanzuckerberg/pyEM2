import numpy
import os
import pytest
import scipy.sparse

import pyEM2

def test_init(tmpdir):
    """Test initialization of EM with all args in order."""
    exp_mat = pyEM2.ExpressionMatrix(
        os.path.join(str(tmpdir), "EM2"),
        1000,
        1000,
        1000,
        1000)

    assert isinstance(exp_mat, pyEM2.ExpressionMatrix)


def test_init_no_optionals(tmpdir):
    """Test initialization of EM with no optional args."""
    exp_mat = pyEM2.ExpressionMatrix(
        os.path.join(str(tmpdir), "EM2"))

    assert isinstance(exp_mat, pyEM2.ExpressionMatrix)


def test_init_named_args(tmpdir):
    """Test initializtion of EM with out of order, names kwargs."""
    exp_mat = pyEM2.ExpressionMatrix(
        os.path.join(str(tmpdir), "EM2"),
        cell_metadata_value_capacity=1000,
        cell_metadata_name_capacity=1000,
        cell_capacity=1000,
        gene_capacity=1000)

    assert isinstance(exp_mat, pyEM2.ExpressionMatrix)

def test_init_from_existing(tmpdir):
    """Test initialization of EM using the from_existing_directory
    staticmethod."""
    dir_path = os.path.join(str(tmpdir), "EM2")
    orig_exp_mat = pyEM2.ExpressionMatrix(dir_path)
    assert isinstance(orig_exp_mat, pyEM2.ExpressionMatrix)
    del orig_exp_mat

    new_exp_mat = pyEM2.ExpressionMatrix.from_existing_directory(dir_path)
    assert isinstance(new_exp_mat, pyEM2.ExpressionMatrix)

def test_init_from_existing_fails_on_missing(tmpdir):
    """Test that from_existing_directory fails when the directory
    doesn't exist.
    """
    dir_path = os.path.join(str(tmpdir), "EM2")

    with pytest.raises(RuntimeError) as excinfo:
        new_exp_mat = pyEM2.ExpressionMatrix.from_existing_directory(dir_path)

    assert "No such file" in str(excinfo.value)

def test_init_from_csc_matrix(tmpdir):
    """Test initialization from a scipy csc matrix."""

    dir_path = os.path.join(str(tmpdir), "EM2")
    row = numpy.array([0, 2, 2, 0, 1, 2])
    col = numpy.array([0, 0, 1, 2, 2, 2])
    data = numpy.array([1, 2, 3, 4, 5, 6])
    csc_mtx = scipy.sparse.csc_matrix((data, (row, col)), shape=(3, 4))

    em = pyEM2.ExpressionMatrix.from_csc_matrix(csc_mtx, dir_path)

    assert set(em.genes) == set(["gene0", "gene1", "gene2"])
