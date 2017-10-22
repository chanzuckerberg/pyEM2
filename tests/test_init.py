import os
import pytest

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
