import pytest

import pyEM2

def test_add_genes(exp_mat):
    """Test that add_genes and genes work correctly together."""
    exp_mat.add_genes(["gene1", "gene2", "gene3"])

    assert len(exp_mat.genes) == 3
    assert exp_mat.genes == ["gene1", "gene2", "gene3"]
