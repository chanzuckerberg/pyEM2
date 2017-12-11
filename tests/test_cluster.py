import numpy
import os
import pytest
import scipy.sparse
import scipy.stats

import pyEM2

@pytest.fixture
def cluster_graph_em(tmpdir):
    """An EM fixture for testing the CellGraph interface."""

    dir_path = os.path.join(str(tmpdir), "EM2")
    mtx = scipy.sparse.random(300, 40, format="csc", density=0.05, dtype=numpy.float,
                              data_rvs=scipy.stats.expon(loc=2, scale=100).rvs)
    em = pyEM2.ExpressionMatrix.from_csc_matrix(mtx, dir_path)
    return em

def test_calculate_clusters(cluster_graph_em):

    sps = cluster_graph_em.create_similar_pairs(0, 10, 512, 42)
    cg = cluster_graph_em.create_cell_graph(0, 5, sps)

    clusters = cluster_graph_em.cluster(cg, min_cluster_size=1)
    
    assert clusters.shape == (40,)
    assert len(numpy.unique(clusters)) < 12
