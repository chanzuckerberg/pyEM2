import numpy
import os
import pytest
import scipy.sparse
import scipy.stats

import pyEM2

@pytest.fixture
def cell_graph_em(tmpdir):
    """An EM fixture for testing the CellGraph interface."""

    dir_path = os.path.join(str(tmpdir), "EM2")
    mtx = scipy.sparse.random(300, 40, format="csc", density=0.05, dtype=numpy.float,
                              data_rvs=scipy.stats.expon(loc=2, scale=100).rvs)
    em = pyEM2.ExpressionMatrix.from_csc_matrix(mtx, dir_path)
    return em

def test_create_cell_graph(cell_graph_em):
    """Test creation of a CellGraph and a few properties."""

    sps = cell_graph_em.create_similar_pairs(0, 10, 512, 42)
    cg = cell_graph_em.create_cell_graph(0, 5, sps)

    assert len(cg.vertices) == 40
    assert len(cg.edges) > 40
    assert len(cg.edges) < 1000

    for edge in cg.edges:
        assert len(edge) == 2
        assert edge[0].name.startswith("cell")
        assert edge[1].name.startswith("cell")

    for vertex in cg.vertices:
        assert vertex.cell.name.startswith("cell")
        assert isinstance(vertex.x, float)
        assert isinstance(vertex.y, float)

def test_cell_graph_cell_set(cell_graph_em):

    cell_set_1 = (c for c in cell_graph_em.cells if int(c.name[4:]) % 2 == 0)
    cell_set_2 = [c for c in cell_graph_em.cells if int(c.name[4:]) % 2 == 1]

    sps = cell_graph_em.create_similar_pairs(0, 10, 512, 42)
    cg = cell_graph_em.create_cell_graph(0, 10, sps, cell_set_1)

    assert len(cg.vertices) == 20
    assert len(cg.edges) > 20

    sps = cell_graph_em.create_similar_pairs(0, 10, 512, 42, cell_set_2)
    cg = cell_graph_em.create_cell_graph(0, 10, sps)

    assert len(cg.vertices) == 20
