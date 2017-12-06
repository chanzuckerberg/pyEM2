import numpy
import os
import pytest
import scipy.sparse

import pyEM2

@pytest.fixture
def similarity_em(tmpdir):
    """Create test matrix where odd columns and even columns are near each other:

    [[ 0, 11,  1,  7,  0, 10],
     [ 5,  0,  7,  0,  6,  1],
     [10,  0,  9,  2, 12,  1],
     [ 8,  0,  8,  0,  7,  1],
     [ 0, 14,  2, 17,  0, 15],
     [ 9,  0,  8,  0, 11,  1]]
    """
    dir_path = os.path.join(str(tmpdir), "EM2")

    mtx = numpy.matrix([[ 0, 11,  1,  7,  0, 10],
                        [ 5,  0,  7,  0,  6,  1],
                        [10,  0,  9,  2, 12,  1],
                        [ 8,  0,  8,  0,  7,  1],
                        [ 0, 14,  2, 17,  0, 15],
                        [ 9,  0,  8,  0, 11,  1]])
    csc_mtx = scipy.sparse.csc_matrix(mtx)
    em = pyEM2.ExpressionMatrix.from_csc_matrix(csc_mtx, dir_path)
    return em

def test_finds_correct_pairs(similarity_em):
    """Test that find_similar_pairs actually reports similar pairs.

    This is just a smoke test really, checking that the basics of the wrapper
    are working as expected.
    """
    sps = similarity_em.create_similar_pairs(0, 2, 1024, 42)

    # Odds go together, evens go together
    for i in range(len(similarity_em.cells)):
        similar_cells = sps.get_similar_cells(i)
        expected_pairs = [k for k in range(len(similarity_em.cells))
                          if k % 2 == i % 2 and k != i]
        assert len(similar_cells) == 2
        assert expected_pairs[0] in similar_cells
        assert expected_pairs[1] in similar_cells

def test_find_all_pairs(similarity_em):
    """Test that find_similar_pairs actually reports similar pairs.

    This is just a smoke test really, checking that the basics of the wrapper
    are working as expected.
    """
    sps = similarity_em.create_similar_pairs(0, 2, 1024, 42)
    nearest_neighbors = sps.get_all_similar_cells()

    assert nearest_neighbors.shape == (len(similarity_em.cells), 2)
    
    for i in range(len(similarity_em.cells)):
        similar_cells = nearest_neighbors[i]
        expected_pairs = [k for k in range(len(similarity_em.cells))
                          if k % 2 == i % 2 and k != i]
        assert len(similar_cells) == 2
        assert expected_pairs[0] in similar_cells
        assert expected_pairs[1] in similar_cells

def test_similar_pairs_cell_set(similarity_em):
    """Test defining a cell set for similar pairs."""

    # Cells have dummy names here, cell0, cell1, ...
    cell_set_1 = filter(lambda c: int(c.name[-1]) < 4, similarity_em.cells)
    cell_set_2 = [c for c in similarity_em.cells if int(c.name[-1]) > 0]

    sps_1 = similarity_em.create_similar_pairs(0, 2, 1024, 42, cell_set=cell_set_1)
    sps_2 = similarity_em.create_similar_pairs(0, 2, 1024, 42, cell_set=cell_set_2)

    # sps_1 has cell0, cell1, cell2, cell3. So the even odd rule should still hold
    assert sps_1.get_similar_cells(0) == [2]
    assert sps_1.get_similar_cells(1) == [3]
    assert sps_1.get_similar_cells(2) == [0]
    assert sps_1.get_similar_cells(3) == [1]

    # sps_2 has cell{1-5}. But, the even odd rule should still hold because we're
    # expecting global indices, not relative to the cell set

    # This guy segfaults right now because cell0 is not in the cell set...that's bad
    # sps_2.get_similar_cells(0))
    assert set(sps_2.get_similar_cells(1)) == set([3, 5])
    assert sps_2.get_similar_cells(2) == [4]
    assert set(sps_2.get_similar_cells(3)) == set([1, 5])
    assert sps_2.get_similar_cells(4) == [2]
    assert set(sps_2.get_similar_cells(5)) == set([1, 3])

def test_similar_pairs_gene_set(similarity_em):
    """Test defining a gene set for similar pairs."""

    # Subsetting genes shouldn't change the results in this case.
    gene_set_1 = reversed([g for g in similarity_em.genes if int(g.name[-1]) in (0,2,4)])
    sps = similarity_em.create_similar_pairs(0, 2, 512, 42, gene_set=gene_set_1)
    for i in range(len(similarity_em.cells)):
        similar_cells = sps.get_similar_cells(i)
        expected_pairs = [k for k in range(len(similarity_em.cells))
                          if k % 2 == i % 2 and k != i]
        assert len(similar_cells) == 2
        assert expected_pairs[0] in similar_cells
        assert expected_pairs[1] in similar_cells

    gene_set_2 = (g for g in similarity_em.genes if int(g.name[-1]) in (1,3,5))
    sps = similarity_em.create_similar_pairs(0, 2, 512, 42, gene_set=gene_set_2)
    for i in range(len(similarity_em.cells)):
        similar_cells = sps.get_similar_cells(i)
        expected_pairs = [k for k in range(len(similarity_em.cells))
                          if k % 2 == i % 2 and k != i]
        assert len(similar_cells) == 2
        assert expected_pairs[0] in similar_cells
        assert expected_pairs[1] in similar_cells

