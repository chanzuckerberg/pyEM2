"""Test some memory management behavior."""


def test_del_em(exp_mat_with_cells):
    """Test that you can del the ExpressionMatrixWrapper and objects it has
    created remain valid.
    """

    sps = exp_mat_with_cells.create_similar_pairs(0, 2, 512, 42)

    # The underlying ExpressionMatrix should survive because SimilarPairsWrapper
    # holds a shared_ptr to it.
    del(exp_mat_with_cells)

    # So, this shouldn't segfault
    sps.get_similar_cells(0)
