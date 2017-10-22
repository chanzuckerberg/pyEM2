import os
import pytest

import pyEM2

@pytest.fixture
def exp_mat(tmpdir):
    return pyEM2.ExpressionMatrix(
        os.path.join(str(tmpdir), "EM2"),
        1000,
        1000,
        1000,
        1000)

