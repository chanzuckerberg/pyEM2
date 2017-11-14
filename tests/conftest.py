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

@pytest.fixture
def exp_mat_with_cells(exp_mat):

    exp_mat.add_cell(
        "cell1",
        {"CellColor": "blue", "Weight": "11.6", "Length": "3"},
        {"gene1": 0.0, "gene2": 1<<32, "gene4": 11})

    exp_mat.add_cell(
        "cell2",
        {"CellShade": "very light", "CellColor": "red", "Weight": "9", "Length": "7.6342"},
        {"gene1": 5, "gene3": 100.1, "gene4": 15, "gene5": 60.123})

    exp_mat.add_cell(
        "cell3",
        {"CellShade": "rather dark", "CellColor": "green", "Weight": "29", "Length": "642"},
        {"gene1": 15, "gene3": 10.1, "gene4": 1.5, "gene5": 6.23})

    return exp_mat
