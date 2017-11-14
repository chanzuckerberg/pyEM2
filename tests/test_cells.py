import pytest

import pyEM2

def test_add_cells(exp_mat_with_cells):
    """Test adding cells via add_cell."""

    cells = exp_mat_with_cells.cells
    assert len(cells) == 3

    names = [c.name for c in cells]
    assert names == ["cell1", "cell2", "cell3"]


def test_add_cells_with_same_name(exp_mat):
    """Test error when adding cell with same name."""
    exp_mat.add_cell("cell1", {}, {})

    with pytest.raises(RuntimeError) as excinfo:
        exp_mat.add_cell("cell1", {}, {})

    assert "already exists" in str(excinfo.value)

def test_edit_metadata(exp_mat_with_cells):
    """Test getting and setting metadata."""

    # Get the cell1 metadata
    cell1_metadata = exp_mat_with_cells.cells[0].metadata
    assert exp_mat_with_cells.cells[0].metadata["CellColor"] == "blue"

    # Update the metadata and write it back
    cell1_metadata["CellType"] = "Blood"
    cell1_metadata["CellColor"] = "yellow"
    exp_mat_with_cells.cells[0].metadata = cell1_metadata

    assert exp_mat_with_cells.cells[0].metadata["CellColor"] == "yellow"
    assert exp_mat_with_cells.cells[0].metadata["CellType"] == "Blood"

def test_get_counts(exp_mat_with_cells):
    """Test getting expression counts from cells."""

    # gene1 shouldn't be present because its expression is zero
    assert "gene1" not in exp_mat_with_cells.cells[0].expression_counts

    # But gene2 should be there
    assert exp_mat_with_cells.cells[0].expression_counts["gene2"] == 1<<32
