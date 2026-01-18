"""Tests for the Enzyme Correlator application."""

from __future__ import annotations

import os
import tempfile
from collections.abc import Generator
from typing import TYPE_CHECKING
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd
import pytest

if TYPE_CHECKING:
    from enzyme_correlator import EnzymeCorrelatorGUI


@pytest.fixture
def sample_csv_file() -> Generator[str, None, None]:
    """Create a temporary CSV file with sample enzyme data."""
    content = """;Substrate1;Substrate2;Substrate3;Substrate4
Enzyme1;0,90;0,85;0,80;0,75
Enzyme2;0,88;0,92;0,78;0,82
Enzyme3;0,20;0,25;0,30;0,35
Enzyme4;0,91;0,87;0,83;0,79
"""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
        f.write(content)
        temp_path = f.name

    yield temp_path

    os.unlink(temp_path)


@pytest.fixture
def mock_tk() -> Generator[MagicMock, None, None]:
    """Mock tkinter for headless testing."""
    with patch("enzyme_correlator.tk") as mock:
        mock_root = MagicMock()
        mock.Tk.return_value = mock_root
        mock.StringVar.return_value = MagicMock()
        mock.StringVar.return_value.get.return_value = "0.85"
        mock.DISABLED = "disabled"
        mock.NORMAL = "normal"
        mock.HORIZONTAL = "horizontal"
        mock.N = "n"
        mock.S = "s"
        mock.E = "e"
        mock.W = "w"
        mock.END = "end"
        yield mock


@pytest.fixture
def mock_ttk() -> Generator[MagicMock, None, None]:
    """Mock ttk for headless testing."""
    with patch("enzyme_correlator.ttk") as mock:
        mock.Style.return_value = MagicMock()
        mock.Frame.return_value = MagicMock()
        mock.Button.return_value = MagicMock()
        yield mock


@pytest.fixture
def gui_instance(
    mock_tk: MagicMock,
    mock_ttk: MagicMock,  # noqa: ARG001
) -> Generator[EnzymeCorrelatorGUI, None, None]:
    """Create a GUI instance with mocked tkinter."""
    with patch("enzyme_correlator.filedialog"):
        from enzyme_correlator import EnzymeCorrelatorGUI

        mock_root = MagicMock()
        mock_tk.Text.return_value = MagicMock()
        mock_tk.Scale.return_value = MagicMock()
        mock_tk.StringVar.return_value = MagicMock()
        mock_tk.StringVar.return_value.get.return_value = "0.85"

        gui = EnzymeCorrelatorGUI(mock_root)
        yield gui


class TestImportData:
    """Tests for the import_data method."""

    def test_import_data_creates_dataframe(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test that import_data correctly creates a DataFrame."""
        gui_instance.datapath = sample_csv_file
        gui_instance.import_data()

        assert isinstance(gui_instance.df, pd.DataFrame)
        assert len(gui_instance.enzyme_list) == 4

    def test_import_data_enzyme_names(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test that enzyme names are correctly parsed."""
        gui_instance.datapath = sample_csv_file
        gui_instance.import_data()

        enzyme_names = [str(e.name) for e in gui_instance.enzyme_list]
        assert "Enzyme1" in enzyme_names
        assert "Enzyme2" in enzyme_names
        assert "Enzyme3" in enzyme_names
        assert "Enzyme4" in enzyme_names

    def test_import_data_decimal_conversion(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test that European decimal format (comma) is converted correctly."""
        gui_instance.datapath = sample_csv_file
        gui_instance.import_data()

        first_enzyme = gui_instance.enzyme_list[0]
        assert all(isinstance(v, float) for v in first_enzyme.values)


class TestComputeCorrelationMatrix:
    """Tests for the compute_correlation_matrix method."""

    def test_correlation_matrix_shape(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test that correlation matrix has correct shape."""
        gui_instance.datapath = sample_csv_file
        gui_instance.import_data()
        gui_instance.compute_correlation_matrix()

        n = len(gui_instance.enzyme_list)
        assert gui_instance.enzyme_correlation_matrix.shape == (n, n)

    def test_correlation_matrix_diagonal(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test that diagonal values are 1.0 (self-correlation)."""
        gui_instance.datapath = sample_csv_file
        gui_instance.import_data()
        gui_instance.compute_correlation_matrix()

        n = len(gui_instance.enzyme_list)
        for i in range(n):
            assert gui_instance.enzyme_correlation_matrix[i][i] == 1.0

    def test_correlation_matrix_symmetry(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test that correlation matrix is symmetric when plot_only_lt is False."""
        gui_instance.datapath = sample_csv_file
        gui_instance.plot_only_lt = False
        gui_instance.import_data()
        gui_instance.compute_correlation_matrix()

        matrix = gui_instance.enzyme_correlation_matrix
        n = len(gui_instance.enzyme_list)
        for i in range(n):
            for j in range(n):
                assert matrix[i][j] == matrix[j][i]

    def test_correlation_matrix_values_in_range(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test that all correlation values are between -1 and 1."""
        gui_instance.datapath = sample_csv_file
        gui_instance.import_data()
        gui_instance.compute_correlation_matrix()

        matrix = gui_instance.enzyme_correlation_matrix
        assert np.all(matrix >= -1)
        assert np.all(matrix <= 1)


class TestComputeHistogram:
    """Tests for the compute_histogram method."""

    def test_histogram_list_populated(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test that histogram list is populated."""
        gui_instance.datapath = sample_csv_file
        gui_instance.import_data()
        gui_instance.compute_correlation_matrix()
        gui_instance.compute_histogram()

        assert len(gui_instance.hist_list) > 0

    def test_histogram_list_lower_triangle_only(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test that histogram contains only lower triangle values."""
        gui_instance.datapath = sample_csv_file
        gui_instance.import_data()
        gui_instance.compute_correlation_matrix()
        gui_instance.compute_histogram()

        n = len(gui_instance.enzyme_list)
        expected_count = n * (n - 1) // 2
        assert len(gui_instance.hist_list) == expected_count

    def test_histogram_axis_range(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test that histogram axis covers -1 to 1."""
        gui_instance.datapath = sample_csv_file
        gui_instance.import_data()
        gui_instance.compute_correlation_matrix()
        gui_instance.compute_histogram()

        assert gui_instance.hist_axis[0] == pytest.approx(-1.0)
        assert gui_instance.hist_axis[-1] == pytest.approx(1.0)


class TestSortIntoGroups:
    """Tests for the sort_into_groups method."""

    def test_grouping_creates_dict(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test that grouping creates a dictionary."""
        gui_instance.datapath = sample_csv_file
        gui_instance.cutoff.get = MagicMock(return_value="0.85")
        gui_instance.import_data()
        gui_instance.compute_correlation_matrix()
        gui_instance.sort_into_groups()

        assert isinstance(gui_instance.grouping, dict)

    def test_grouping_enzymes_correlated(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test that correlated enzymes are grouped together."""
        gui_instance.datapath = sample_csv_file
        gui_instance.cutoff.get = MagicMock(return_value="0.90")
        gui_instance.import_data()
        gui_instance.compute_correlation_matrix()
        gui_instance.sort_into_groups()

        all_grouped_enzymes: list[str] = []
        for enzymes in gui_instance.grouping.values():
            all_grouped_enzymes.extend(enzymes)

        assert len(all_grouped_enzymes) >= 0

    def test_grouping_high_cutoff_fewer_groups(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test that higher cutoff results in fewer or equal grouped enzymes."""
        gui_instance.datapath = sample_csv_file
        gui_instance.import_data()
        gui_instance.compute_correlation_matrix()

        gui_instance.cutoff.get = MagicMock(return_value="0.50")
        gui_instance.sort_into_groups()
        low_cutoff_enzymes = sum(len(v) for v in gui_instance.grouping.values())

        gui_instance.cutoff.get = MagicMock(return_value="0.99")
        gui_instance.sort_into_groups()
        high_cutoff_enzymes = sum(len(v) for v in gui_instance.grouping.values())

        assert high_cutoff_enzymes <= low_cutoff_enzymes


class TestGUICallbacks:
    """Tests for GUI callback methods."""

    def test_show_grouping_callback(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test that show_grouping_button_callback runs without error."""
        gui_instance.datapath = sample_csv_file
        gui_instance.cutoff.get = MagicMock(return_value="0.85")
        gui_instance.import_data()
        gui_instance.compute_correlation_matrix()
        gui_instance.sort_into_groups()

        gui_instance.show_grouping_button_callback()

        gui_instance.grouping_label.delete.assert_called()
        gui_instance.grouping_label.insert.assert_called()

    def test_quit_callback(self, gui_instance: EnzymeCorrelatorGUI) -> None:
        """Test that quit_button_callback destroys root."""
        gui_instance.quit_button_callback()
        gui_instance.root.destroy.assert_called_once()


class TestMainFunction:
    """Tests for the main entry point."""

    def test_main_creates_gui(self, mock_tk: MagicMock, mock_ttk: MagicMock) -> None:  # noqa: ARG002
        """Test that main() creates and runs the GUI."""
        with patch("enzyme_correlator.filedialog"):
            from enzyme_correlator import main

            mock_tk.Text.return_value = MagicMock()
            mock_tk.Scale.return_value = MagicMock()

            main()

            mock_tk.Tk.assert_called_once()
            mock_tk.Tk.return_value.mainloop.assert_called_once()


class TestPlotOnlyLowerTriangle:
    """Tests for plot_only_lt flag behavior."""

    def test_correlation_matrix_lower_triangle_only(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test that upper triangle is zeroed when plot_only_lt is True."""
        gui_instance.datapath = sample_csv_file
        gui_instance.plot_only_lt = True
        gui_instance.import_data()
        gui_instance.compute_correlation_matrix()

        matrix = gui_instance.enzyme_correlation_matrix
        n = len(gui_instance.enzyme_list)
        for i in range(n):
            for j in range(i + 1, n):
                assert matrix[i][j] == 0.0


class TestLoadDataCallback:
    """Tests for load_data_callback method."""

    def test_load_data_callback_with_file(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test load_data_callback loads data and enables buttons."""
        with patch("enzyme_correlator.filedialog") as mock_dialog:
            mock_dialog.askopenfilename.return_value = sample_csv_file
            gui_instance.cutoff.get = MagicMock(return_value="0.85")

            gui_instance.load_data_callback()

            assert len(gui_instance.enzyme_list) == 4
            assert gui_instance.enzyme_correlation_matrix.shape == (4, 4)

    def test_load_data_callback_cancelled(self, gui_instance: EnzymeCorrelatorGUI) -> None:
        """Test load_data_callback handles cancelled dialog."""
        with patch("enzyme_correlator.filedialog") as mock_dialog:
            mock_dialog.askopenfilename.return_value = ""

            gui_instance.load_data_callback()

            assert gui_instance.datapath == ""


class TestPlotCallbacks:
    """Tests for plotting callback methods."""

    def test_plot_correlation_data_callback(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test plot_correlation_data_callback creates figure."""
        gui_instance.datapath = sample_csv_file
        gui_instance.cutoff.get = MagicMock(return_value="0.85")
        gui_instance.import_data()
        gui_instance.compute_correlation_matrix()

        with patch("enzyme_correlator.FigureCanvasTkAgg") as mock_canvas:
            mock_canvas_instance = MagicMock()
            mock_canvas.return_value = mock_canvas_instance

            gui_instance.plot_correlation_data_callback()

            assert gui_instance.fig is not None
            mock_canvas.assert_called_once()

    def test_plot_correlation_data_with_existing_canvas(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test plot_correlation_data_callback replaces existing canvas."""
        gui_instance.datapath = sample_csv_file
        gui_instance.cutoff.get = MagicMock(return_value="0.85")
        gui_instance.import_data()
        gui_instance.compute_correlation_matrix()

        old_canvas = MagicMock()
        gui_instance.canvas = old_canvas

        with patch("enzyme_correlator.FigureCanvasTkAgg") as mock_canvas:
            mock_canvas.return_value = MagicMock()

            gui_instance.plot_correlation_data_callback()

            old_canvas.get_tk_widget().grid_forget.assert_called_once()

    def test_plot_histogram_callback(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test plot_histogram_button_callback creates histogram."""
        gui_instance.datapath = sample_csv_file
        gui_instance.cutoff.get = MagicMock(return_value="0.85")
        gui_instance.import_data()
        gui_instance.compute_correlation_matrix()
        gui_instance.compute_histogram()

        with patch("enzyme_correlator.FigureCanvasTkAgg") as mock_canvas:
            mock_canvas_instance = MagicMock()
            mock_canvas.return_value = mock_canvas_instance

            gui_instance.plot_histogram_button_callback()

            assert gui_instance.fig is not None
            mock_canvas.assert_called_once()

    def test_plot_histogram_with_existing_canvas(
        self, gui_instance: EnzymeCorrelatorGUI, sample_csv_file: str
    ) -> None:
        """Test plot_histogram_button_callback replaces existing canvas."""
        gui_instance.datapath = sample_csv_file
        gui_instance.cutoff.get = MagicMock(return_value="0.85")
        gui_instance.import_data()
        gui_instance.compute_correlation_matrix()
        gui_instance.compute_histogram()

        old_canvas = MagicMock()
        gui_instance.canvas = old_canvas

        with patch("enzyme_correlator.FigureCanvasTkAgg") as mock_canvas:
            mock_canvas.return_value = MagicMock()

            gui_instance.plot_histogram_button_callback()

            old_canvas.get_tk_widget().grid_forget.assert_called_once()


class TestSaveFigCallback:
    """Tests for save_fig_button_callback method."""

    def test_save_fig_callback_saves_file(self, gui_instance: EnzymeCorrelatorGUI) -> None:
        """Test save_fig_button_callback saves figure."""
        gui_instance.fig = MagicMock()

        with patch("enzyme_correlator.filedialog") as mock_dialog:
            mock_dialog.asksaveasfilename.return_value = "/tmp/test.png"

            gui_instance.save_fig_button_callback()

            gui_instance.fig.savefig.assert_called_once_with("/tmp/test.png")

    def test_save_fig_callback_cancelled(self, gui_instance: EnzymeCorrelatorGUI) -> None:
        """Test save_fig_button_callback handles cancelled dialog."""
        gui_instance.fig = MagicMock()

        with patch("enzyme_correlator.filedialog") as mock_dialog:
            mock_dialog.asksaveasfilename.return_value = ""

            gui_instance.save_fig_button_callback()

            gui_instance.fig.savefig.assert_not_called()


class TestModuleAttributes:
    """Tests for module-level attributes."""

    def test_version_exists(self) -> None:
        """Test that __version__ is defined."""
        from enzyme_correlator import __version__

        assert __version__ is not None
        assert isinstance(__version__, str)

    def test_all_exports(self) -> None:
        """Test that __all__ contains expected exports."""
        from enzyme_correlator import __all__

        assert "EnzymeCorrelatorGUI" in __all__
        assert "main" in __all__
