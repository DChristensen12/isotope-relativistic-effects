""" This is a Refactored port of Foundation.ipynb (fall_2024/). Same public names as the
notebook, just split one class/function group per file instead of one giant cell --
see each module's docstring for what moved where and why.

The one behavioral change is calculate_areas_and_isotopes's bin-summation loop
(see peak_areas.py), which now calls into the spectral_kernels Rust extension. It's
a bug-for-bug port of the original loop, not a cleanup -- see
src/rust/spectral_kernels/src/lib.rs for the details on what that means.
"""

from .state import (
    clear_GV_PD,
    clearGV,
)
from . import state

from .formatting import (
    Get_chemical_formula_in_Latex,
    format_chemical_formula,
)
from .io_utils import (
    file_number,
    import_google_spreadsheet_file,
    export_dataframe,
)
from .mass_spec_data import Mass_Spectrometer_Data
from .peak_areas import The_Method
from .error_propagation import Error_Propagation
from .graphs import (
    Graphs,
    αVIsol,
    αVNCE,
    Combined_Data_Line_Plot_EVI,
    Combined_Data_Line_Plot_AVI,
    overlaid_line_plot_EVI,
    overlaid_line_plot_AVI,
)
from .pdf_export import MakePDF, Add_Combined_Data_Table_To_PDF
from .pipeline import all_calculations_and_graphs

__all__ = [
    "state",
    "clear_GV_PD",
    "clearGV",
    "Get_chemical_formula_in_Latex",
    "format_chemical_formula",
    "file_number",
    "import_google_spreadsheet_file",
    "export_dataframe",
    "Mass_Spectrometer_Data",
    "The_Method",
    "Error_Propagation",
    "Graphs",
    "αVIsol",
    "αVNCE",
    "Combined_Data_Line_Plot_EVI",
    "Combined_Data_Line_Plot_AVI",
    "overlaid_line_plot_EVI",
    "overlaid_line_plot_AVI",
    "MakePDF",
    "Add_Combined_Data_Table_To_PDF",
    "all_calculations_and_graphs",
]
