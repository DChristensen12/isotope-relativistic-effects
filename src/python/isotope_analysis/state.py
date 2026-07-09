# Global state, split out of the single Foundation.ipynb cell it used to live in.
#
# These are the same module-level dicts/lists/tuples the notebook had at the top of
# its one big code cell, just moved here so every other module can share them by
# importing this module instead of everything living in one notebook namespace.
# Other modules should do `from . import state` and read/write `state.X` rather than
# `from .state import X`, since the latter would just copy the reference at import
# time and not see later reassignments (e.g. after clearGV()/clear_GV_PD()).

# These dicts all share the same key format:
# {
#   'Isolation': {
#     'Isotope': {
#       'Parent': {
#         'Daughter': {
#           'Normalized Collision Energy Value': <value>
#         }
#       },
#       'Daughter': {
#         'Parent': {
#           'Normalized Collision Energy Value': <value>
#         }
#       }
#     }
#   }
# }
All_alphas_dict_PD = {}       # from all_calculations_and_graphs, updated every call
All_enrichments_dict_PD = {}  # from all_calculations_and_graphs, updated every call
All_Errors_dict_PD = {}       # from all_calculations_and_graphs, updated every call
All_Parent_Areas_dict_PD = {}   # from all_calculations_and_graphs, updated every call
All_Daughter_Areas_dict_PD = {} # from all_calculations_and_graphs, updated every call
Daughter_NCE_Values_PD = []   # from all_calculations_and_graphs, updated every call
Parent_NCE_Values_PD = []     # from all_calculations_and_graphs, updated every call

# These dicts all share the date the data was taken as their keys.
All_alphas = {}          # from Error_Propagation.calculate_δα_indv, updated every call
All_enrichments = {}     # from Error_Propagation.calculate_δα_indv, updated every call
All_Errors = {}          # from Error_Propagation.calculate_δα_indv, updated every call
All_Parent_Areas = {}    # from Error_Propagation.calculate_δα_indv, updated every call
All_Daughter_Areas = {}  # from Error_Propagation.calculate_δα_indv, updated every call
All_Isotopes = ()   # set once -- every file is assumed to share the same isotopes
Parent_masses = ()  # set once -- every file is assumed to share the same parent masses
Daughter_masses = () # set once -- every file is assumed to share the same daughter masses
Atom_Names = " "     # set once -- this code doesn't support graphing multiple elements at once


def clear_GV_PD():
    """Clears everything except what clearGV() handles. The _PD suffix marks
    Parent/Daughter comparison state, used for graphs that compare across multiple
    isolations once all of those isolations have already been analyzed -- pair this
    with clearGV() if you want a full reset."""
    global All_alphas_dict_PD, All_enrichments_dict_PD, All_Errors_dict_PD
    global All_Parent_Areas_dict_PD, All_Daughter_Areas_dict_PD
    global Daughter_NCE_Values_PD, Parent_NCE_Values_PD

    All_alphas_dict_PD = {}
    All_enrichments_dict_PD = {}
    All_Errors_dict_PD = {}
    All_Parent_Areas_dict_PD = {}
    All_Daughter_Areas_dict_PD = {}
    Daughter_NCE_Values_PD = []
    Parent_NCE_Values_PD = []


def clearGV():
    """Clears the state that clear_GV_PD() doesn't touch."""
    global All_alphas, All_enrichments, All_Errors, All_Parent_Areas, All_Daughter_Areas
    global All_Isotopes, Parent_masses, Daughter_masses, Atom_Names

    All_alphas = {}
    All_enrichments = {}
    All_Errors = {}
    All_Parent_Areas = {}
    All_Daughter_Areas = {}
    All_Isotopes = ()
    Parent_masses = ()
    Daughter_masses = ()
    Atom_Names = " "
