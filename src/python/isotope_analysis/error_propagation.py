"""Error_Propagation, ported from Foundation.ipynb. Operates on arrays sized by
number_of_peaks (single digits), so nothing here scales with data size.

Global dict access (All_alphas, All_Isotopes, etc.) goes through the `state` module
instead of `global` statements, since those variables now live in a different file --
see state.py for why callers should go through `state.X` rather than importing the
names directly."""

import math

import numpy as np
import pandas as pd

from . import state
from .mass_spec_data import Mass_Spectrometer_Data
from .peak_areas import The_Method


class Error_Propagation(Mass_Spectrometer_Data):
    def __init__(self, Mass_Spectrometer_Data_Instance):
        self.Data = Mass_Spectrometer_Data_Instance
        self.The_Method_Function_Access = The_Method()
        self.The_Method_Function_Access.copy_from(self.Data)
        self.parent_total_areas, self.daughter_total_areas = [], []
        self.Error_D, self.Error_P = [], []
        self.δ_PA, self.δ_DA = [], []
        self.alphas = []
        self.enrichments = []
        self.δα = []
        self.Isotopes = []
        self.masses_P = []
        self.masses_D = []
        self.Parent_ratios = []
        self.Daughter_ratios = []
        self.print_information = False
        self.file_names_P = []
        self.file_names_D = []
        self.Preference_list = []

    def Individual_δ(self, *Areas):
        """Calculates the individual δ of every Area in the list or individual area depending what is passed in. This allows you to individually select Areas to find the error of."""
        if len(Areas) == 1 and isinstance(Areas[0], list):
            return np.sqrt(np.array(Areas[0])).tolist()
        else:
            return math.sqrt(Areas[0])

    def calculate_δα_indv(self, row_starts_and_ends_parent, row_starts_and_ends_daughter, Bin_size, number_of_peaks, normalization, parent_daughter, ROWSKIP=7, NBTA=False, rows_are_masses=False):
        """Calculates δα for each individual area as a method of getting error bars. Then combines them. This method allows you to keep track of each individual error before the Areas are added up and keep track of all data being
        produced up to the calculation of the errors.

        Specifications:
        row_starts_and_ends_parent = [row_start_parent, row_end_parent] if you have evenly spaced peaks, otherwise you can put [[row_start_parent_1, row_end_parent_1], [row_start_parent_2, row_end_parent_2]] and so on for any unevenly spaced peaks.
        row_starts_and_ends_daughter = same as row_starts_and_ends_parent, but for daughters
        Bin_Size = integer of how many values you want for each rows in the spreedsheet
        number_of_peaks = integer of how many peaks corresponding to the lighter molecules you found.
        normalization = integer corresponding to number you want to use to normalize the data to.
        parent_daughter = = ['parent_NCE_Value', 'daughter_NCE_Value'] to compare between two NCE (of consistent isolations) data.
        NBTA is True if you want to normalize by the total areas, otherwise it'll use the normalization integer you put.
        """

        # Uses the specified NCE values as keys
        parent = parent_daughter[0]
        daughter = parent_daughter[1]

        # Gathers the file paths for the specified keys. Checks if the specified parent_daughter pair works, if not it'll attempt to use the file_paths that would be gathered from seperate_files.
        Parent_Files = self.Data.NCE_File_Paths[parent]
        Daughter_Files = self.Data.NCE_File_Paths[daughter]

        # Create all variables to store information in
        parent_areas, parent_error, daughter_areas, daughter_error, parent_percent_error, daughter_percent_error, δ_PA, δ_DA, δα = ([] for _ in range(9))
        xx, yy = np.zeros(number_of_peaks), np.zeros(number_of_peaks)  # Total parent Areas array, Total daughter Areas array

        # To count through the file names for labeling purposes
        file_tracker_D = -1
        file_tracker_P = -1
        file_number_P, file_number_D = self.Data.numeric_file_indices[:len(Parent_Files)], self.Data.numeric_file_indices[len(Parent_Files):]  # Actual number of file
        self.file_names_P = file_number_P
        self.file_names_D = file_number_D

        # Summing areas of the same peak areas across all specified files
        for file_P in Parent_Files:
            file_tracker_P += 1
            self.The_Method_Function_Access.raw_data_1 = pd.read_csv(file_P, skiprows=ROWSKIP)
            x = np.array([])
            for row_starts_and_ends in row_starts_and_ends_parent:
                start_P, end_P = row_starts_and_ends[0], row_starts_and_ends[1]
                if rows_are_masses:
                    area_list = self.The_Method_Function_Access.calculate_areas_and_isotopes(start_P, end_P, Bin_size, 1, ROWSSKIPPED=ROWSKIP, rows_are_masses=True)
                else:
                    area_list = self.The_Method_Function_Access.calculate_areas_and_isotopes(start_P - ROWSKIP, end_P - ROWSKIP, Bin_size, 1, ROWSSKIPPED=ROWSKIP)
                x = np.concatenate([x, area_list])
            print(f"Areas for P file {file_number_P[file_tracker_P]}: {x}") if self.print_information else None
            print(f" δ_PA file {file_number_P[file_tracker_P]}: {np.sqrt(x)}") if self.print_information else None
            xx = xx + x
        for file_D in Daughter_Files:
            file_tracker_D += 1
            self.The_Method_Function_Access.raw_data_2 = pd.read_csv(file_D, skiprows=ROWSKIP)
            y = np.array([])
            for row_starts_and_end in row_starts_and_ends_daughter:
                start, end = row_starts_and_end[0], row_starts_and_end[1]
                if rows_are_masses:
                    area_lst = self.The_Method_Function_Access.calculate_areas_and_isotopes(start, end, Bin_size, 2, ROWSSKIPPED=ROWSKIP, rows_are_masses=True)
                else:
                    area_lst = self.The_Method_Function_Access.calculate_areas_and_isotopes(start - ROWSKIP, end - ROWSKIP, Bin_size, 2, ROWSSKIPPED=ROWSKIP)
                y = np.concatenate([y, area_lst])
            print(f"Areas for D file {file_number_D[file_tracker_D]}: {y}") if self.print_information else None
            print(f" δ_DA file {file_number_D[file_tracker_D]}: {np.sqrt(y)}") if self.print_information else None
            yy = yy + y
        self.parent_total_areas = xx
        if self.Data.date_data_was_taken in state.All_Parent_Areas:
            key = self.Data.date_data_was_taken
            new_key = key
            # Looping to check if the key exists and to keep appending a prime until it's a unique key
            while new_key in state.All_Parent_Areas:
                new_key = new_key + "'"
            state.All_Parent_Areas[new_key] = tuple(xx)
        else:
            state.All_Parent_Areas[self.Data.date_data_was_taken] = tuple(xx)
        self.daughter_total_areas = yy
        if self.Data.date_data_was_taken in state.All_Daughter_Areas:
            key = self.Data.date_data_was_taken
            new_key = key
            while new_key in state.All_Daughter_Areas:
                new_key = new_key + "'"
            state.All_Daughter_Areas[new_key] = tuple(yy)
        else:
            state.All_Daughter_Areas[self.Data.date_data_was_taken] = tuple(yy)
        print(f"total areas for P: {xx}") if self.print_information else None
        print(f"total areas for D: {yy}") if self.print_information else None
        δ_PA = np.sqrt(xx)
        δ_DA = np.sqrt(yy)
        self.δ_PA = δ_PA
        self.δ_DA = δ_DA
        print(f"δ_PA: {δ_PA}") if self.print_information else None
        print(f"δ_DA: {δ_DA}") if self.print_information else None
        Error_P = δ_PA / xx
        Error_D = δ_DA / yy
        print(f"% error for P: {Error_P}") if self.print_information else None
        print(f"% error for D: {Error_D}") if self.print_information else None
        self.Error_D = Error_D
        self.Error_P = Error_P
        daughter_ratios = yy / np.sum(yy) if NBTA else yy / normalization
        parent_ratios = xx / np.sum(xx) if NBTA else xx / normalization
        print(f"Parent Ratios: {parent_ratios}") if self.print_information else None
        print(f"Daughter Ratios: {daughter_ratios}") if self.print_information else None
        self.Parent_ratios = parent_ratios
        self.Daughter_ratios = daughter_ratios

        alpha = daughter_ratios / parent_ratios
        self.alphas = alpha
        print(f"Alpha Values: {alpha}") if self.print_information else None
        if self.Data.date_data_was_taken in state.All_alphas:
            key = self.Data.date_data_was_taken
            new_key = key
            while new_key in state.All_alphas:
                new_key = new_key + "'"
            state.All_alphas[new_key] = tuple(alpha)
        else:
            state.All_alphas[self.Data.date_data_was_taken] = tuple(alpha)  # Puts the alpha values in the global variable that is a dictionary of all alphas corresponding to a key which is the date the data was taken.

        Enrichment = alpha - 1
        self.enrichments = Enrichment
        print(f"Enrichment Values: {Enrichment}") if self.print_information else None
        if self.Data.date_data_was_taken in state.All_enrichments:
            key = self.Data.date_data_was_taken
            new_key = key
            while new_key in state.All_enrichments:
                new_key = new_key + "'"
            state.All_enrichments[new_key] = tuple(Enrichment)
        else:
            state.All_enrichments[self.Data.date_data_was_taken] = tuple(Enrichment)  # Puts the alpha values in the global variable that is a dictionary of all alphas corresponding to a key which is the date the data was taken.

        δα = (alpha * np.sqrt((Error_D) ** 2 + (Error_P) ** 2))
        self.δα = δα
        print(f"δα: {δα}") if self.print_information else None
        if self.Data.date_data_was_taken in state.All_Errors:
            key = self.Data.date_data_was_taken
            new_key = key
            while new_key in state.All_Errors:
                new_key = new_key + "'"
            state.All_Errors[new_key] = tuple(δα)
        else:
            state.All_Errors[self.Data.date_data_was_taken] = tuple(δα)  # Puts the Error values in the global variable that is a dictionary of all Errors corresponding to a key which is the date the data was taken.

        # Getting the isotopes and masses
        self.The_Method_Function_Access.raw_data_1 = pd.read_csv(Parent_Files[0], skiprows=ROWSKIP)
        ISOTOPES = []
        PARENT_MASSES = []
        DAUGHTER_MASSES = []
        self.The_Method_Function_Access.raw_data_1 = pd.read_csv(Parent_Files[0], skiprows=ROWSKIP)
        self.The_Method_Function_Access.raw_data_2 = pd.read_csv(Daughter_Files[0], skiprows=ROWSKIP)
        for row_starts_and_ends in row_starts_and_ends_parent:
            start_P, end_P = row_starts_and_ends[0], row_starts_and_ends[1]
            if rows_are_masses:
                self.The_Method_Function_Access.calculate_areas_and_isotopes(start_P, end_P, Bin_size, 1, ROWSSKIPPED=ROWSKIP, Calculate_Isotopes_and_Masses=True, rows_are_masses=True)
            else:
                self.The_Method_Function_Access.calculate_areas_and_isotopes(start_P - ROWSKIP, end_P - ROWSKIP, Bin_size, 1, ROWSSKIPPED=ROWSKIP, Calculate_Isotopes_and_Masses=True)
            for i in self.The_Method_Function_Access.isotopes_1:
                ISOTOPES.append(i)
            for j in self.The_Method_Function_Access.masses_P:
                PARENT_MASSES.append(j)
        print(f"Isotopes: {ISOTOPES}") if self.print_information else None
        if state.All_Isotopes == ():  # This is to make the All_Isotopes global variable have the isotopes. This only occurs once when the tuple is empty.
            state.All_Isotopes = tuple(ISOTOPES)
        print(f"Masses Parent: {PARENT_MASSES}") if self.print_information else None
        for row_starts_and_ends in row_starts_and_ends_daughter:
            start_D, end_D = row_starts_and_ends[0], row_starts_and_ends[1]
            if rows_are_masses:
                self.The_Method_Function_Access.calculate_areas_and_isotopes(start_D, end_D, Bin_size, 2, ROWSSKIPPED=ROWSKIP, Calculate_Isotopes_and_Masses=True, rows_are_masses=True)
            else:
                self.The_Method_Function_Access.calculate_areas_and_isotopes(start_D - ROWSKIP, end_D - ROWSKIP, Bin_size, 2, ROWSSKIPPED=ROWSKIP, Calculate_Isotopes_and_Masses=True)
            for k in self.The_Method_Function_Access.masses_D:
                DAUGHTER_MASSES.append(k)
        print(f"Masses Daughter: {DAUGHTER_MASSES}") if self.print_information else None
        self.Isotopes = ISOTOPES
        self.masses_P = PARENT_MASSES
        self.masses_D = DAUGHTER_MASSES
        # For storing parent and daughter masses in the Parent_masses and Daughter_Masses global variables only once when they are empty.
        if state.Parent_masses == ():
            Parent_Masses = tuple(PARENT_MASSES)
        if state.Daughter_masses == ():
            Daughter_Masses = tuple(DAUGHTER_MASSES)

        if state.Atom_Names == " ":
            state.Atom_Names = self.The_Method_Function_Access.atom_name

        self.Preference_list = np.array(['Parent' if ε < 0 else 'Neither' if ε == 0 else 'Daughter' for ε in self.enrichments])
        # If you notice any data you want to remove, use the drop_files function in the parent class Mass_Spectrometer_Data, run the get_file_paths function again, then perform all calculations
