"""The_Method, ported from Foundation.ipynb.

calculate_areas_and_isotopes is the one method in this whole codebase that scales
with row count instead of just number_of_peaks, so its bin-summation loop is
the part that now calls into the spectral_kernels Rust extension
(src/rust/spectral_kernels) instead of looping in Python. Everything else here is
unchanged bookkeeping over small lists.

The per-bin "Warning: Adjusting bin size..." / "Warning: bin_data is empty..."
prints fired from inside that loop and are gone now that it's in Rust; areas and
masses are unaffected/function same as the notebook :).
"""

import numpy as np
import pandas as pd

from .mass_spec_data import Mass_Spectrometer_Data

try:
    import spectral_kernels
except ImportError as _e:
    spectral_kernels = None
    _spectral_kernels_import_error = _e


class The_Method(Mass_Spectrometer_Data):
    """Draws comparisions among two different files of different NCE only"""

    def __init__(self):
        self.raw_data_1 = None
        self.raw_data_2 = None
        self.isotopes_1 = []
        self.isotopes_2 = []  # You get Isotopes from the parent, but I included this in case it's desired.
        self.masses_P = []
        self.masses_D = []
        self.areas_1 = []
        self.areas_2 = []
        self.ratios_1 = []
        self.ratios_2 = []
        self.alphas = []
        self.enrichments = []
        self.YN_list = []
        self.preference_list = []
        self.print_information = False

    def copy_from(self, Mass_Spectrometer_Data_instance):
        # Copy attributes from another Mass_Spectrometer_Data instance, assuming you want to use The_Method with the files from Mass_Spectrometer_Data.
        self.molecule_name = Mass_Spectrometer_Data_instance.molecule_name
        self.atom_name = Mass_Spectrometer_Data_instance.atom_name
        self.date_data_was_taken = Mass_Spectrometer_Data_instance.date_data_was_taken
        self.molecular_weight_of_molecule = Mass_Spectrometer_Data_instance.molecular_weight_of_molecule
        self.AMU_Element = Mass_Spectrometer_Data_instance.AMU_Element
        self.all_files_path = Mass_Spectrometer_Data_instance.all_files_path
        self.info_on_files_path_or_URL = Mass_Spectrometer_Data_instance.info_on_files_path_or_URL
        self.info_on_files_data = Mass_Spectrometer_Data_instance.info_on_files_data
        self.file_names = Mass_Spectrometer_Data_instance.file_names
        self.NCE = Mass_Spectrometer_Data_instance.NCE
        self.NCE_Labels_and_File_Numbers = Mass_Spectrometer_Data_instance.NCE_Labels_and_Files
        self.numeric_file_indices = Mass_Spectrometer_Data_instance.numeric_file_indices
        self.NCE_File_Paths = Mass_Spectrometer_Data_instance.NCE_File_Paths
        self.print_information = Mass_Spectrometer_Data_instance.print_information
        self.errors = Mass_Spectrometer_Data_instance.errors
        self.warnings = Mass_Spectrometer_Data_instance.warnings
        self.notices = Mass_Spectrometer_Data_instance.notices
        self.NCE_dataframes = Mass_Spectrometer_Data_instance.NCE_dataframes

    def load_data(self, csv_file_1, csv_file_2, row_skip):
        """Load the data from the two CSV files you are going to compare into the Pandas Dataframe. Takes in two strings of the file_paths for the files and an integer for how many rows to skip. This is for the event that you don't
        want to use the parent class Mass_Spectrometer_Data."""
        self.raw_data_1 = pd.read_csv(csv_file_1, skiprows=row_skip)  # If the file has text other than the data, put the number of rows needed to skip to get to the data. Ensure all files are uniform
        self.raw_data_2 = pd.read_csv(csv_file_2, skiprows=row_skip)

    def calculate_areas_and_isotopes(self, row_start, row_end, bin_size, file_number, ROWSSKIPPED=7, Calculate_Isotopes_and_Masses=False, rows_are_masses=False):
        """This will calculate areas and isotopes from the mass spectrometer data given the specifications. bin = int, row_start = int, row_end = int where the bin is a group of values you want to group together for each peak,
        row_start and row_end are respectively where you'll start grouping data and where the last row will be. It is inclusive of the value you put in. Finally, you need to put which file you're calculating
        the areas for, file = 1 or 2 for the respective files since the parameters of each is different. One note is that your rows should correspond to the row number in either excel or google sheets after deleting any rows before the headers and data. The indexing and reported bins values are specifically tailored for this.
        Calculate Isotopes = False if you want the function to calculate areas only, set it to True if you need the isotope values. If you index by mass, ROWSSKIPPED is needed to specify how
        many rows you'd skip in the original spreadsheet to get to the data and rows_are_mass must = True"""

        data = self.raw_data_1 if file_number == 1 else self.raw_data_2

        pandas_difference_index = 1  # Pandas uses 0-Based indexing whereas google sheets uses 1-Based indexing, so the answer would be one index off. This accounts for that
        header_affect_on_index = 1  # The index is also offset by the header, so this accounts for it being off by 1

        if rows_are_masses:
            masses_to_find = [row_start, row_end]

            # Define tolerance for approximate matching due to floating decimals between files. Change if necessary for more precise matching
            if "." in str(masses_to_find[0]):
                decimal_part = str(masses_to_find[0]).split(".")[1]
                if len(decimal_part) >= 5:
                    tolerance = 1e-5
                else:
                    n = len(decimal_part)
                    tolerance = 10 ** -n
            else:
                tolerance = 1e-5

            for m in masses_to_find:
                print(f"Mass {m}: Matches: {data['Mass'].apply(lambda x: abs(x - m) < tolerance).sum()}")

            row_indice = data[data['Mass'].apply(lambda x: any(abs(x - m) < tolerance for m in masses_to_find))].index.tolist()
            row_indices = []
            for index in row_indice:
                row_indices.append(index + ROWSSKIPPED + header_affect_on_index + pandas_difference_index)
            print(f"rows_indicies: {row_indices}")  # Debugger
            row_start = row_indices[0]
            row_end = row_indices[1]
            print(f"row_start: {row_start}")
            print(f"row_end: {row_end}")

        # Ensure row_start and row_end are within valid bounds
        if row_start < 0 or row_end >= len(data) or row_start > row_end:
            print(f"Invalid row range: row_start={row_start}, row_end={row_end}, max_row={len(data)}") if self.warnings else None
            return []

        if spectral_kernels is None:
            raise ImportError(
                "spectral_kernels Rust extension isn't built. From "
                "src/rust/spectral_kernels, run `pip install maturin && maturin develop`."
            ) from _spectral_kernels_import_error

        intensity = data['Intensity'].to_numpy(dtype=np.float64)
        mass_column = data['Mass'].to_numpy(dtype=np.float64) if Calculate_Isotopes_and_Masses else None

        areas_arr, masses_arr = spectral_kernels.bin_areas(
            intensity, row_start, row_end, bin_size,
            mass=mass_column,
            pandas_difference_index=pandas_difference_index,
            header_affect_on_index=header_affect_on_index,
        )
        areas = areas_arr.tolist()  # List to hold area sums
        masses = masses_arr.tolist() if Calculate_Isotopes_and_Masses else None  # Lists to hold isotope and masses values
        isotopes = [round(self.AMU_Element - (self.molecular_weight_of_molecule - x)) for x in masses] if Calculate_Isotopes_and_Masses else []

        # Store areas and isotopes in class attributes
        if file_number == 1:
            self.areas_1 = areas
            self.isotopes_1 = isotopes if Calculate_Isotopes_and_Masses else None
            self.masses_P = masses if Calculate_Isotopes_and_Masses else None
        else:
            self.areas_2 = areas
            self.isotopes_2 = isotopes if Calculate_Isotopes_and_Masses else None
            self.masses_D = masses if Calculate_Isotopes_and_Masses else None

        if self.print_information:
            print(f"Calculated Areas for file {file_number}: {areas}")
            print(f"Calculated Isotopes values for file {file_number}: {isotopes}") if Calculate_Isotopes_and_Masses else None
            print(f"Calculated Masses values for file {file_number}: {masses}") if Calculate_Isotopes_and_Masses else None
        return areas  # Returns the areas list, the isotope values are just stored.

        # For reference, the molar mass of Nd(NO3)4 is 392.2616 g/mol.

    def add_areas(self, row_start, row_end, bin_size, file_number, Calculate_Isotopes_and_Masses=False):
        "In some event where you have data files you need to add up but the peaks themselves are not in the same row index, you can use this function to use calculate_areas_and_isotopes in the exact same manner, except for different row start and row end"
        Areas = np.array(self.areas_1 if file_number == 1 else self.areas_2)
        self.calculate_areas_and_isotopes(row_start, row_end, bin_size, file_number, Calculate_Isotopes_and_Masses)
        New_areas = np.array(self.areas_1 if file_number == 1 else self.areas_2)
        Calculate_New_Areas = Areas + New_areas
        if file_number == 1:
            self.areas_1 = Calculate_New_Areas.tolist()
        else:
            self.areas_2 = Calculate_New_Areas.tolist()
        print(f"Added Areas: {Calculate_New_Areas.tolist()}") if self.print_information else None

        self.calculate_areas_and_isotopes(row_start, row_end, bin_size, file_number, Calculate_Isotopes_and_Masses)

    def calculate_ratio(self, Normalization_Number, *areas, file_number):
        """Calculates either the ratio of either a list of areas, several individual areas, or one individual area as inputted from *areas. The number will either assign the results to ratios_1 or ratios_2,
        Normalization_Number is an integer of the number you'll be using to normalize the data for the ratios as an integer value and file_number is also an integer value as to which file you want to calculate the ratios for."""
        if file_number != 1 and file_number != 2:
            print("An error will occur, you need to make the input 'number' either 1 or 2 to store the information in self.ratios_1 or self.ratios_2") if self.warnings else None
        else:
            ratio_lst = []
            if len(areas) == 1 and isinstance(areas[0], list):
                for area in areas[0]:  # loops over the list inside areas[0] which is needed since the list in *args will become a tuple
                    rat = area / Normalization_Number
                    ratio_lst.append(rat)  # Append each ratio between two peaks to the list
                if file_number == 1:
                    self.ratios_1 = ratio_lst
                else:
                    self.ratios_2 = ratio_lst
            else:
                for area in areas:
                    rat = area / Normalization_Number
                    ratio_lst.append(rat)
                if file_number == 1:
                    self.ratios_1 = ratio_lst
                else:
                    self.ratios_2 = ratio_lst
            print(f"Calculated Ratios for file {file_number}: {ratio_lst}") if self.print_information else None
            return ratio_lst

    # Test cases for clarity as to how calculate_ratio works:
    # 1.) For a list of areas: print(calculate_ratio(10, [20, 30, 40])) -> [2.0, 3.0, 4.0]
    # 2.) for Individual areas: print(calculate_ratio(10, 20, 30, 40)) -> [2.0, 3.0, 4.0]
    # 3.) For a single area: print(calculate_ratio(10, 20)) -> [2.0] which is the singular ratio

    def calculate_alpha(self, ratio_1, ratio_2):
        """Calculate the alpha of a specific ratio between the two files calculated ratios"""
        return ratio_2 / ratio_1 if ratio_1 and ratio_2 is not None else None

    def make_alphas_list(self):
        if len(self.ratios_1) != len(self.ratios_2):
            print("Type Error: The amount of ratios from the two files are not equal, you should use calculate alpha to manually assign the peaks you want to compare") if self.warnings else None
        alphas_lst = []
        list_length_differences = max(len(self.ratios_1), len(self.ratios_2)) - min(len(self.ratios_1), len(self.ratios_2))
        for i in range(min(len(self.ratios_1), len(self.ratios_2))):
            current_alpha = self.ratios_2[i] / self.ratios_1[i]
            alphas_lst.append(current_alpha)
        if list_length_differences != 0:
            while list_length_differences != 0:  # This ensures it won't error even with unequal list sizes.
                list_length_differences = list_length_differences - 1
                alphas_lst.append(None)
        self.alphas = alphas_lst
        print(f"Calculated Alphas: {alphas_lst}") if self.print_information else None
        return self.alphas

    def calculate_enrichment(self, alpha):
        return 1 - alpha

    def enrichment_list(self):
        enrichment_list = [self.calculate_enrichment(alpha) for alpha in self.alphas]
        self.enrichments = enrichment_list
        print(f"Calculated Enrichment: {enrichment_list}") if self.print_information else None
        return enrichment_list

    def isotope_effect(self):
        YN_lst = self.YN_list
        for alpha in self.alphas:
            if alpha == 1:
                YN_lst.append('No')
            else:
                YN_lst.append('Yes')
        self.YN_list = YN_lst
        print(f"Calculated YN: {YN_lst}") if self.print_information else None
        return YN_lst

    def preferences(self):
        preference_lst = self.preference_list
        for enrichment in self.enrichments:
            if enrichment > 1:
                preference_lst.append("D")  # preference for daughter species
            elif enrichment < 1:
                preference_lst.append("P")  # preference for parent species
            else:
                preference_lst.append("=")  # equal preference
        self.preference_list = preference_lst
        print(f"Calculated Preferences: {preference_lst}") if self.print_information else None
        return self.preference_list

    def drop(self, edit=False, *what_indexes, **what_list):
        """Drop specific values from the self.class attributes.

        Args:
            what_indexes (int): Indexes to drop from the list.
            edit (bool): If True, modifies the actual attribute by removing values.
            what_list (str): The name of the attribute to drop from.
        """
        for list_name, _ in what_list.items():
            if hasattr(self, list_name):
                target_list = getattr(self, list_name)
                if edit:
                    for i in sorted(what_indexes, reverse=True):  # Sort to avoid index shift
                        if i < len(target_list):
                            target_list.pop(i)
                else:
                    print([target_list[i] for i in what_indexes if i < len(target_list)]) if self.print_information else None
            else:
                print(f"Attribute '{list_name}' not found in the class.") if self.notices else None
