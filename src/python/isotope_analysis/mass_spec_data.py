"""Mass_Spectrometer_Data, ported from Foundation.ipynb. Logic is unchanged from the
notebook version, since this class is all file-path/dict bookkeeping over the file list,
not the spectral data itself."""

import os
import pandas as pd

from .formatting import format_chemical_formula
from .io_utils import import_google_spreadsheet_file


class Mass_Spectrometer_Data:
    """Allows you to connect all the classes and data_files"""

    def __init__(self):
        self.molecule_name = "[Dy(NO_{3})_{4}]^{-}"  # Example, the format doesn't matter, I just prefer it this way.
        self.atom_name = 'Dy'
        self.date_data_was_taken = "11/14/1987"
        self.molecular_weight_of_molecule = 410.516  # g/mol. This is just an example
        self.AMU_Element = 162.5
        self.all_files_path = "path"  # A string of the Path to folder containing every CSV file
        self.info_on_files_path_or_URL = "pathh"  # A string of the path to a CSV file of all info about the files. This is specific to the project. If google sheets,this should be a URL.
        self.info_on_files_data = None
        self.file_names = []
        self.NCE = []
        self.Isolations = []
        self.NCE_Labels_and_Files = {}  # A dictionary that shows which files belong to which NCE values. This does not specify isolation.
        self.numeric_file_indices = []
        self.NCE_File_Paths = {}  # A dictionary containing the key of the NCE values and the value pair of a list of file path names corresponding to each file.
        self.print_information = False  # Set to True if you want to print out relevant information you have. This does not include error messages.
        self.errors = True  # Set to False if you don't want to print out any Errors.
        self.warnings = True  # Set to False if you don't want to print out any warnings.
        self.notices = True  # Set to False if you don't want to print out any notices
        self.NCE_dataframes = {}  # For the event that you want to just store all data into a data frame
        self.nce_id = 'ISOUPTOPE'  # A list of the unique NCE IDs
        self.The_Lorax = []  # Has all the unique NCE IDs in the order they were found
        self.unsorted_NCE = []  # A list of the NCE value for every file.
        self.unsorted_Isolations = []  # A list of the Isolation for every file.
        self.unsorted_NCE_ID = []  # A list of the NCE ID for every file.
        self.numeric_file_indices_dict = {}  # Puts the number of the files

    def define_molecule_characteristics(self, molecule_name, name_of_atom_of_interest, date_data_was_taken, molecular_weight, AMU):
        """Takes the of the molecule that you analyzed in the Mass spectrometer as a string, the name of the specific atom of the molecule that is of interest/all ligands are
        bound to as a string, the date of the analysis in the format of MM/DD/YEAR as a string, the molecular weight of the molecule as a integer, and the AMU of the element
        (such as Dy or Nd) from the periodic table as a integer."""
        try:
            self.molecule_name = format_chemical_formula(molecule_name)
        except:  # noqa: E722
            self.molecule_name = molecule_name
        self.atom_name = name_of_atom_of_interest
        self.date_data_was_taken = date_data_was_taken
        self.molecular_weight_of_molecule = molecular_weight
        self.AMU_Element = AMU

    def define_all_file_path(self, file_path):
        """Takes in a string of the file_path that contains every single file you want to use for data analysis. This path should contain every single CSV file you want to use"""
        self.all_files_path = file_path

    def define_NCE(self, info_on_files_path_or_URL, **Manually):
        """For this you can either put in a string of the file_path of a file with every NCE labled in a column with the associated file numbers
        and let the code sort it out or add an additional arguments Manual = True and NCE_values = [put values in here]. The intended file_path has columns labeled
        file number, isolation, NCE, and comment. file number is the number of the associated file such as file 0, file 1, etc. Isolation is the associated isolation
        for that file in the format of isolated mass(bin width). NCE is the associated normalized collision energy used in the mass spectrometer. The comment is not
        needed for the code, but you should examine it to deicde whether to that file later on in events such as running out of spray or forgetting to change the range you used;
        the comment column is not needed for the code, but you could have one to decide whether the data isn't good to use for statisitcal analysis or not."""

        CSV, Excel, GoogleSheets, Manual = (Manually.get(key, False) for key in ['CSV', 'Excel', 'GoogleSheets', 'Manual'])  # Extracts optional parameters or set default False values
        self.info_on_files_path_or_URL = info_on_files_path_or_URL

        # Load the appropriate file type
        if CSV:
            self.info_on_files_data = pd.read_csv(info_on_files_path_or_URL)
        elif Excel:
            self.info_on_files_data = pd.read_excel(info_on_files_path_or_URL)
        elif GoogleSheets:
            self.info_on_files_data = import_google_spreadsheet_file(info_on_files_path_or_URL)
        else:
            print('You need to specify what type of file it is (CSV, Excel, or GoogleSheets)') if self.errors else None
            return  # Exit the function if file type is not specified
        NCE_values = []
        if Manual:  # If Manual is True, retrieve the provided NCE values
            NCE_values = Manually.get('NCE_values', [])
        else:
            # Otherwise, extract NCE values from the data
            data = self.info_on_files_data
            nce_columns = {'NCE', 'nce', 'nce ', ' nce', 'NCE ', ' NCE'}
            file_columns = {'file', 'file ', ' file', 'file #', 'file#', 'file number', 'File', 'file_number', 'File_Number', 'File_number', 'File_#', 'File#', 'File #', 'File ', ' File'}
            isolations_and_width = {'Isolation', 'isolation', ' isolation', ' isolation', ' Isolation', 'Isolation '}
            data.columns = ['NCE' if col in nce_columns else 'file number' if col in file_columns else 'isolation' if col in isolations_and_width else col
                            for col in data.columns]

            if not any(col in data.columns for col in nce_columns) or not any(col in data.columns for col in file_columns) or not any(col in data.columns for col in isolations_and_width):
                print("Error: 'NCE' or 'file number' or 'isolation' column missing from data.") if self.error else None
                return {'>:('}

            self.info_on_files_data = data  # To save the changed names for consistency
            if 'NCE' in data.columns:
                x = data['NCE'].tolist()
                y = data['isolation'].tolist()
                self.unsorted_NCE = x  # To store NCE for each file so it can be used for indexing in different functions
                self.unsorted_Isolations = y  # To store NCE for each file so it can be used for indexing later in different functions
                nce_id = [lambda t=t: str(x[t]) + "_" + str(y[t]) for t in range(len(x))]
                nce_ids = [func() for func in nce_id]  # NCE_IDs for every single file
                self.unsorted_NCE_ID = nce_ids  # To store NCE for each file so it can be used for indexing later in different functions
                seen = set()
                nce_ids = [x for x in nce_ids if not (x in seen or seen.add(x))]  # Unique NCE_IDs
                print(f"All unique NCE and Isolations: {nce_ids}") if self.print_information else None
                # Ensures unique values while preserving order
                The_Lorax = []  # Protects the unique NCE values like he protected the trees
                Isolations = []  # Holds Isolations for any unique value
                # If a NCE value is the first instance of itself, it'll be just that number. If it repeats and the isolations are not the same, it'll add a new NCE value as NCE_ISOLATION
                for i in range(len(x)):
                    nce_isolation_id = str(x[i]) + "_" + str(y[i])  # Unique identifier as 'NCE_isolation'
                    if str(x[i]) not in NCE_values:
                        NCE_values.append(str(x[i]))  # Append as a string
                        The_Lorax.append(nce_isolation_id)
                        Isolations.append(str(y[i]))
                    elif (str(x[i]) + ' ' + str(y[i])) not in NCE_values and nce_isolation_id not in The_Lorax:
                        NCE_values.append(nce_isolation_id)
                        The_Lorax.append(nce_isolation_id)
                        Isolations.append(str(y[i]))
                # Removes duplicate entries from the final list
                NCE_values = list(dict.fromkeys(NCE_values))  # Keeps the order and removes duplicates

        # Optionally print the NCE values if requested
        print(f"NCE Values: {NCE_values}") if self.print_information else None
        self.nce_id = nce_ids  # NCE IDs
        self.NCE = NCE_values  # NCE Values where if it is the first instance of the NCE value it'll be that value, but if it's not it'll be NCE_ISOLATION. This also has no duplicates and is ordered in the order they're found
        self.The_Lorax = The_Lorax  # NCE IDs with preserved order and holds unique values only
        self.Isolations = Isolations  # Isolations for each Unique NCE value. This is ordered

    def only_include_CSV_files(self):
        """This function may not be necessary for all, however if your folder contains other files that aren't CSV files of the data you want to use for the data analysis, this function
        will sort through the files from self.all_files_path and know to only use the CSV files"""
        try:
            all_files = os.listdir(self.all_files_path)
            self.file_names = [file for file in all_files if file.endswith('.csv')]
            print("Included CSV files:", self.file_names) if self.print_information else None
        except FileNotFoundError:
            print(f"The directory '{self.all_files_path}' does not exist.") if self.errors else None
        except Exception as womp:
            print(f"An error occurred: {womp}") if self.warnings else None

    def define_NCE_Labels_and_Associated_Files(self):
        """Associates file names with NCE values (0 or 70 or 0(isolation)) based on the 'NCE' column in the info_on_files_data.
        Creates a numeric list from file names for proper indexing, ensuring that missing files are handled appropriately.
        ensure the file is named in the following format: filename_####.csv, as to properly map file names to file numbers,
        it'll drop the .csv and look at the digits after the _ ."""

        # Create a numeric list from file names by extracting numbers before the first underscore
        numeric_file_indices = []
        for filename in self.file_names:
            num_part = filename.split('_')[1].replace('.csv', '') if '_' in filename else ''
            numeric_file_indices.append(int(num_part))
        self.numeric_file_indices = numeric_file_indices

        # Initialize dictionaries for storing files based on plain NCE and formatted NCE
        initial_nce_files = {nce: [] for nce in self.NCE}
        Same_NCE_different_ISOLATION = []
        for ID in self.nce_id:
            if ID in self.NCE:
                Same_NCE_different_ISOLATION.append(ID)
        data = self.info_on_files_data
        x = data['NCE'].tolist()
        y = data['isolation'].tolist()

        # Populate NCE_ID with files, using formatted `nce_key` as keys
        for idx, row in self.info_on_files_data.iterrows():
            nce_value = str(row['NCE'])
            isolation = str(row['isolation'])
            nce_key = f"{nce_value}_{isolation}"
            file_number = int(row['file number'])

            # Match file_number with numeric_file_indices and retrieve the filename
            if file_number in numeric_file_indices:
                index_in_list = numeric_file_indices.index(file_number)
                filename = self.file_names[index_in_list]
                if nce_key in Same_NCE_different_ISOLATION:
                    initial_nce_files[nce_key].append(filename)
                else:
                    initial_nce_files[nce_value].append(filename)
            else:
                print(f"Notice: File number {file_number} not found in defined file names.") if self.notices else None

        self.NCE_Labels_and_Files = initial_nce_files

        # Create a dictionary for numeric file indices associated with each NCE key
        numeric_indices_per_nce = {}
        for nce_key, filenames in self.NCE_Labels_and_Files.items():
            numeric_indices = []
            for filename in filenames:
                num_part = filename.split('_')[1].replace('.csv', '') if '_' in filename else ''
                numeric_indices.append(int(num_part))
            numeric_indices_per_nce[nce_key] = numeric_indices

        self.numeric_file_indices_dict = numeric_indices_per_nce
        print(f"NCE Labels and Files: {self.NCE_Labels_and_Files}") if self.print_information else None
        return self.NCE_Labels_and_Files

    def drop_files(self, NCE_values_and_which_files_to_drop):
        """Removes specified files from the dataset based on NCE values and file numbers provided in a dictionary format.
        Args:
            NCE_values_and_which_files_to_drop (dict): Dictionary where keys are NCE values (e.g., '0', '70')
            and values are lists of file numbers to drop (e.g., {'0': [3, 6], '70': [5, 7]})."""
        # Create a mapping from filename to the extracted file index
        file_index_map = {file_name: int(file_name.split('_')[1].split('.')[0]) for file_name in self.file_names}

        for NCE_value, file_numbers_to_drop in NCE_values_and_which_files_to_drop.items():
            if NCE_value in self.NCE_Labels_and_Files:
                files_to_drop = [file_name for file_name, index in file_index_map.items() if index in file_numbers_to_drop]
                updated_files = [file for file in self.NCE_Labels_and_Files[NCE_value] if file not in files_to_drop]
                self.numeric_file_indices = [index for index in self.numeric_file_indices if index not in file_numbers_to_drop]
                self.NCE_Labels_and_Files[NCE_value] = updated_files

                if self.print_information:
                    print(f"Updated NCE {NCE_value}: {updated_files}")
            else:
                print(f"NCE {NCE_value} not found in the dictionary.") if self.warnings else None

    def get_file_paths(self, reverse_order=False):
        """Associates each file number with a string of a file path from info_on_files_data
        and returns a dictionary with file paths instead of file numbers. This function should be
        used after NCE_Labels_and_Files is defined and drop_files is used if necessary.

        Parameters:
        reverse_order (bool): If True, the dictionary value is reversed to correctly map the file_paths.
        It will correctly reverse the order back after mapping the corresponding file_paths.
        Use this if the folder has the files from the last file at the top to the first file at the bottom."""
        NCE_File_Paths = {}
        files = os.listdir(self.all_files_path)
        files.sort()  # Ensure files are sorted correctly

        for NCE_value, file_names in self.NCE_Labels_and_Files.items():
            file_paths = []
            for file_name in file_names:
                if file_name in files:
                    file_path = os.path.join(self.all_files_path, file_name)
                    file_paths.append(file_path)
                else:
                    print(f"Notice: File '{file_name}' not found in the directory.") if self.notices else None

            if reverse_order:
                file_paths.reverse()  # This restores the original order if needed

            NCE_File_Paths[NCE_value] = file_paths
            print(f"NCE File Paths for {NCE_value}: {file_paths}") if self.print_information else None

        self.NCE_File_Paths = NCE_File_Paths
        return NCE_File_Paths

    def merge__and_store_files_by_nce(self, nce_values=None, skip_rows=7):
        """Merges all the file_paths from Mass_Spectrometer_Data.NCE_File_Paths and returns a dictionary containing the Pandas data frames.
        Takes in NCE values in a list of the values you want to merge:

        # Option 1: Merges files for all NCEs (default behavior)
        merged_dataframes_all = the_method_over_all_files.merge_files_by_nce()

        # Option 2: Merges files only for specified NCE values (e.g., NCE '0' and '70')
        merged_dataframes_selected = the_method.merge_files_by_nce(nce_values=['0', '70'])

        And you can access it by doing this:
        Access merged DataFrame for NCE '0'
        nce_0_df = NCE_dataframes_selected.get('0')"""
        merged_dataframes = {}

        if not hasattr(self, 'NCE_File_Paths'):
            raise AttributeError("self does not have an attribute 'NCE_File_Paths'.") if self.errors else None

        if nce_values is None:
            nce_values = list(self.NCE_File_Paths.keys())

        for nce in nce_values:
            if nce in self.NCE_File_Paths:
                df_list = []
                for file_path in self.NCE_File_Paths[nce]:
                    try:
                        df = pd.read_csv(file_path, on_bad_lines='skip', skiprows=skip_rows)
                        df_list.append(df)
                    except pd.errors.ParserError as e:
                        print(f"ParserError: Problem reading {file_path}: {e}") if self.warnings else None
                    except Exception as e:
                        print(f"Error reading {file_path}: {e}") if self.warnings else None
                if df_list:
                    merged_dataframes[nce] = pd.concat(df_list, ignore_index=True)
                else:
                    print(f"No valid data for NCE {nce}") if self.warnings else None
            else:
                print(f"NCE value '{nce}' not found in the NCE_File_Paths dictionary.") if self.warnings else None

        self.NCE_dataframes = merged_dataframes
        return merged_dataframes
