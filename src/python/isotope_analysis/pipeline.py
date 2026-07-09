"""This is purely orchestration/dict-merging glue code around the other classes"""

import numpy as np

from . import state
from .mass_spec_data import Mass_Spectrometer_Data
from .error_propagation import Error_Propagation
from .graphs import Graphs


def all_calculations_and_graphs(molecule_name, name_of_atom_of_interest, date_data_was_taken, molecular_weight, AMU, all_file_path, info_on_file_path,
                                 row_starts_and_ends_parent, row_starts_and_ends_daughter, Bin_size, number_of_peaks, Normalization_Number,
                                 parent_daughter, ROWSKIP=7, Drop_files=None, Print_Info=True, Print_Warnings=True, Print_Notices=True, Print_Errors=True, professional=False,
                                 Normalization_Number_is_total_areas=False, rows_are_masses=False, graphs=False, save_PDF_Graph=False, save_PDF_Table=False, PDF_MAKER=None, Number_of_Decimals=None, Separate_files=False,
                                 file_path_for_separate_file=None, separate_file_info_on_file_path=None, Debug=False):
    """This function will produce all the intended graphs and perform the calculations needed to generate graphs listed within this function for data taken on the same day. You'd use this multiple times
    to add data to the global variables then use the functions above to generate graphs for all data all together. These are graphs where you can only make graphs under the same isolation. To use the charts that
    use data from multiple different isolations, you'll need to use

    What this function does and explanation of the arguments:

    It first makes the class instance, then takes a string of the molecule name, a string of the name of the atom you're looking at, a string of the data it was taken, an integer of the molecular weight
    of the molecule, and an integer value of the AMU of the atom, and defines it for later use. Next, it takes the string of the file_path that has all the data files and takes the string of the link to a public google sheets
    or other method of getting the information pertaining to the data and defines the characteristics of that data. Third, it takes the rows and ends of parents and daughters in lists of where each peak starts and ends
    (can be only one list if the data is evenly spaced by bins), an integer of the bin size of how many rows from row starts (inclusive) it'll inlcude to row end (inclusive), an integer of the number_of peaks you are including in the
    calcultions, a integer Normalization_Number to what you'd want to normalize the data to (1 if none), parent_daughter which is the two NCE values you're comparing as a string in a list ['parent','daughter'], and ROWSKIP which
    you can specify to be equal to the amount of rows you want to skip to get to the headers of the data in the event that your data doesn't start at the first row (in this function you don't need to account for any
    indexing differences, just put the rows as they appear in the documents based on indexes or masses depending). Additionally, if you have bad data files you can always specify which files to drop
     like {'0': [1,2], '70': [8,9]} which works even if theres only one specified. You can also specify Print_Info = False if you do not want to print out every final calculation, Print_Warnings = False if you want to not
    print out any warnings (Things that do not cause errors, but may cause an error in later functions), Print_Errors = False if you don't want to print any errors (Things that do not work and break functionality due to inputs),
    Print_Notices = False if you do not want to print any notices (things that do not cause errors or require a warning, but may be unintended on the Users part), professional = True if you want the graphs to not have color themes and all be the same color,
    Normalization_Number_is_total_areas = True (if you are normalizing by total area then this will automatically normalize it by this for you), rows_are_masses = True (if you want to bin by starting masses and ending masses instead
    of indexing by row), and graphs = False (if you don't want to make any graphs).

    The PDF_MAKER argument is for an instance of MakePDF class of which you'd save any data to such as tables or graphs. If save_PDF_Graph or Save_PDF_Table is True, the graphs and/or tables will save respectively. For this function,
    it will output the tables and graphs I set it to within the function, for more customization you'll have to change it yourself. You can also specify to what decimal place the values will be rounded with Number_of_Decimals = integer of how many decimals you want in the tables.

    In the event you need to use calculate_δα_indv, but your Parent and Daughter files are in seperate folders with different NCE_Values, then you'll have to Put Seperate_files = True and include the new file path to the folder containing said value in file_path_for_seperate_file. Additionally, recall
    that if it's not the first and only instance of that NCE you'll need to also specify the isolation in ['Parent', 'Daughter'] like '20_391(28)' for example.
    You'll also need to provide a google_sheet/spreadsheet in the same format as before for the seperate_file in seperate_file_info_on_file_path.

    Set Debug = True to run any information during intermediate steps meant for debugging.

    The names can be edited in multiple ways such as making the data just a name if you want to manually take it apart for specific analyses, you can just copy paste
    the workflow and fill out the data and change variable names. Everything is also specified in more detail in the actual functions themselves.

    NOTE: There is a case where you'd need to do something easy that I did not implement. When you run this function it stores total areas in the global variable for future use with a key being the dates the data was taken.
    If you took data for parent or daughter on multiple days and the parent or daughter (whichever you didn't use first) is less than or greater than the number of days you took data,
    you'd still need to have a parent and daughter for this function to work. It would have repeats of the same data. In this event, you should manually clear the global variable for parent or daughter masses and total areas
    after using this function n amount of times, then run the function after you have an m number of files. By this I mean if you have 3 days where you took the data relating to the parent data and only took data for this daughter on one day, you'll
    run the function with separate files arguments twice, then set All_Daughter_Areas = {} and Daughter_Masses = {}. After this you run the last function still with the separate files arguments and now you'll
    have the data for the parent NCE from the three seperate days and the daughter NCE from the one day for data analysis. To reiterate, in this case n = 3 and m = 1, you clear the respective global variable when m = n dates.

    For further clarification, I am defining Error, Warning, and Notice as follows:
    Error: Anything that breaks the function it is in, where it would not be able to continue so it stops where the Error occurs
    Warning: Anything that does not break the function it is in, but it will likely break a diffrent function that relies on anything that is changed or outputted within that function
    Notice: Anything that does not break the function it is in and will likely not break any of the other functions, but may be unintended by the user that this occured or it won't function as intended due to it
    All of these assume that there is an issue with what the User puts in and not the functions themselves.

    Quick Guide to Arguments:
    molecule_name: A string of the name of the molecule. It will be formatted if entered a format like Nd(NO_{3})_{4}^{-}.
    name_of_atom_of_interest: A string of the name of the Atom you're observing, so for Nd(NO_{3})_{4}^{-} this would be 'Nd'
    date_data_was_taken: A string of the date you took the data in the format of 'MM/DD/YYYY'. If the Parent and Daughter were taken on seperate dates, you can enter both dates or only one.
    molecular_weight: An integer containing the molecular weight of the molecule before it was broken up in a mass spectrometer.
    AMU: An integer of mass in AMU for your atom_of_interest
    all_file_path: A string containing the path that takes you to the collected data. Formatting is important for this function, so ensure the columns have the specified names used in other functions.
    info_on_file_path: A string of the path that takes you to a googlesheet (change function if spreadsheet) containing the information formatted in the specified format in other functions.
    row_starts_and_ends_parent. A list of lists containg the index of the rows where you want to start collecting data and stop collection data such as [[start_1,end_1], [start_2,end_2]] for the 'Parent' molecule
    row_starts_and_ends_daughterA list of lists containg the index of the rows where you want to start collecting data and stop collection data such as [[start_1,end_1], [start_2,end_2]] for the 'Daughter' molecule
    Bin_size: A integer of how many rows will be included in one bin, you should account for this in your row_starts_and_ends.
    number_of_peaks: An integer of how many peaks are you looking at in the spectrum that you are including in the data (account for this in Bin_size and row_starts_and_ends.)
    Normalization_Number: An inyeger of what you want to normalize the peaks by, if you choose to normalize by total areas then this can be an arbitary integer and you set the optional argument to True
    parent_daughter: An list containing two strings containing what the parent and daughter NCE values are. For example, ['0','70']
    ROWSKIP: An integer of how many rows you'd need to skip in all_data_files to get to the header of the data. If none then put zero.

    Optional Arguments:
    Drop_files = None # Put in a dictionary with keys and and lists of what file numbers to drop
    Print_Info = True # To print calculations and general information
    Print_Warnings = True # Print anything that may cause Errors later
    Print_Notices = True # Print any information someone may want to know that isn't causing any errors, but may be an unintended use of the function on the users part
    Print_Errors = True
    professional = True
    Normalization_Number_is_total_areas = False
    rows_are_masses = False
    graphs = False # Show the graphs you generate. Otherwise it'll skip it
    save_PDF_Graph = False
    save_PDF_Table = False
    PDF_MAKER = None
    Number_of_Decimals = None
    Seperate_files = False.  If you want an NCE value that's in a different folder, set as True
    file_path_for_seperate_file = None
    seperate_file_info_on_file_path = None
    Debug = False
    """

    MSD = Mass_Spectrometer_Data()  # Initialize a Mass_Spectrometer_Data class
    # Any attribute in MSD will be copied by The_Method class instance that is made in the Error_of_Data instance
    MSD.warnings = Print_Warnings  # True by default
    MSD.notices = Print_Notices  # True by default
    MSD.errors = Print_Errors  # True by default
    MSD.define_molecule_characteristics(molecule_name, name_of_atom_of_interest, date_data_was_taken, molecular_weight, AMU)  # Define molecule characteristics
    if Separate_files:
        MSD_2 = Mass_Spectrometer_Data()  # Make a second instance of this data class to get the files for the NCE in a Seperate_File
        MSD_2.warnings = Print_Warnings  # True by default
        MSD_2.notices = Print_Notices  # True by default
        MSD_2.errors = Print_Errors  # True by default
        # Ensuring That MSD and MSD_2 Functions without Breaking due to missing values for the other NCE.
        MSD.define_all_file_path(all_file_path)  # Gather the individual file paths from the path to the folder with all the files
        MSD_2.define_all_file_path(file_path_for_separate_file)  # Gather the individual file paths from the path to the folder with all the files
        MSD.define_NCE(info_on_file_path, GoogleSheets=True)  # Change this if it's not true for you
        MSD_2.define_NCE(separate_file_info_on_file_path, GoogleSheets=True)  # Change this if it's not true for you
        MSD.only_include_CSV_files()  # Only gathers CSV files
        MSD_2.only_include_CSV_files()  # Only gathers CSV files
        MSD.define_NCE_Labels_and_Associated_Files()  # Gathers the NCE labels from the google sheet/other spreadsheet and associates them to the file paths by indexes
        MSD_2.define_NCE_Labels_and_Associated_Files()  # Gathers the NCE labels from the google sheet/other spreadsheet and associates them to the file paths by indexes
        MSD.get_file_paths()
        MSD_2.get_file_paths()
        try:
            # Case where Parent key is in MSD and Daughter key is in MSD_2. It identifies this by checking the specific NCE_ID made the define_NCE function in the MSD class.
            # Separate files will assume that the specific NCE wouldn't be in both as that is its purpose. There is no print but this would be self.warnings if it were in both since this works if NCE_ID is in MSD.
            MSD_KEY = parent_daughter[0]  # Assumes MSD holds parent key
            MSD_2_KEY = parent_daughter[1]  # Assumes MSD_2 holds daughter key
            Error = False  # True if Parent key is not in MS
            # Compare Parent key to NCE_IDs and NCE. If it were '0' that means it's the first instance of '0' otherwise it'd be '0_391(28)', for example, which is 'NCE_ISOLATION'. However both files could have a zero, so isolation may matter
            if MSD_KEY in MSD.NCE:  # Check if the key is in the NCE values
                if MSD_KEY in MSD_2.NCE:  # Check if the NCE value also happens to be in MSD_2. If this is the case we have to check the isolation with the NCE value stored in the IDs
                    # At this point in the code we can see that MSD 2 and MSD 1 have an NCE value like '0' with no specific isolation specified. So now we need to see if '0' for Parent key exists in MSD
                    Index = MSD.NCE.index(MSD_KEY)  # Get the appropriate index
                    if MSD_KEY + "_" + MSD.Isolations[Index] == MSD.nce_id[Index]:  # Check if the isolation of the Key matches the NCE_id made using MSD instance data.
                        pass  # The Parent key is in MSD.NCE values and the isolation matches the NCE_ISOLATION IDs created for MSD. Therefore we can proceed with this try block
                    else:
                        Error = True  # The Parent key is in MSD.NCE values, but after we compared the NCE_ID we find that they are not the same Isolation. So this means the parent key is in MSD_2
                else:
                    pass  # The Parent key is in MSD.NCE, but not MSD_2.NCE, so the Parent key is in MSD and we can continue with the try block
            else:
                Error = True  # The Parent key is not in MSD.NCE values in any form
            if Error:
                raise ValueError("Parent key is not in MSD.NCE values")
            print('Parent NCE value is in the first folder, the daughter NCE value is in the separate folder') if Print_Notices or Debug else None
            if Drop_files is not None:
                if MSD_KEY in Drop_files:
                    MSD.drop_files({MSD_KEY: Drop_files[MSD_KEY]})
                if MSD_2_KEY in Drop_files:
                    MSD_2.drop_files({MSD_2_KEY: Drop_files[MSD_2_KEY]})
            if Debug:
                print(f"NCE for MSD before combining it: {MSD.NCE}")
                print(f"NCE for MSD_2 before combining it: {MSD_2.NCE}")
                print(f"numeric_file_indices for MSD beofre combining: {MSD.numeric_file_indices}")
                print(f"numeric_file_indices for MSD_2 beofre combining: {MSD_2.numeric_file_indices}")
                print(f"File names in MSD.file_names before combining: {MSD.file_names}")
                print(f"File names in MSD_2.file_names before combining: {MSD_2.file_names}")
                print(f"Dictionary with NCE keys and corresponding files in MSD before combining: {MSD.NCE_Labels_and_Files}")
                print(f"Dictionary with NCE keys and corresponding files in MSD_2 before combining: {MSD_2.NCE_Labels_and_Files}")
                print(f"The Dictionary containing the NCE keys and their corresponding file paths in MSD before combination: {MSD.NCE_File_Paths}")
                print(f"The Dictionary containing the NCE keys and their corresponding file paths in MSD_2 before combination: {MSD_2.NCE_File_Paths}")
            # Check and remove duplicates from MSD_2 before merging. We have confirmed the parent NCE is in MSD, but if there were a key of the same NCE and isolation in MSD 2, it would cause errors. To deal with this we do a check to remove the shared key from the daughter dictionary

            # Removes any similar keys from other since it would lead to duplications. However we also want to remove MSD_2 from MSD for that reason as if we didn't we'd get the wrong values when we combine them
            # Remove entire key from MSD_2.numeric_file_indices_dict if key exists in MSD.numeric_file_indices_dict, but keep the specified key MSD_2_KEY
            for key in list(MSD_2.numeric_file_indices_dict.keys()):
                if key in MSD.numeric_file_indices_dict and key != MSD_2_KEY:
                    print(f"Removing duplicate key '{key}' from MSD_2.numeric_file_indices_dict") if Debug else None
                    del MSD_2.numeric_file_indices_dict[key]
                    if MSD_2_KEY in MSD.numeric_file_indices_dict:
                        print(f"Also removing {MSD_2_KEY} from MSD.numeric_file_indices_dict") if Debug else None
                        del MSD.numeric_file_indices_dict[MSD_2_KEY]

            # Remove entire key from MSD_2.NCE_Labels_and_Files if key exists in MSD.NCE_Labels_and_Files, but keep the specified key MSD_2_KEY
            for key in list(MSD_2.NCE_Labels_and_Files.keys()):
                if key in MSD.NCE_Labels_and_Files and key != MSD_2_KEY:
                    print(f"Removing duplicate key '{key}' from MSD_2.NCE_Labels_and_Files") if Debug else None
                    del MSD_2.NCE_Labels_and_Files[key]
                    if MSD_2_KEY in MSD.NCE_Labels_and_Files:
                        print(f"Also removing {MSD_2_KEY} from MSD.NCE_Labels_and_Files") if Debug else None
                        del MSD.NCE_Labels_and_Files[MSD_2_KEY]

            # Remove entire key from MSD_2.NCE_File_Paths if key exists in MSD.NCE_File_Paths, but keep the specified key MSD_2_KEY
            for key in list(MSD_2.NCE_File_Paths.keys()):
                if key in MSD.NCE_File_Paths and key != MSD_2_KEY:
                    print(f"Removing duplicate key '{key}' from MSD_2.NCE_File_Paths") if Debug else None
                    del MSD_2.NCE_File_Paths[key]
                    if MSD_2_KEY in MSD.NCE_File_Paths:
                        print(f"Also removing {MSD_2_KEY} from MSD.NCE_File_Paths") if Debug else None
                        del MSD.NCE_File_Paths[MSD_2_KEY]

            # Merging MSD.NCE_Labels_and_Files
            merged_dict = {}
            for key in set(MSD.NCE_Labels_and_Files) | set(MSD_2.NCE_Labels_and_Files):  # Union of keys. Which is fine because they should be the exact same isolation
                if key in MSD.NCE_Labels_and_Files and key in MSD_2.NCE_Labels_and_Files:
                    merged_dict[key] = MSD.NCE_Labels_and_Files[key] + MSD_2.NCE_Labels_and_Files[key]
                else:
                    merged_dict[key] = MSD.NCE_Labels_and_Files.get(key, MSD_2.NCE_Labels_and_Files.get(key))
            MSD.NCE_Labels_and_Files = merged_dict
            print(f"MSD.NCE_Labels_and_Files: {MSD.NCE_Labels_and_Files}") if Debug else None

            # Merging MSD.NCE_File_Paths
            merged_dict = {}  # Clear merged_dict from previous function
            for key in set(MSD.NCE_File_Paths) | set(MSD_2.NCE_File_Paths):  # Union of keys. As they have to be the same isolation so all files are valid to be put under the same key.
                if key in MSD.NCE_File_Paths and key in MSD_2.NCE_File_Paths:
                    merged_dict[key] = MSD.NCE_File_Paths[key] + MSD_2.NCE_File_Paths[key]
                else:
                    merged_dict[key] = MSD.NCE_File_Paths.get(key, MSD_2.NCE_File_Paths.get(key))
            MSD.NCE_File_Paths = merged_dict
            merged_numeric_indices = {}
            for key in set(MSD.numeric_file_indices_dict) | set(MSD_2.numeric_file_indices_dict):  # Union of keys
                if key in MSD.numeric_file_indices and key in MSD_2.numeric_file_indices:
                    merged_numeric_indices[key] = MSD.numeric_file_indices_dict[key] + MSD_2.numeric_file_indices_dict[key]
                else:
                    merged_numeric_indices[key] = MSD.numeric_file_indices_dict.get(key, MSD_2.numeric_file_indices_dict.get(key))
            MSD.numeric_file_indices_dict = merged_numeric_indices
            MSD.numeric_file_indices = MSD.numeric_file_indices_dict[MSD_KEY] + MSD.numeric_file_indices_dict[MSD_2_KEY]
            print(f"new numeric file indices: {MSD.numeric_file_indices}") if Debug else None
            print(f"Combined NCE File Path dictionary: {MSD.NCE_File_Paths}") if Debug else None
        except:  # noqa: E722
            print('Parent NCE value is in in the separate folder') if Print_Notices or Debug else None
            # Case where Daughter key is in MSD and Parent key is in MSD_2
            MSD_KEY = parent_daughter[1]
            MSD_2_KEY = parent_daughter[0]
            if Drop_files is not None:
                if MSD_KEY in Drop_files:
                    MSD.drop_files({MSD_KEY: Drop_files[MSD_KEY]})
                if MSD_2_KEY in Drop_files:
                    MSD_2.drop_files({MSD_2_KEY: Drop_files[MSD_2_KEY]})
            if Debug:
                print(f"NCE for MSD before combining it: {MSD.NCE}")
                print(f"NCE for MSD_2 before combining it: {MSD_2.NCE}")
                print(f"numeric_file_indices for MSD before combining: {MSD.numeric_file_indices}")
                print(f"numeric_file_indices for MSD_2 beofre combining: {MSD_2.numeric_file_indices}")
                print(f"File names in MSD.file_names before combining: {MSD.file_names}")
                print(f"File names in MSD_2.file_names before combining: {MSD_2.file_names}")
                print(f"Dictionary with NCE keys and corresponding files in MSD before combining: {MSD.NCE_Labels_and_Files}")
                print(f"Dictionary with NCE keys and corresponding files in MSD_2 before combining: {MSD_2.NCE_Labels_and_Files}")
                print(f"The Dictionary containing the NCE keys and their corresponding file paths in MSD before combination: {MSD.NCE_File_Paths}")
                print(f"The Dictionary containing the NCE keys and their corresponding file paths in MSD_2 before combination: {MSD_2.NCE_File_Paths}")
            # Check and remove duplicates from MSD_2 before merging. Although we confirmed parent is in MSD_2, there is still the possibility of the same NCE and isolation in MSD which would cause errors. To deal with this we do a check to that key from the daughter dictionary

            # Removes any similar keys from other since it would lead to duplications. However we also want to remove MSD from MSD_2 for that reason as if we didn't we'd get the wrong values when we combine them
            for key in list(MSD.numeric_file_indices_dict.keys()):
                if key in MSD_2.numeric_file_indices_dict and key != MSD_KEY:
                    print(f"Removing duplicate key '{key}' from MSD.numeric_file_indices_dict") if Debug else None
                    del MSD.numeric_file_indices_dict[key]
                    if MSD_KEY in MSD_2.numeric_file_indices_dict:
                        print(f"Also removing {MSD_KEY} from MSD_2.numeric_file_indices_dict") if Debug else None
                        del MSD_2.numeric_file_indices_dict[MSD_KEY]

            for key in list(MSD.NCE_Labels_and_Files.keys()):
                if key in MSD_2.NCE_Labels_and_Files and key != MSD_KEY:
                    print(f"Removing duplicate key '{key}' from MSD.NCE_Labels_and_Files") if Debug else None
                    del MSD.NCE_Labels_and_Files[key]
                    if MSD_KEY in MSD_2.NCE_Labels_and_Files:
                        print(f"Also removing {MSD_KEY} from MSD_2.NCE_Labels_and_Files") if Debug else None
                        del MSD_2.NCE_Labels_and_Files[MSD_KEY]

            for key in list(MSD.NCE_File_Paths.keys()):
                if key in MSD_2.NCE_File_Paths and key != MSD_KEY:
                    print(f"Removing duplicate key '{key}' from MSD.NCE_File_Paths") if Debug else None
                    del MSD.NCE_File_Paths[key]
                    if MSD_KEY in MSD_2.NCE_File_Paths:
                        print(f"Also removing {MSD_KEY} from MSD_2.NCE_File_Paths") if Debug else None
                        del MSD_2.NCE_File_Paths[MSD_KEY]

            # Merging MSD_2.NCE_File_Paths
            merged_dict = {}
            for key in set(MSD_2.NCE_File_Paths) | set(MSD.NCE_File_Paths):  # Union of keys
                if key in MSD_2.NCE_File_Paths and key in MSD.NCE_File_Paths:
                    merged_dict[key] = MSD_2.NCE_File_Paths[key] + MSD.NCE_File_Paths[key]
                else:
                    merged_dict[key] = MSD_2.NCE_File_Paths.get(key, MSD.NCE_File_Paths.get(key))
            MSD.NCE_File_Paths = merged_dict

            # Merging Numeric_File_Indices
            merged_numeric_indices = {}
            for key in set(MSD_2.numeric_file_indices_dict) | set(MSD.numeric_file_indices_dict):  # Union of keys
                if key in MSD_2.numeric_file_indices and key in MSD.numeric_file_indices:
                    merged_numeric_indices[key] = MSD_2.numeric_file_indices_dict[key] + MSD.numeric_file_indices_dict[key]
                else:
                    merged_numeric_indices[key] = MSD_2.numeric_file_indices_dict.get(key, MSD.numeric_file_indices_dict.get(key))
            MSD.numeric_file_indices_dict = merged_numeric_indices
            MSD.numeric_file_indices = MSD.numeric_file_indices_dict[MSD_2_KEY] + MSD.numeric_file_indices_dict[MSD_KEY]
            print(f"new numeric file indices: {MSD.numeric_file_indices}") if Debug else None
            print(f"Combined NCe file path dictionary: {MSD.NCE_File_Paths}") if Debug else None
    else:
        # All Files for specified NCEs are in the same folder
        print('Parent and Daughter NCE values are in the same folder') if Print_Notices or Debug else None
        MSD.define_all_file_path(all_file_path)  # Gather the individual file paths from the path to the folder with all the files
        MSD.define_NCE(info_on_file_path, GoogleSheets=True)  # Change this if it's not true for you
        MSD.only_include_CSV_files()  # Only gathers CSV files
        MSD.define_NCE_Labels_and_Associated_Files()  # Gathers the NCE labels from the google sheet/other spreadsheet and associates them to the file paths by indexes
        if Drop_files is not None:
            MSD.drop_files(Drop_files)  # Drops files for any of the NCEs if specifed
        MSD.get_file_paths()  # Gathers all the file paths for the defined NCEs and makes a dictionary where each NCE has defined file paths
    MSD.merge__and_store_files_by_nce(MSD.NCE)  # Makes a big data frame by combining all data with same values together, but its use here is to ensure functionality as it has multiple uses and needs to be defined for some of the functions used here to work
    Error_of_Data = Error_Propagation(MSD)  # Initalize a Error class taking in the data from the MSD instance
    Error_of_Data.print_information = Print_Info  # True by default, will print out all calculations
    if rows_are_masses:  # If you chose to index by mass instead of the row indexes directly, otherwise it'll continue with indexing by rows
        Error_of_Data.calculate_δα_indv(row_starts_and_ends_parent, row_starts_and_ends_daughter, Bin_size, number_of_peaks, Normalization_Number, parent_daughter, ROWSKIP=7, NBTA=Normalization_Number_is_total_areas, rows_are_masses=True)
    else:
        Error_of_Data.calculate_δα_indv(row_starts_and_ends_parent, row_starts_and_ends_daughter, Bin_size, number_of_peaks, Normalization_Number, parent_daughter, ROWSKIP=7, NBTA=Normalization_Number_is_total_areas)
    if save_PDF_Table:
        custom_labels = []
        for isotope in state.All_Isotopes:
            isotope_label = f"$\\mathrm{{^{{{isotope}}}{state.Atom_Names}}}$"
            custom_labels.append(isotope_label)  # Custom labels for Isotopes
        if Separate_files:
            # Make title distinction because the default title assumes all the data was in the same folder taken on the same day
            Current_Table_Title = f"Tabulated Computational Outcomes for {MSD.molecule_name} from Data taken on Different Dates"
        else:
            Current_Table_Title = f"Tabulated Computational Outcomes for {MSD.molecule_name} Based on {MSD.date_data_was_taken} Data"
        PDF_MAKER.add_table(PDF_MAKER.define_table(custom_labels, np.array(['Alpha', 'Enrichment', 'Error Bar', 'Parent_Ratio', 'Daughter_Ratio']), [Error_of_Data.alphas, Error_of_Data.enrichments, Error_of_Data.δα, Error_of_Data.Parent_ratios, Error_of_Data.Daughter_ratios], pivot=True, number_of_decimals=Number_of_Decimals), title=Current_Table_Title,
                          Additional_Information=f"The NCEs are Parent: {parent_daughter[0]} {'|'} Daughter: {parent_daughter[1]}")  # Table totals
    if graphs:
        # To make the graphs if graphs is True
        Graphing_Molecule_Data = Graphs(MSD, Error_of_Data)  # Intalize a Graphing class that takes in the MSD and Error_of_Data instances
        Graphing_Molecule_Data.generate_all_IVM(387.5, 400, 340, 355, ['0', '70'])  # You can change the bounds of the data with xlim left and right as well as change which NCE they are for within this function
    elif save_PDF_Graph:
        # Save the Graphs to the PDF if save_PDF_Graph is True. Same Descriptions as in Graphs
        Graphing_Molecule_Data = Graphs(MSD, Error_of_Data)  # Intalize a Graphing class that takes in the MSD and Error_of_Data instances
        Graphing_Molecule_Data.professional = professional  # False by default, will make all graphs the same color
        PDF_MAKER.add_chart(Graphing_Molecule_Data.αVIs, save_pdf=save_PDF_Graph)
        PDF_MAKER.add_chart(Graphing_Molecule_Data.εVIs, save_pdf=save_PDF_Graph)
        Graphing_Molecule_Data.RIVM(PDF_Maker=PDF_MAKER)
