"""MakePDF and Add_Combined_Data_Table_To_PDF, ported from Foundation.ipynb. Table
sizes here are bounded by number_of_isotopes, and PDF writing is matplotlib/IO-bound, so nothing here scales with data size."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from . import state


class MakePDF:
    """Example Use of MakePDF:

    Use the MakePDF class with a 'with' statement:
    with MakePDF("example_with_context.pdf", "path/to/save") as pdf_creator:
    pdf_creator.add_chart(generate_line_chart, data)
    Add a table to the PDF:
    pdf_creator.add_table(table_data, column_names=column_names, additional_rows=additional_rows)

    Alternatively, you can use each function individually if the table requires some changes that you don't want to manually type out
    """

    def __init__(self, pdf_filename, file_path=None):
        """Initializes the PDF generator.

        Args:
            pdf_filename (str): The name of the PDF file where content will be saved.
        """
        if not pdf_filename.endswith('.pdf'):
            pdf_filename += '.pdf'
        self.pdf_path = pdf_filename if not file_path else file_path + '/' + pdf_filename
        self.pdf_pages = PdfPages(self.pdf_path)
        self.table_data = None  # Holds the table itself in it's list of list form.
        self.debug = False  # Set to True if you want to have debug statements on.
        self.warnings = False
        self.Errors = False
        self.notice = False

    def pivot_table_data(self, *args):
        """Pivots columns of lists of data into rows for the table.

        Args:
            *args: Multiple lists of data to be arranged vertically into table rows.

        Returns:
            list: A list of rows (lists), each containing values from the input lists.
        """
        table_data = list(zip(*args))  # This effectively transposes the data
        return table_data

    def handle_floats_and_transpose(self, data_sets, Number_of_decimals=None):
        float_counter = 1
        placeholder_map = {}  # to store float->placeholder mappings
        transformed_data_sets = []

        # Step 1: Replace float values with unique codes
        for data in data_sets:
            print(f"Processing data set: {data}") if self.debug else None
            transformed_data = []

            # Ensure all rows have the same length
            row_lengths = [len(row) for row in data]
            if len(set(row_lengths)) > 1:
                print("Warning: Rows have unequal lengths in data set. Skipping transposition for inconsistent data.") if self.warnings else None
                continue

            for row in data:
                print(f"Processing row: {row}") if self.debug else None
                transformed_row = []
                for value in row:
                    print(f"Processing value: {value} (type: {type(value)})") if self.debug else None
                    if isinstance(value, np.float64):
                        # Assign a unique code to the floating point value
                        placeholder = f'FLOAT{float_counter}'
                        transformed_row.append(placeholder)
                        placeholder_map[placeholder] = value
                        float_counter += 1
                    else:
                        transformed_row.append(value)
                transformed_data.append(transformed_row)
            transformed_data_sets.append(transformed_data)

        # Step 2: Transpose the data (pivot)
        print("Transposing data sets...") if self.debug else None
        pivoted_data_sets = [list(zip(*data)) for data in transformed_data_sets]

        # Step 3: Restore the floats from the placeholders
        restored_data_sets = []
        for data in pivoted_data_sets:
            restored_data = []
            for row in data:
                restored_row = []
                for value in row:
                    if isinstance(value, str) and value.startswith('FLOAT'):
                        original_value = placeholder_map.get(value, value)
                        if Number_of_decimals is not None:
                            original_value = round(original_value, Number_of_decimals)
                        restored_row.append(original_value)
                    else:
                        restored_row.append(value)
                restored_data.append(restored_row)
            restored_data_sets.append(restored_data)

        return restored_data_sets

    def define_table(self, column_names, row_names, *data_sets, pivot=False, number_of_decimals=None):
        """Defines the table with the specified row names, column headers, and multiple data sets.
        Optionally pivots each data set (switches rows and columns).

        Args:
            column_names (list): List of column names (e.g., ['File 0', 'File 1', 'File 2']).
            row_names (list): List of row names (e.g., ['Alpha', 'Enrichment', 'Error']).
            *data_sets (lists): Multiple lists of data (e.g., [enrichment_data], [alpha_data], [ratio_data]).
            pivot (bool): If True, pivots the data sets (rows become columns and columns become rows).
            number_of_decimals (int): Number of decimal places to include when displaying data. It will round the last digit to the number of decimals specified.
        Returns:
            list: A 2D list representing the table with row names, column headers, and data.
        """
        # Step 1: Handle floaters and transpose the data
        if pivot:
            data_sets = self.handle_floats_and_transpose(data_sets, Number_of_decimals=number_of_decimals)
        else:
            if number_of_decimals is not None:
                rounded_data_sets = []
                for data in data_sets:
                    rounded_data = [
                        [round(value, number_of_decimals) if isinstance(value, (float, np.float64)) else value for value in row]
                        for row in data
                    ]
                    rounded_data_sets.append(rounded_data)
                data_sets = rounded_data_sets

        # Step 2: Create the final table with an empty cell in the first column and column headers
        table_data = []
        table_data.append([''] + row_names.tolist())

        for i, file_name in enumerate(column_names):
            row = [file_name]
            for data in data_sets:
                row.extend(data[i])
            table_data.append(row)

        self.table_data = table_data
        return table_data

    def truncate_text(cell_text, max_char):
        """Truncates text to fit within a specified character limit per cell."""
        if isinstance(cell_text, str) and len(cell_text) > max_char:
            return cell_text[:max_char] + "..."
        return cell_text

    def add_table(self, table_data, title='I owe you a title', Additional_Information=None):
        """Adds a table to the PDF.
        Args:
            table_data (list): A 2D list representing the table data. You can pass in define_table() here
            title (str): A title for the table.
            Additional_Information (str): Additional information to be displayed below the table. Such as Isolations and NCE values.
        """
        ax = plt.gca()
        ax.axis('tight')
        ax.axis('off')

        table = plt.table(
            cellText=table_data,
            loc='center',
            colWidths=[0.2] * (len(table_data[0]) if table_data else 1),
            cellLoc='center'
        )

        table.auto_set_font_size(False)
        table.set_fontsize(12)
        table.scale(1.5, 1.5)

        if title:
            plt.title(title, fontsize=14, weight='bold', pad=20)

        table.auto_set_font_size(False)
        table.set_fontsize(10)

        if Additional_Information is not None:
            ax.text(0.5, -0.5, Additional_Information,
                     fontsize=10, ha='center', transform=ax.transAxes)

        self.pdf_pages.savefig(bbox_inches='tight', dpi=300)
        plt.close()

    def add_chart(self, chart_function, **additional_arguments):
        """Adds a chart to the PDF using an existing chart function.

        Args:
            chart_function (function): A pre-existing function that generates a chart (e.g., `overlaid_line_plot_AVI,`).
            **additional_arguments: Additional arguments to pass to the chart function.
        """
        plt.figure(figsize=(10, 5))
        chart_function()
        chart_function(**additional_arguments)  # Pass kwargs to the chart function
        self.pdf_pages.savefig(dpi=300)
        plt.close()

    def save_pdf(self):
        """Closes the PDF file and saves all pages."""
        self.pdf_pages.close()

    def __enter__(self):
        """Start the context and open the PDF for writing."""
        self.pdf_pages = PdfPages(self.pdf_path)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """End the context and close the PDF."""
        if self.pdf_pages:
            self.pdf_pages.close()


def Add_Combined_Data_Table_To_PDF(Normalization_Number=1, Normalization_Number_is_total_areas=True, save_Table=False, PDF_MAKER=None, Pivot=False, Number_of_Decimals=None, additional_Information=None):
    """Generates a table with columns I defined for combined data for final graphs after using all_calculations_and_graphs multiple times or calculate_δα_indv multiple times.
    The table will have columns Total Areas Parent, Total Areas Daughter, α ε, δα, and preference. The rows will correspond to the Isotopes. Make an instance/put an instance of
    MakePDF in the PDF_Maker argument to put where it'll store the table. Put pivot as true if the data in the lists will be put in the table vertically as a column."""
    custom_labels = []
    for isotope in state.All_Isotopes:
        isotope_label = f"$\\mathrm{{^{{{isotope}}}{state.Atom_Names}}}$"
        custom_labels.append(isotope_label)

    # Sums values across all the tuples resulting in one tuple of all the summed areas.
    summed_tuple_Parent = tuple(sum(values) for values in zip(*state.All_Parent_Areas.values()))
    summed_tuple_Daughter = tuple(sum(values) for values in zip(*state.All_Daughter_Areas.values()))

    Combined_Error_Parent = np.sqrt(np.array(summed_tuple_Parent))
    Combined_Error_Daughter = np.sqrt(np.array(summed_tuple_Daughter))

    # Turns the tuples into an array to get the ratios which normalizes it
    Parent_Ratios = np.array(summed_tuple_Parent) / np.array(summed_tuple_Parent).sum() if Normalization_Number_is_total_areas else np.array(summed_tuple_Parent) / Normalization_Number
    Daughter_Ratios = np.array(summed_tuple_Daughter) / np.array(summed_tuple_Daughter).sum() if Normalization_Number_is_total_areas else np.array(summed_tuple_Daughter) / Normalization_Number

    # Calculating the errors and alphas of the newly combined data
    Alphas_of_Combined_Data = tuple(Daughter_Ratios / Parent_Ratios)
    Enrichments_of_Combined_Data = tuple((Daughter_Ratios / Parent_Ratios) - 1)
    Preferences = ['Parent' if ε < 0 else 'Neither' if ε == 0 else 'Daughter' for ε in Enrichments_of_Combined_Data]
    Errors_of_Combined_Data = tuple((Daughter_Ratios / Parent_Ratios) * np.sqrt(
        (np.sqrt(np.array(summed_tuple_Daughter)) / np.array(summed_tuple_Daughter)) ** 2 +
        (np.sqrt(np.array(summed_tuple_Parent)) / np.array(summed_tuple_Parent)) ** 2))

    # Round off Data to the decimal place specified in Number_of_Decimals if Number_of_Decimals is an integer
    if Number_of_Decimals is not None:
        summed_tuple_Parent = tuple(round(value, Number_of_Decimals) for value in summed_tuple_Parent)
        summed_tuple_Daughter = tuple(round(value, Number_of_Decimals) for value in summed_tuple_Daughter)
        Combined_Error_Parent = np.round(Combined_Error_Parent, Number_of_Decimals)
        Combined_Error_Daughter = np.round(Combined_Error_Daughter, Number_of_Decimals)
        Parent_Ratios = np.round(Parent_Ratios, Number_of_Decimals)
        Daughter_Ratios = np.round(Daughter_Ratios, Number_of_Decimals)
        Alphas_of_Combined_Data = tuple(round(value, Number_of_Decimals) for value in Alphas_of_Combined_Data)
        Enrichments_of_Combined_Data = tuple(round(value, Number_of_Decimals) for value in Enrichments_of_Combined_Data)
        Errors_of_Combined_Data = tuple(round(value, Number_of_Decimals) for value in Errors_of_Combined_Data)

    # Make the Table
    Table_of_Ultimate_Power = PDF_MAKER.define_table(custom_labels, np.array(['α', 'ε', 'Parent Error', 'Daughter Error', 'Error Bar', 'Preference']), np.array([Alphas_of_Combined_Data, Enrichments_of_Combined_Data, Combined_Error_Parent, Combined_Error_Daughter, Errors_of_Combined_Data, Preferences]), pivot=Pivot)
    if save_Table:
        Current_Table_Title = "Integrated Overview of All Data linked to the Same Isolation"
        PDF_MAKER.add_table(Table_of_Ultimate_Power, title=Current_Table_Title, Additional_Information=additional_Information)
    if save_Table:
        return Table_of_Ultimate_Power
