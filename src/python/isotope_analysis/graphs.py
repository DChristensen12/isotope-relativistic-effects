"""Graphs class and the standalone cross-run plotting functions, ported from
Foundation.ipynb. All matplotlib calls on arrays sized by number_of_peaks /
number_of_isotopes, so nothing here scales with data size. 
Global dict reads (All_alphas, All_Isotopes, etc.) go through `state`, same as
in error_propagation.py."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from . import state


class Graphs:
    """Used to generate all types of Graphs for the data between two specific values after all the data has been produced in the other Classes. There are graphs that are in the global function section.
    The save_pdf in each graph is used to not do plt.show() so that it can be saved as a figure instead in the Make_PDF class (True to make_pdfs, False to plot the graph)."""

    def __init__(self, Mass_Spectrometer_Data_Instance, Error_instance):
        self.MSDI = Mass_Spectrometer_Data_Instance  # Specify for which graphs you want to create.
        self.error = Error_instance  # Specify Errors for which graphs you'd want to create
        self.professional = False

    def IVM(self, file_path, xlim_left, xlim_right, plot_title, save_pdf=False):
        """Generates a Intensity Vs Mass graph from the data
        Enter how many rows you need to skip to access the data into row_skip as an integer
        Enter where you want the graphs to start and end to zoom in as integers, into xlim_left and xlim_right respectively
        Enter the plot_title as a string of what you want to name the graph into plot_title"""
        data = pd.read_csv(file_path, skiprows=7)
        Column_labels = ['Mass', 'Intensity']
        data.columns = Column_labels
        x_axis = 'Mass'
        y_axis = 'Intensity'
        plt.plot(data[x_axis], data[y_axis], label=f'{y_axis} vs {x_axis}')
        plt.xlabel(x_axis)
        plt.ylabel(y_axis)
        plt.title(plot_title)
        plt.xlim(xlim_left, xlim_right)
        if not save_pdf:
            plt.show()

    def generate_all_IVM(self, xlim_left_P, xlim_right_P, xlim_left_D, xlim_right_D, Parent_Daughter, save_pdf=False):
        """Generates all Intensity versus Masses graphs for every file for the parent and daughter species. I assume you only have two for this"""
        Parent_key = Parent_Daughter[0]  # To get the file_paths for the specifc NCE values in NCE_File_Paths
        Daughter_key = Parent_Daughter[1]
        P = 0  # Initalize counters to keep number the files in order
        D = 0
        for i in self.MSDI.NCE_File_Paths[Parent_key]:
            plot_title_P = "Intensity Vs Mass" + " " + "For Parent File:" + " " + str(P) + " " + "on" + " " + self.MSDI.date_data_was_taken
            self.IVM(i, xlim_left_P, xlim_right_P, plot_title_P) if not save_pdf else self.IVM(i, xlim_left_P, xlim_right_P, plot_title_P, save_pdf=True)
            P += 1
        for j in self.MSDI.NCE_File_Paths[Daughter_key]:
            plot_title_D = "Intensity Vs Mass" + " " + "for Daughter File:" + " " + str(D) + " " + "on" + " " + self.MSDI.date_data_was_taken
            self.IVM(j, xlim_left_D, xlim_right_D, plot_title_D) if not save_pdf else self.IVM(j, xlim_left_D, xlim_right_D, plot_title_D, save_pdf=True)
            D += 1

    def Big_IVM(self, Parent_Daughter, save_pdf=False):
        """IVM showing the total combined data frames from all the files for each NCE Parent_Daughter is list of strings for NCE values you are representing as parent and daughter such as ['0','70']"""
        Parent_key = Parent_Daughter[0]
        Daughter_key = Parent_Daughter[1]
        Parent = self.MSDI.NCE_dataframes[Parent_key]
        Daughter = self.MSDI.NCE_dataframes[Daughter_key]
        Column_labels = ['Mass', 'Intensity']
        Parent.columns = Column_labels
        Daughter.columns = Column_labels
        x_axis = 'Mass'
        y_axis = 'Intensity'
        plt.plot(Parent[x_axis], Parent[y_axis], label=f'{y_axis} vs {x_axis}')
        plt.xlabel(x_axis)
        plt.ylabel(y_axis)
        plt.title("Fragmentation Analysis: Intensity Patterns Across Parent Molecule Masses" + self.MSDI.date_data_was_taken)
        plt.show()
        plt.plot(Daughter[x_axis], Daughter[y_axis], label=f'{y_axis} vs {x_axis}')
        plt.xlabel(x_axis)
        plt.ylabel(y_axis)
        plot_title_D = "Fragmentation Analysis: Intensity Patterns Across Daughter Molecule Masses" + self.MSDI.date_data_was_taken
        plt.title(plot_title_D)
        if not save_pdf:
            plt.show()

    def εVIs(self, save_pdf=False):
        """Generates an Basic Enrichment vs Isotope Graph with error bars from data calculated in a method. Make the plot title a string of the title you want."""
        Column_labels = ['Isotope', 'ε']
        x_axis = self.error.Isotopes
        y_axis = self.error.enrichments
        yerr = self.error.δα
        custom_labels = []
        for isotope in self.error.Isotopes:
            isotope_label = f"$\\mathrm{{^{{{isotope}}}{self.MSDI.atom_name}}}$"
            custom_labels.append(isotope_label)  # Custom labels
        plt.plot(x_axis, y_axis, color='#1a2b6d') if self.professional else plt.plot(x_axis, y_axis, color=(0.6, 0.8, 1.0))
        plt.xticks(x_axis, custom_labels)
        plt.xlabel('Isotope')
        plt.ylabel('ε')
        plt.title(("Exploring the Relationship Between Enrichment Levels and Isotopic Composition" + " " + self.MSDI.date_data_was_taken))
        plt.errorbar(x_axis, y_axis, yerr=yerr, fmt='s', color='#8b0000', ecolor='#888888', capsize=15) if self.professional else plt.errorbar(x_axis, y_axis, yerr=yerr, fmt='s', color=(1.0, 0.5, 0.4), ecolor=(0.5, 0.0, 0.5), capsize=15)  # A soft serenity color theme
        graph_name = "Exploring the Relationship Between Enrichment Levels and Isotopic Composition" + " " + self.MSDI.date_data_was_taken
        for i, txt in enumerate(self.error.enrichments.tolist()):
            plt.annotate(txt, (x_axis[i], y_axis[i]))
        if not save_pdf:
            plt.show()

    def αVIs(self, save_pdf=False):
        """Generates a Alpha Vs Isotope graph from the data"""
        # Purely for Decoration
        if not self.professional:
            Color_Themes_Picker = np.array([0, 1, 2, 3])
            Random_Theme = np.random.choice(Color_Themes_Picker)
            Line_Colors = [(0.0, 0.2, 0.8), '#6dbb22', '#d52d00', '#3b9d1e']
            Dot_Colors = [(0.7, 0.5, 0.9), '#8c5c2f', '#0033cc', '#2e8b57']
            Error_Bar_colors = ['green', '#4a7023', '#003399', '#001f4d']
        Column_labels = ['Isotope', 'α']
        x_axis = self.error.Isotopes
        custom_labels = []
        for isotope in self.error.Isotopes:
            isotope_label = f"$\\mathrm{{^{{{isotope}}}{self.MSDI.atom_name}}}$"
            custom_labels.append(isotope_label)  # Custom labels
        y_axis = self.error.alphas
        yerr = self.error.δα
        plt.plot(x_axis, y_axis, color='#1a2b6d') if self.professional else plt.plot(x_axis, y_axis, color=Line_Colors[Random_Theme])
        plt.xticks(x_axis, custom_labels)
        plt.xlabel('Isotope')
        plt.ylabel('α')
        plt.title("Alpha Coefficient Trends as a Function of Isotopic Composition" + " " + self.MSDI.date_data_was_taken)
        plt.errorbar(x_axis, y_axis, yerr=yerr, fmt='s', color='#8b0000', ecolor='#888888', capsize=15) if self.professional else plt.errorbar(x_axis, y_axis, yerr=yerr, fmt='s', color=Dot_Colors[Random_Theme], ecolor=Error_Bar_colors[Random_Theme], capsize=15)
        graph_name = "Alpha Coefficient Trends as a Function of Isotopic Composition" + " " + self.MSDI.date_data_was_taken
        for i, txt in enumerate(self.error.alphas.tolist()):
            plt.annotate(txt, (x_axis[i], y_axis[i]))
        if not save_pdf:
            plt.show()

    def RIVM(self, PDF_Maker=None):
        """Relative intensity vs masses for parent and daughter graphs"""
        # Purely for Decoration
        if not self.professional:
            Color_Themes_Picker_1 = np.array([0, 1, 2, 3, 4, 5, 6, 7])
            Color_Themes_Picker_2 = np.array([0, 1, 2, 3, 4, 5, 6, 7])
            Random_Theme_1 = np.random.choice(Color_Themes_Picker_1)
            Random_Theme_2 = np.random.choice(Color_Themes_Picker_2)
            Line_Colors_1 = [(0.7, 0.1, 0.1), '#a23b3b', '#d16002', '#ff6f61', '#FFD700', "#8C1515", "#FFA500", "#D4A017"]
            Dot_Colors_1 = [(1.0, 0.85, 0.0), '#d9a673', '#9e2a2b', '#f28a30', '#FF6F61', "#B3995D", "#FF69B4", "#3366CC"]
            Error_Bar_colors_1 = [(1.0, 0.6, 0.2), '#ffb087', '#b17a34', '#f7c6a3', '#FFB6A0', "#4D4F53", "#FFD700", "#FFCC33"]
            Line_Colors_2 = [(0.2, 0.4, 0.2), '#2a4d69', '#006d77', '#191970', '#2e5d55', '#0f4c81', "#003262", "#0099FF", "#87CEEB"]
            Dot_Colors_2 = [(0.5, 0.25, 0.1), '#81894e', '#ff6b6b', '#7a7bb0', '#4682b4', '#84c0c6', "#FDB515", "#CC6600", "#8B4513"]
            Error_Bar_colors_2 = [(0.6, 1.0, 0.6), '#4b8e8d', '#83c5be', '#a4c3e2', '#b0c4de', '#d8e2dc', "#BDC3C7", "#66BB66", "#228B22"]
        # Actual plotting
        plt.figure(figsize=(10, 5))
        Column_labels = ['Mass', 'Relative Intensity']
        x_axis = self.error.masses_P  # rounded
        y_axis = self.error.parent_total_areas.tolist()
        yerr = self.error.Error_P
        plt.plot(x_axis, y_axis, color='#1a2b6d') if self.professional else plt.plot(x_axis, y_axis, color=Line_Colors_1[Random_Theme_1])
        ratio_label_P = self.error.Parent_ratios.tolist()
        plt.xlabel('Mass')
        plt.ylabel('Relative Intensity')
        plt.title("Relative Intensity by Mass in Parent Molecule Fragmentation" + " " + self.MSDI.date_data_was_taken)
        plt.errorbar(x_axis, y_axis, yerr=yerr, fmt='s', color='#8b0000', ecolor='#888888', capsize=15) if self.professional else plt.errorbar(x_axis, y_axis, yerr=yerr, fmt='s', color=Dot_Colors_1[Random_Theme_1], ecolor=Error_Bar_colors_1[Random_Theme_1], capsize=15)
        graph_name = "Relative Intensity by Mass in Parent Molecule Fragmentation" + " " + self.MSDI.date_data_was_taken
        for i, txt in enumerate(self.error.Parent_ratios):
            plt.annotate(txt, (x_axis[i], y_axis[i]))
        if PDF_Maker is not None:
            PDF_Maker.pdf_pages.savefig(dpi=300)
            plt.close()
        else:
            plt.show()

        plt.figure(figsize=(10, 5))
        Column_labels = ['Mass', 'Relative Intensity']
        x_axis = self.error.masses_D  # rounded
        y_axis = self.error.daughter_total_areas.tolist()
        yerr = self.error.Error_D
        plt.plot(x_axis, y_axis, color='#1a2b6d') if self.professional else plt.plot(x_axis, y_axis, color=Line_Colors_2[Random_Theme_2])
        plt.errorbar(x_axis, y_axis, yerr=yerr, fmt='s', color='#8b0000', ecolor='#888888', capsize=15) if self.professional else plt.errorbar(x_axis, y_axis, yerr=yerr, fmt='s', color=Dot_Colors_2[Random_Theme_2], ecolor=Error_Bar_colors_2[Random_Theme_2], capsize=15)
        plt.xlabel('Mass')
        plt.ylabel('Relative intensity')
        plt.title("Relative Intensity by Mass in Daughter Molecule Fragmentation" + " " + self.MSDI.date_data_was_taken)
        graph_name = "Relative Intensity by Mass in Daughter Molecule Fragmentation" + " " + self.MSDI.date_data_was_taken
        for i, txt in enumerate(self.error.Daughter_ratios):
            plt.annotate(txt, (x_axis[i], y_axis[i]))
        if PDF_Maker is not None:
            PDF_Maker.pdf_pages.savefig(dpi=300)
            plt.close()
        else:
            plt.show()


def αVIsol(Isotope, constant_NCE='P', NCE_Value_Key_Parent=None, NCE_Value_Key_Daughter=None, save_pdf=False):
    """Plot Alpha vs. Isolation for a given Isotope and constant NCE (Parent or Daughter).

    Parameters:
        Isotope (str): The isotope key to focus on.
        constant_NCE (str): 'P' for constant Parent, 'D' for constant Daughter.
        NCE_Value_Key_Parent (str or float): The NCE value for the Parent.
        NCE_Value_Key_Daughter (str or float): The NCE value for the Daughter.
        save_pdf (bool): Whether to save the plot as a PDF.
    """
    x_vals = []  # Isolation keys (top-level keys in the dictionary)
    y_vals = []  # Alpha values
    y_errs = []  # Error bars

    # Iterate over the isolation keys
    for isolation, isotopes in state.All_alphas_dict_PD.items():
        if Isotope not in isotopes:
            raise ValueError(f"Isotope '{Isotope}' not found in the specified isolation '{isolation}'.")

        # Process alphas and errors for the specified isotope
        isotope_alphas = state.All_alphas_dict_PD[isolation][Isotope]
        isotope_errors = state.All_Errors_dict_PD[isolation][Isotope]

        if constant_NCE == 'P':
            # Parent is constant; gather alphas and errors for Daughters
            for parent, daughters in isotope_alphas.items():
                for daughter, alpha_value in daughters.items():
                    x_vals.append(isolation)
                    y_vals.append(alpha_value)
                    y_errs.append(isotope_errors[parent][daughter])
        elif constant_NCE == 'D':
            # Daughter is constant; gather alphas and errors for Parents
            for daughter, parents in isotope_alphas.items():
                for parent, alpha_value in parents.items():
                    x_vals.append(isolation)
                    y_vals.append(alpha_value)
                    y_errs.append(isotope_errors[daughter][parent])
        else:
            raise ValueError("Invalid constant_NCE value. Must be 'P' or 'D'.")

    # Determine legend label
    constant_label = f"{constant_NCE}_NCE"  # P_NCE or D_NCE
    legend_label = f"N = {Isotope}, {constant_label} = {NCE_Value_Key_Parent if constant_NCE == 'P' else NCE_Value_Key_Daughter}"

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.errorbar(
        x_vals, y_vals, yerr=y_errs, fmt='o', capsize=5, label=legend_label
    )
    plt.xlabel("Isolation")
    plt.ylabel("Alpha Value")
    plt.title("Alpha vs. Isolation")
    plt.xticks(rotation=45)
    plt.legend()
    plt.grid(True)

    custom_labels = []
    for isotope in state.All_Isotopes:
        isotope_label = f"$\\mathrm{{^{{{isotope}}}{state.Atom_Names}}}$"
        custom_labels.append(isotope_label)  # Custom labels
    x_axis = state.All_Isotopes
    plt.xticks(x_axis, custom_labels)

    # Save or display the plot
    if not save_pdf:
        plt.show()


def αVNCE(Isotope, x_axis_NCE, constant_NCE=None, save_pdf=False):
    """Generates an alpha vs normalized collision energy graph where the x_axis NCE is something that is changed and the constant_NCE is a NCE value the data shares.
    If constant_NCE value is equal to None, it will assume alpha is at a constant NCE value and graph the corresponding alpha vs NCE plot. If it is not none, it assumes the
    plotted alphas are not under the same NCE values."""
    return 3


def Combined_Data_Line_Plot_EVI(Normalization_Number=1, Normalization_Number_is_total_areas=True, save_pdf=False):
    """Makes an Enrichment vs Isotope graph of all combined data after combined the same data stored in the global variables.
    It does need a Normalization number which is an integer of the value you want to normalize the data by, put 1 if none. If you want to normalize by total areas,
    just put the integer 1 into the function, otherwise put Normalization_Number_is_total_areas = False. Put save_pdf = True if you are in the process of adding this chart to a pdf"""

    # Sums values across all the tuples resulting in one tuple of all the summed areas.
    summed_tuple_Parent = tuple(sum(values) for values in zip(*state.All_Parent_Areas.values()))
    summed_tuple_Daughter = tuple(sum(values) for values in zip(*state.All_Daughter_Areas.values()))

    # Turns the tuples into an array to get the ratios which normalizes it
    Parent_Ratios = np.array(summed_tuple_Parent) / np.array(summed_tuple_Parent).sum() if Normalization_Number_is_total_areas else np.array(summed_tuple_Parent) / Normalization_Number
    Daughter_Ratios = np.array(summed_tuple_Daughter) / np.array(summed_tuple_Daughter).sum() if Normalization_Number_is_total_areas else np.array(summed_tuple_Daughter) / Normalization_Number

    # Calculating the errors and enrichments
    Enrichments_of_Combined_Data = tuple((Daughter_Ratios / Parent_Ratios) - 1)
    Errors_of_Combined_Data = tuple((Daughter_Ratios / Parent_Ratios) * np.sqrt(
        (np.sqrt(np.array(summed_tuple_Daughter)) / np.array(summed_tuple_Daughter)) ** 2 +
        (np.sqrt(np.array(summed_tuple_Parent)) / np.array(summed_tuple_Parent)) ** 2))

    # Making the graph
    plt.errorbar(state.All_Isotopes, Enrichments_of_Combined_Data, yerr=Errors_of_Combined_Data, capsize=3, marker='o', linestyle='-', color=(0.6, 0.8, 1.0), ecolor=(0.5, 0.0, 0.5), markerfacecolor=(1.0, 0.5, 0.4), markeredgecolor=(0.9, 0.7, 0.3))
    plt.xlabel('Isotope')
    plt.ylabel('ε')
    plt.title("ε Trends vs. Isotopic Composition for Combined Data")
    plt.grid(True)
    custom_labels = []
    for isotope in state.All_Isotopes:
        isotope_label = f"$\\mathrm{{^{{{isotope}}}{state.Atom_Names}}}$"
        custom_labels.append(isotope_label)  # Custom labels
    x_axis = state.All_Isotopes
    plt.xticks(x_axis, custom_labels)
    if not save_pdf:
        plt.show()


def Combined_Data_Line_Plot_AVI(Normalization_Number=1, Normalization_Number_is_total_areas=True, save_pdf=False):
    """Makes an Alpha vs Isotope graph after combined the data stored in the global variables.
    It does need a Normalization number which is an integer of the value you want to normalize the data by, put 1 if none. If you want to normalize by total areas,
    just put the integer 1 into the function, otherwise put Normalization_Number_is_total_areas = False. Put save_pdf = True if you are in the process of adding this chart to a pdf"""

    # Sums values across all the tuples resulting in one tuple of all the summed areas.
    summed_tuple_Parent = tuple(sum(values) for values in zip(*state.All_Parent_Areas.values()))
    summed_tuple_Daughter = tuple(sum(values) for values in zip(*state.All_Daughter_Areas.values()))

    # Turns the tuples into an array to get the ratios which normalizes it
    Parent_Ratios = np.array(summed_tuple_Parent) / np.array(summed_tuple_Parent).sum() if Normalization_Number_is_total_areas else np.array(summed_tuple_Parent) / Normalization_Number
    Daughter_Ratios = np.array(summed_tuple_Daughter) / np.array(summed_tuple_Daughter).sum() if Normalization_Number_is_total_areas else np.array(summed_tuple_Daughter) / Normalization_Number

    # Calculating the errors and alphas of the newly combined data
    Alphas_of_Combined_Data = tuple(Daughter_Ratios / Parent_Ratios)
    Errors_of_Combined_Data = tuple((Daughter_Ratios / Parent_Ratios) * np.sqrt(
        (np.sqrt(np.array(summed_tuple_Daughter)) / np.array(summed_tuple_Daughter)) ** 2 +
        (np.sqrt(np.array(summed_tuple_Parent)) / np.array(summed_tuple_Parent)) ** 2))

    # Color themes for the Graph! Put """ """ around it if you want the regular default colors
    Line_Colors = ('#6dbb22', "#87CEEB")
    Dot_Colors = ('#8c5c2f', "#8B4513")
    Error_Bar_Colors = ('#4a7023', "#228B22")
    Color_Themes_Picker = np.array([0, 1])
    Random_Theme = np.random.choice(Color_Themes_Picker)

    # Making the graph
    plt.errorbar(state.All_Isotopes, Alphas_of_Combined_Data, yerr=Errors_of_Combined_Data, capsize=3, marker='o', linestyle='-', color=Line_Colors[Random_Theme], ecolor=Error_Bar_Colors[Random_Theme], markerfacecolor=Dot_Colors[Random_Theme])
    plt.xlabel('Isotope')
    plt.ylabel('α')
    plt.title("Combined Data α Trends as a Function of Isotopic Composition")
    plt.grid(True)
    custom_labels = []
    for isotope in state.All_Isotopes:
        isotope_label = f"$\\mathrm{{^{{{isotope}}}{state.Atom_Names}}}$"
        custom_labels.append(isotope_label)  # Custom labels
    x_axis = state.All_Isotopes
    plt.xticks(x_axis, custom_labels)
    if not save_pdf:
        plt.show()


def overlaid_line_plot_EVI(save_pdf=False):
    """Makes an graph of all enrichment vs isotope data. This assumes you have used all_calculations_and_graphs multiple times or calculate_δα_indv multiple times and now want to plot all the data
    Put save_pdf = True if you are in the process of adding this chart to a pdf"""
    for key in state.All_enrichments:
        plt.errorbar(state.All_Isotopes, state.All_enrichments[key], yerr=state.All_Errors[key], label=key, capsize=3, marker='o', linestyle='-')
    plt.xlabel('Isotope')
    plt.ylabel('ε')
    plt.title("ε Trends vs. Isotopic Composition for Each Date")
    plt.legend(title="Dates", fontsize='small', borderpad=0.3, labelspacing=0.2, loc='best')
    plt.grid(True)
    custom_labels = []
    for isotope in state.All_Isotopes:
        isotope_label = f"$\\mathrm{{^{{{isotope}}}{state.Atom_Names}}}$"
        custom_labels.append(isotope_label)  # Custom labels
    x_axis = state.All_Isotopes
    plt.xticks(x_axis, custom_labels)
    if not save_pdf:
        plt.show()


def overlaid_line_plot_AVI(save_pdf=False):
    """Makes an graph of all alpha vs isotope data. This assumes you have used all_calculations_and_graphs multiple times or calculate_δα_indv multiple times and now want to plot all the data
    Put save_pdf = True if you are in the process of adding this chart to a pdf"""
    for key in state.All_alphas:
        plt.errorbar(state.All_Isotopes, state.All_alphas[key], yerr=state.All_Errors[key], label=key, capsize=3, marker='o', linestyle='-')
    plt.xlabel('Isotope')
    plt.ylabel('α')
    plt.title("α Coefficient Trends as a Function of Isotopic Composition for Each Date")
    plt.legend(title="Dates", fontsize='small', borderpad=0.3, labelspacing=0.2, loc='best')
    plt.grid(True)
    custom_labels = []
    for isotope in state.All_Isotopes:
        isotope_label = f"$\\mathrm{{^{{{isotope}}}{state.Atom_Names}}}$"
        custom_labels.append(isotope_label)  # Custom labels
    x_axis = state.All_Isotopes
    plt.xticks(x_axis, custom_labels)
    if not save_pdf:
        plt.show()
