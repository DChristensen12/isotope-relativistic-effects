"""Chemical formula formatting helpers. These are just string/regex transforms, nothing here scales with data size."""

import re


def Get_chemical_formula_in_Latex(chemical_formula):
    """Takes a chemical formula and returns the LaTeX-formatted string. Handles
    subscripts, superscripts, and charges."""
    # Step 1: Handle groups with subscripts (e.g., NO3 -> NO_3, SO4 -> SO_4)
    formula = re.sub(r'([A-Za-z]+)(\d+)', r'\1_\2', chemical_formula)
    formula = re.sub(r'(NO3|SO4|CO3|PO4|Cl|Br|I|NH4)(\d*)', r'\1_\2', formula)

    # Step 2: Handle superscript for charges (e.g., - -> ^{-})
    formula = re.sub(r'\^-(?=\])', r'^{-}', formula)
    formula = re.sub(r'\^(\d+)', r'^{\1}', formula)

    # Step 3: Wrap the formula with \mathrm{} for proper chemical notation
    return r'$\mathrm{' + formula + r'}$'


def format_chemical_formula(chemical_formula):
    """Takes a LaTeX-style formula using _{...} and ^{...} for subscripts and
    superscripts and returns it with real Unicode sub/superscript characters.

    'H_{2}O' -> H₂O
    'NO_{3}^{-}' -> NO₃⁻
    'Nd(NO_{3})_{4}^{-}' -> Nd(NO₃)₄⁻
    '[Nd(NO_{3})_{4}]^{-}' -> [Nd(NO₃)₄]⁻
    """
    subscript_map = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    superscript_map = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")

    # Convert _{number} to subscripts
    formula = re.sub(r'_{(\d+)}', lambda m: m.group(1).translate(subscript_map), chemical_formula)

    # Convert ^{number} or ^{-} to superscripts
    formula = re.sub(r'\^{(-?\d+)}', lambda m: m.group(1).translate(superscript_map), formula)
    formula = re.sub(r'\^{(-)}', '⁻', formula)

    # Remove leftover curly braces (NO_{3} -> NO₃)
    formula = re.sub(r'\{|\}', '', formula)

    # Catch the ^- case that isn't wrapped in braces
    formula = formula.replace('^-', '⁻')

    return formula
