# isotope-rel-effects

This is a repository for the code I developed during my research in the Heavy Elements Group at the Lawrence Berkeley National Lab under the direct supervision of Dr. Jennifer Pore. 

## About the code

There will be two aspects to this repository:
1. Code I developed for the project in 2024
2. Code I added or refactored after 2024

`fall_2024_isotope_project.py` and `Foundation.ipynb` contain the code I originally developed during my time in the lab. I am leaving it untouched in the `fall_2024` folder. 

The `fall_2024 folder also contains ` the flow chart I made showing how that code worked (in `fall_2024/charts`), `fall_2024/outputs/` which contains example PDFs produced by the original code, included to show what the analyses look like on the produced data, and `run_logs/`, which contains meta_data and notes about the sample collections we did.


Everything else is new work. The refactored code lives in `src/`.

I am adding new code because the original code was made when I was new to coding, so some choices in how I organized it were not the best. The notebook in particular became laggy, with processes stacking on each other and time complexity growing accordingly.

The underlying work itself is very useful and relevant for research, but the organization makes it impractical to build on. Now that I am much more experienced, I am going to fix the code to function the same as before, but now with a better structure (it would be a shame to let this code go to waste)

## A note on data

The mass spectra themselves are not included. They belong to the Heavy Elements Group and remain part of ongoing research at LBNL. What is included, in `run_logs/`, are the acquisition records for each file: isolation settings, collision energies, sample counts, and notes taken during measurement. The code is written to run on spectra in the original format, and the outputs shown here were produced from those datasets.

## Theory

I'll write the detailed overview of the theory later on. 

The main molecules this part of the experiment were done with were 

$[Nd(NO_{3})_{4}]^{-}$ and $[Dy(NO_{3})_{4}]^{-}$

