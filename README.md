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

<picture>
  <source media="(prefers-color-scheme: dark)"
          srcset="https://latex.codecogs.com/png.image?\bg_black%20%5BNd%28NO_%7B3%7D%29_4%5D%5E%7B-%7D" />
  <img src="https://latex.codecogs.com/png.image?\bg_white%20%5BNd%28NO_%7B3%7D%29_4%5D%5E%7B-%7D" alt="[Nd(NO3)4]-" />
</picture>
and
 <picture>
  <source media="(prefers-color-scheme: dark)"
          srcset="https://latex.codecogs.com/png.image?\bg_black%20%5BDy%28NO_%7B3%7D%29_4%5D%5E%7B-%7D" />
  <img src="https://latex.codecogs.com/png.image?\bg_white%20%5BDy%28NO_%7B3%7D%29_4%5D%5E%7B-%7D" alt="[Dy(NO3)4]-" />
</picture>

## Current Status/Plans

Current Status: I am refactoring the code that doesn't need to be changed or can be easily changed (getting rid of nested for loops) in python. So I'm in the process of doing 2.) below. I've already been working on this for a while before publically showing the changes, so I'll be done with 2, moving onto 3 pretty quickly. 

Plans:
1.) I want to hold off on writing the theory until I find my notes I took during this research, it'll make everything much more clear if I reference everything directly. So the theory to top off the repository itself will likely be added sometime Decemeber 2026).

2.) I will first refactor all the code in basic ways (for loops become vectorized operations, etc) and keep all functionality exactly the same

3.) After 2.), I will likely add new features and code some processes in Rust for a speedup where I think optimization would be useful (to any degree, it doesn't have to be the largest speedup, I just want to do this)

4.) I'll tie it all together as a repository of the newly refactored code, now usable for mass spectrometry data analysis of Isotope species.
