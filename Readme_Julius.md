# Description
Multi-stage Model Predictive Control algorithm and simulation for a biogas plant
(anaerobic digestion) with an optional external gas storage. Considers 
uncertainties in macronutrients of the fed substrates.

# File Structure
- ad_meal_prep_control: folder that contains the actual repo and all the python code based on do-mpc
- Matlab: contains Matlab-based code that was necessary to fix all the bugs in the AD model
- results: contains intermediate results Julius needed in the process of compiling his thesis
- tests: some test files

# Install
## Regular install
```
$ pip install .
```
## Dev install
```
$ pip install -e .
```

# Troubleshooting
## Plotting
If the Matplotlib based plot does not show up, try installing the GUI backend tk, e.g.
```
$ sudo apt-get install python3-tk
```
or refer to this Stackoverflow thread: https://stackoverflow.com/questions/56656777/userwarning-matplotlib-is-currently-using-agg-which-is-a-non-gui-backend-so