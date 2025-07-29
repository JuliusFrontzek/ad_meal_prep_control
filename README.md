# Readme
## Short Info
With this project one can simulate the control of a biogas plant using the ADM1-R3-Frac model as well as an optional gas storage model. Both a nominal and robust (multi-stage) control can be set up using a Model Predictive Control (MPC) scheme. The project is based on the Python library 'do-mpc' which allows to robustly control the biogas plant while taking into account uncertain influent concentrations of the substrates. It is based on the work of Julius Frontzek, who started this repo, and David M. Zarate, both of whom were supervised by Simon Hellmann.

## File structure
- ad_meal_prep_control: here the actual code for simulation/optimization and plotting is saved
- data: contains input (in) and generated output (out) data required for computing inlet concentrations 
- data_preparation: contains computation of influent concentrations (2 versions: v2023 (with data used during Julius' MA), and v2025_overkill (data tried for do-MPC paper but then discarded))
- documentation: contains paper, notes, some further analyses worth documenting
- 

## Install

After cloning/downloading this project, you can install it as a Python library using one of the two following commands inside its main directory:

### Regular install
```
$ pip install .
```
### Dev install
```
$ pip install -e .
```

## Usage instructions
### Running simulations
- to run simulations, open '\ad_meal_prep_control\scenarios\' and choose the case scenario of your preference. The files should run as-are, provided the installation was done correctly.
- _note_: getting the HSL MA27 solver running on mac is challenging, see documentation directory!
- in the directory 'No feedback' you will find the sensitivity analysis script. The other ones include feedback as the default option

### Plotting
- for plotting, run the respective scripts in \ad_meal_prep_control\ad_meal_prep_control\plot_generation. 
- _note_: 
	- Verify that the number of standard deviations in the controller_params and the kwargs are the same
	- extensions to these plots were made by David, see scripts in subfolder \Controller_output_plotting\.
	- If you want to include the MPC predictions into the plot, open the directory 'Controller output plotting' and choose the scenario you want to plot
- the resulting plots will be stored in ad_meal_prep_control/scenarios/results/plots/{Scenario 1;Scenario 2; Sensitivity}. 
- further plots of the do-MPC paper with Bioresource Technology are also created from within \ad_meal_prep_control\ad_meal_prep_control\plot_generation 
	- _note:_ mind that additional plots separate from do-MPC's scenarios are saved in ad_meal_prep_control\results !

## Additional information
### Hard- and software specs for David's internship
- Executed on Microsoft Windows 10 Home Single Language Versión	10.0.19045 compilación 19045 
- Executed using WSl and Ubuntu 22.04.5 LTS
- Executed using PyCharm Professional 2024.1.4
- Executed using Python 3.10.12 (on WSl and Ubuntu 22.04.5 LTS)
- Executed using IPopt 3.14.16 with solver MA27 of the Harwell Subroutine Library 
### Hard- and software specs Simon's work
- macBook pro M1 (8GB RAM)
- Editor: PyCharm 2024.1.3 (Professional Edition)
- HSL MA27 solver: coinhsl-archive-2022.12.02

## Troubleshooting
### run code with HSL solvers
- David specified the exact location of the HSL solver in mpc.py, see 'nlpsol_opts', and then path under 'ipopt.hsllib'
- this might produce errors, then try leaving out the 'ipopt.hsllib' statement
- the solver is then searched in the default directory.

### Plotting
If the Matplotlib based plot does not show up, try installing the GUI backend tk, e.g.
```
$ sudo apt-get install python3-tk
```
or refer to this Stackoverflow thread: https://stackoverflow.com/questions/56656777/userwarning-matplotlib-is-currently-using-agg-which-is-a-non-gui-backend-so

## Modifications since Julius Frontzek's Master's thesis

- Plant-model mismatch is modified in order to have carbohydrates as the only uncertain parameter
- Implementation of an extra simulator Simulator_plant for sensitivity analysis
- Addition of MPC predictions into the plots
- New scenarios for ideal, nominal and robust MPC
