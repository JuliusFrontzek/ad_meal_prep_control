# Readme
## Description

This project offers the possibility to simulate the control of a biogas plant using the ADM1-R3-Frac model as well as an optional gas storage model. Both a nominal and robust (multi-stage) control can be set up using a Model Predictive Control (MPC) scheme. The project is based on the Python library 'do-mpc' which allows us to robustly control the biogas plant while taking into account uncertain carbohydrate concentration in the substrates. It is based on the work of Julius Frontzek, who started this repo, and David M. Zarate both of whom were supervised by Simon Hellmann.

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

- Open the directory 'ad_meal_prep_control'
- Open the directory 'scenarios'
- Choose the case scenario of your preference, in the directory 'No feedback' you will find the sensitivity analysis script. The other ones include feedback as the default option
- Verify that the number of standard deviations in the controller_params and the kwargs are the same
- Run the script
- Once completed, open the directory 'plot_generation', choose the scenario you want to plot
- If you want to include the MPC predictions into the plot, open the directory 'Controller output plotting' and choose the scenario you want to plot
- Run the script 

## Additional information

- Executed on Microsoft Windows 10 Home Single Language Versión	10.0.19045 compilación 19045
- Executed using WSl and Ubuntu 22.04.5 LTS
- Executed using PyCharm Professional 2024.1.4
- Executed using Python 3.10.12 (on WSl and Ubuntu 22.04.5 LTS)
- Executed using IPopt 3.14.16 with solver MA27 of the Harwell Subroutine Library 

## Troubleshooting
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
