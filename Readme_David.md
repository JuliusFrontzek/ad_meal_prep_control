# Readme
By David M. Zarate
## This is a How-to for running the scripts.

### Description

This project offers the possibility to simulate the control of a biogas plant using the ADM1-R3-Frac model. It is based on the work of Frotzek (2024), who started this repo

### Modifications since Julius' master's thesis

- Plant-model mismatch is modified in order to have carbohydrates as the only uncertain parameter
- Implementation of an extra simulator Simulator_plant for sensitivity analysis
- Addition of MPC predictions into the plots
- New scenarios for ideal, nominal and robust MPC

### Instructions

- Open the folder 'ad_meal_prep_control'
- Open the folder 'scenarios'
- Choose the case scenario of your preference, in the folder 'No feedback' you will find the sensitivity analysis script. The other ones include feedback as the default option
- Verify that the number os standard deviations in the controller_params and the kwargs are the same
- Run the script
- Once completed, open the foler 'plot_generation', choose the scenario you want to plot
- If you want to include the MPC predictions into the plot, open the folder 'Controller output plotting' and choose the scenario you want to plot
- Run the script 

## Additional information

- Exectued on Microsoft Windows 10 Home Single Language Versión	10.0.19045 compilación 19045
- Executed using WSl and Ubuntu 22.04.5 LTS
- Executed using PyCharm Professional 2024.1.4
- Executed using Python 3.10.12 (on WSl and Ubuntu 22.04.5 LTS)
- Executed using solver MA27 of IPopt 3.14.16 
