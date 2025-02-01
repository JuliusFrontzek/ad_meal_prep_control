from ad_meal_prep_control.postprocessing import (
    PlotVarProperty,
    MPLProperties,
    PostProcessing,
)
import matplotlib.pyplot as plt
import sys


plot_save_name = "Scenario_1a_comparison"

try:
    dpi = int(sys.argv[1])
    show_plot = int(sys.argv[2])
except IndexError:
    dpi = 600
    show_plot = True

post_processing = PostProcessing(
    result_directory="/home/julius/Projects/ad_meal_prep_control/results",
    scenario_names=["Scenario_1a_linear", "Scenario_1a_quadratic"],
)

post_processing.plot(
    [
        (
            r"$\dot V_{CH_4}$" + "\n" + r"$[m^3/d]$",
            {
                "v_ch4_dot_tank_in_setpoint": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="orange", linewidth=1.5, linestyle="-."
                    ),
                    label="Reference",
                ),
                "v_ch4_dot_tank_in": [
                    PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="green",
                            linewidth=1.5,
                            linestyle="-",
                        ),
                        label="linear substrate cost",
                    ),
                    PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="orange",
                            linewidth=1.5,
                            linestyle="-",
                        ),
                        label="quadratic substrate cost",
                    ),
                ],
            },
        ),
    ],
    plot_inputs=False,
    plot_save_name=plot_save_name,
    time_range=(0.0, 15.0),
    dpi=dpi,
    show_plot=show_plot,
)
