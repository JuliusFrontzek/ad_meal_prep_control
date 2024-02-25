from ad_meal_prep_control.postprocessing import (
    PlotVarProperty,
    MPLProperties,
    PostProcessing,
)
import sys


plot_save_name = "Scenario_1b_comparison"

try:
    dpi = int(sys.argv[1])
    show_plot = int(sys.argv[2])
except IndexError:
    dpi = 600
    show_plot = True

post_processing = PostProcessing(
    result_directory="/home/julius/Projects/ad_meal_prep_control/results",
    scenario_names=["Scenario_1b_non_robust", "Scenario_1b_multi_stage"],
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
                        label=r"$N_r=0$",
                    ),
                    PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="orange",
                            linewidth=1.5,
                            linestyle="-",
                        ),
                        label=r"$N_r=1$",
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
