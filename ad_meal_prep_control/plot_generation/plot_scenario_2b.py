from ad_meal_prep_control.postprocessing import (
    PlotVarProperty,
    MPLProperties,
    PostProcessing,
    Constraint,
)
from ad_meal_prep_control import params_R3
import sys

post_processing = PostProcessing(
    result_directory="/home/julius/Projects/ad_meal_prep_control/results",
    scenario_name="Scenario_2b_dynamic",
)

try:
    dpi = int(sys.argv[1])
    show_plot = int(sys.argv[2])
except IndexError:
    dpi = 600
    show_plot = True

post_processing.plot(
    [
        (
            r"$\dot d_{feed}$" + "\n" + r"$[m^3/d]$",
            {
                "dictated_sub_feed_1": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color=post_processing.MEASUREMENT_COLOR,
                        linewidth=1.5,
                        linestyle="-",
                    ),
                    label="cattle manure\nwith large uncertainty",
                ),
            },
        ),
        (
            r"$x$" + "\n" + r"$[1]$",
            {
                "x_19": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="green", linewidth=1.5, linestyle="-"
                    ),
                    label=r"$V_{CH_4,tank}$",
                ),
                "x_20": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="black", linewidth=1.5, linestyle="-"
                    ),
                    label=r"$V_{CO_2,tank}$",
                ),
            },
        ),
        (
            r"$V_{g, tank}$" + "\n" + r"$[m^3]$",
            {
                "v_gas_storage": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color=post_processing.MEASUREMENT_COLOR,
                        linewidth=1.5,
                        linestyle="-",
                    )
                )
            },
        ),
        (
            r"$\dot V_{CH_4,AD}$" + "\n" + r"$[m^3/d]$",
            {
                "v_ch4_dot_tank_in": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color=post_processing.MEASUREMENT_COLOR,
                        linewidth=1.5,
                        linestyle="-",
                    ),
                    label="",
                ),
            },
        ),
        (
            r"$\dot V_{g, AD}$" + "\n" + r"$[m^3/d]$",
            {
                "y_1": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color=post_processing.MEASUREMENT_COLOR,
                        linewidth=1.5,
                        linestyle="-",
                    ),
                    label="",
                )
            },
        ),
        (
            r"$pH$" + "\n" + r"$[1]$",
            {
                "y_4": PlotVarProperty(
                    mpl_properties=MPLProperties(linewidth=1.5, linestyle="-"),
                    label="",
                )
            },
        ),
    ],
    plot_save_name="scenario_2b",
    constraints=[
        Constraint(value=0.0, ax_idx=2),
        Constraint(value=0.0, ax_idx=3),
        Constraint(value=params_R3.V_GAS_STORAGE_MAX, ax_idx=3),
    ],
    dpi=dpi,
    show_plot=show_plot,
)
