from ad_meal_prep_control.postprocessing import (
    PlotVarProperty,
    MPLProperties,
    PostProcessing,
    Constraint,
)
from ad_meal_prep_control import params_R3
import sys

post_processing = PostProcessing(
    result_directory="../scenarios/results",
    scenario_name="Scenario_2c_dynamic",
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
            r"$d_{feed}$" + "\n" + r"$[m^3/d]$",
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
            "Gas storage  \nlevel" + r" $[\%]$",
            {
                "x_19": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="green", linewidth=1.5, linestyle="-"
                    ),
                    label=r"$V_{CH_4,tank, plant}$",
                    scaling=100.0,
                ),
                "x_20": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="black", linewidth=1.5, linestyle="-"
                    ),
                    label=r"$V_{CO_2,tank, plant}$",
                    scaling=100.0,
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
                    ),
                    label=r"$\dot V_{g, tank, plant}$",
                )
            },
        ),
        (
            r"$\dot V$" + "\n" + r"$[m^3/d]$",
            {
                "v_ch4_dot_tank_in": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color=post_processing.MEASUREMENT_COLOR,
                        linewidth=1.5,
                        linestyle="-",
                    ),
                    label=r"$\dot V_{CH_4,AD, plant}$",
                ),
                "y_1": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="blue",
                        linewidth=1.5,
                        linestyle="-",
                    ),
                    label=r"$\dot V_{g, AD, plant}$",
                ),
            },
        ),
        (
            r"$pH$" + "\n" + r"$[-]$",
            {
                "y_4": PlotVarProperty(
                    mpl_properties=MPLProperties(linewidth=1.5, linestyle="-"),
                    label="$pH_{plant}$",
                )
            },
        ),
        (
            r"$Inhibition$" + "\n" + r"$[-]$",
            {
                "inhibition_1": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="blue", linewidth=1, linestyle="-"
                    ),
                    label=r"$Inhibition_1$",
                ),
                "inhibition_2": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="green", linewidth=1, linestyle="-"
                    ),
                    label=r"$Inhibition_2$",
                ),
                "inhibition_3": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="red", linewidth=1, linestyle="-"
                    ),
                    label=r"$Inhibition_3$",
                ),
            },
        ),
    ],
    plot_save_name="scenario_2c",
    constraints=[
        Constraint(value=0.0, ax_idx=2),
        Constraint(value=0.0, ax_idx=3),
        Constraint(value=params_R3.V_GAS_STORAGE_MAX, ax_idx=3),
        Constraint(value=0.05 * params_R3.V_GAS_STORAGE_MAX, ax_idx=3, color="blue"),
        Constraint(value=0.95 * params_R3.V_GAS_STORAGE_MAX, ax_idx=3, color="blue"),
    ],
    dpi=dpi,
    show_plot=show_plot,
    height_ratios=[2, 1, 2, 2, 2, 1, 1],
)
