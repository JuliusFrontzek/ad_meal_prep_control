from ad_meal_prep_control.postprocessing import (
    PlotVarProperty,
    MPLProperties,
    PostProcessing,
    Constraint,
)
from ad_meal_prep_control import params_R3
import sys
from ad_meal_prep_control.plot_generation.Controller_output_plotting.Output_Scenario_2c import controller_plotting_2c

post_processing = PostProcessing(
    result_directory="../scenarios/results",
    scenario_name="Scenario_2c_dynamic",
)

try:
    dpi = int(sys.argv[1])
    show_plot = int(sys.argv[2])
except IndexError:
    dpi = 600
    show_plot = False

post_processing.plot(
    [
        (
            r"$d_{feed}$" + "\n" + r"$[m^3/d]$",
            {
                "dictated_sub_feed_1": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color='saddlebrown',
                        linewidth=1,
                        linestyle="-",
                    ),
                    label="cattle manure\nwith large uncertainty",
                ),
            },
        ),
        # (
        #    "Gas storage  \nlevel" + r" $[\%]$",
        #    {
        #        "x_19": PlotVarProperty(
        #            mpl_properties=MPLProperties(
        #                color="black", linewidth=1, linestyle="-"
        #            ),
        #            label=r"$V_{CH_4,tank, plant}$",
        #            scaling=100.0,
        #        ),
        #        "x_20": PlotVarProperty(
        #            mpl_properties=MPLProperties(
        #                color="blue", linewidth=1, linestyle="-"
        #            ),
        #            label=r"$V_{CO_2,tank, plant}$",
        #            scaling=100.0,
        #        ),
        #    },
        # ),
        (
            r"$V_{g, tank}$" + "\n" + r"$[m^3]$",
            {
                "v_gas_storage": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color='black',
                        linewidth=1,
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
                        color='blue',
                        linewidth=1,
                        linestyle="-",
                    ),
                    label=r"$\dot V_{CH_4,AD, plant}$",
                ),
                "y_1": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="black",
                        linewidth=1,
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
                    mpl_properties=MPLProperties(linewidth=1, linestyle="-", color="black"),
                    label="$pH_{plant}$",
                )
            },
        ),
        (
            r"$Inhibition$" + "\n" + r"$[-]$",
            {
                "inhibition_1": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="black", linewidth=1, linestyle="-"
                    ),
                    label=r"$Inhibition_1$",
                ),
                "inhibition_2": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="mediumblue", linewidth=1, linestyle="-"
                    ),
                    label=r"$Inhibition_2$",
                ),
                "inhibition_3": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="cornflowerblue", linewidth=1, linestyle="-"
                    ),
                    label=r"$Inhibition_3$",
                ),
            },
        ),
    ],
    plot_save_name="Scenario_2c_dynamic",
    constraints=[
        Constraint(value=0.0, ax_idx=2),
        Constraint(value=params_R3.V_GAS_STORAGE_MAX, ax_idx=2),
        Constraint(value=0.05 * params_R3.V_GAS_STORAGE_MAX, ax_idx=2, color="grey"),
        Constraint(value=0.95 * params_R3.V_GAS_STORAGE_MAX, ax_idx=2, color="grey"),
        # adapt ylim of plots by adding invisible horizontal lines:
        Constraint(value=0, ax_idx=3, color="white"),  # gas production lower bound
        Constraint(value=7.6, ax_idx=4, color="white"),  # pH upper bound
        Constraint(value=0, ax_idx=5, color="white"),  # inhibition lower bound
    ],
    dpi=dpi,
    show_plot=show_plot,
    height_ratios=[2, 1, 2, 2, 1, 1],
    other_inset_axes=[
        {
            "plot_idx": 2,
            "days": (6, 7.5),
            "ylimit": (200, 225),
            "inset_axis_specs": (0.12, 0.3, 0.1, 0.2),
        },
    ],
    color_background_indices=(2,),
)
controller_plotting_2c()
