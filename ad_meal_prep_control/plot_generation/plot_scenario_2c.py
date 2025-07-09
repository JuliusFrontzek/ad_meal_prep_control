from ad_meal_prep_control.postprocessing import (
    PlotVarProperty,
    MPLProperties,
    PostProcessing,
    Constraint,
)
from ad_meal_prep_control import params_R3
import sys
from ad_meal_prep_control.plot_generation.Controller_output_plotting.Output_Scenario_2c import controller_plotting_2c

scenario_name = "Scenario_2c_dynamic"
# scenario_name = "Scenario_2c_uninhibited"

post_processing = PostProcessing(
    result_directory="../scenarios/results",
    scenario_name=scenario_name,
)

try:
    dpi = int(sys.argv[1])
    show_plot = int(sys.argv[2])
except IndexError:
    dpi = 1000
    show_plot = False

post_processing.plot(
    [
        (
            r"$d_{feed} \; [m^3/d]$",
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
            r"$V_{GS} \; [m^3]$",
            {
                "v_gas_storage": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color='black',
                        linewidth=1,
                        linestyle="-",
                    ),
                    # label=r"$\dot V_{g, tank, plant}$",
                )
            },
        ),
        (
            r"$\dot V \; [m^3/d]$",
            {
                "v_ch4_dot_tank_in": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color='blue',
                        linewidth=1,
                        linestyle="-",
                    ),
                    label=r"$\dot V_{CH_4}$",
                ),
                "y_1": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="black",
                        linewidth=1,
                        linestyle="-",
                    ),
                    label=r"$\dot V_{g}$",
                ),
            },
        ),
        (
            r"$pH \; [-]$",
            {
                "y_4": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        linewidth=1,
                        linestyle="-",
                        color="black"
                    ),
                    # label="$pH_{plant}$",
                )
            },
        ),
        # (
        #     r"$Inhibition$" + "\n" + r"$[-]$",
        #     {
        #         "inhibition_1": PlotVarProperty(
        #             mpl_properties=MPLProperties(
        #                 color="black", linewidth=1, linestyle="-"
        #             ),
        #             label=r"$I_{N-lim}$",
        #         ),
        #         "inhibition_2": PlotVarProperty(
        #             mpl_properties=MPLProperties(
        #                 color="mediumblue", linewidth=1, linestyle="-"
        #             ),
        #             label=r"$I_{pH}$",
        #         ),
        #         "inhibition_3": PlotVarProperty(
        #             mpl_properties=MPLProperties(
        #                 color="cornflowerblue", linewidth=1, linestyle="-"
        #             ),
        #             label=r"$I_{nh3}$",
        #         ),
        #     },
        # ),
    ],
    plot_save_name=scenario_name,
    constraints=[
        Constraint(value=0.0, ax_idx=3),
        Constraint(value=params_R3.V_GAS_STORAGE_MAX, ax_idx=3),
        Constraint(value=0.05 * params_R3.V_GAS_STORAGE_MAX, ax_idx=3, color="grey"),
        Constraint(value=0.95 * params_R3.V_GAS_STORAGE_MAX, ax_idx=3, color="grey"),
        # adapt ylim of plots by adding invisible horizontal lines:
        Constraint(value=2, ax_idx=1, color="white"),    # dictated feeding upper bound
        Constraint(value=0, ax_idx=4, color="white"),    # gas production lower bound
        Constraint(value=7.5, ax_idx=5, color="white"),  # pH upper bound
        Constraint(value=7.3, ax_idx=5, color="white"),  # pH lower bound
        # Constraint(value=0, ax_idx=6, color="white"),    # inhibition lower bound
    ],
    dpi=dpi,
    show_plot=show_plot,
    height_ratios=[3, 1, 1, 3, 2, 1],
    # input_inset_axes=[
    #     {
    #         "days": (12, 13),
    #         "ylimit": (0, 1),
    #         "inset_axis_specs": (0.4, 0.7, 0.05, 0.15),
    #     }],
    # other_inset_axes=[
    #     {
    #         "plot_idx": 3,
    #         "days": (6, 7.5),
    #         "ylimit": (320, 340),
    #         "inset_axis_specs": (0.3, 0.35, 0.1, 0.2),
    #     },
    # ],
    color_background_indices=(0,3,4),
    plot_olr=True,
)
#controller_plotting_2c()
