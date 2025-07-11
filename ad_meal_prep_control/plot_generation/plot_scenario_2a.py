from ad_meal_prep_control.postprocessing import (
    PlotVarProperty,
    MPLProperties,
    PostProcessing,
    Constraint,
)
from ad_meal_prep_control import params_R3
import sys
from ad_meal_prep_control.plot_generation.Controller_output_plotting.Output_Scenario_2a import (
    controller_plotting_2a,
)

scenario_names = [
    "Scenario_2a_dynamic_nominal_feedback_mismatch_2std_12tsap",
    # "Scenario_2a_dynamic_nominal_ideal_feedback_8tsap",
    "Scenario_2a_dynamic_robust_feedback_mismatch_2std_12tsap",
]
plot_save_names = scenario_names

for scenario_name, plot_save_name in zip(scenario_names, plot_save_names):
    post_processing = PostProcessing(
        result_directory="../scenarios/results",
        scenario_name=plot_save_name,
    )

    try:
        dpi = int(sys.argv[1])
        show_plot = int(sys.argv[2])
    except IndexError:
        dpi = 1000
        show_plot = False

    post_processing.plot(
        [
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
                            color="black",
                            linewidth=1,
                            linestyle="-",
                        ),
                        label=r"$\dot V_{GS}$",
                    )
                },
            ),
            (
                r"$\dot V \; [m^3/d]$",
                {
                    "v_ch4_dot_tank_in": PlotVarProperty(
                        mpl_properties=MPLProperties(
                            color="blue",
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
                            linewidth=1, linestyle="-", color="black"
                        ),
                        label="$pH$",
                    )
                },
            ),
            # (
            #     r"$Inhibition \; [-]$",
            #     {
            #         "inhibition_1": PlotVarProperty(
            #             mpl_properties=MPLProperties(
            #                 color="black", linewidth=1, linestyle="-"
            #             ),
            #             label=r"$Inhibition_1$",
            #         ),
            #         "inhibition_2": PlotVarProperty(
            #             mpl_properties=MPLProperties(
            #                 color="mediumblue", linewidth=1, linestyle="-"
            #             ),
            #             label=r"$Inhibition_2$",
            #         ),
            #         "inhibition_3": PlotVarProperty(
            #             mpl_properties=MPLProperties(
            #                 color="cornflowerblue", linewidth=1, linestyle="-"
            #             ),
            #             label=r"$Inhibition_3$",
            #         ),
            #     },
            # ),
        ],
        plot_save_name=plot_save_name,
        constraints=[
            # Gas storage soft and hard constraints:
            Constraint(value=0.0, ax_idx=2),
            Constraint(value=params_R3.V_GAS_STORAGE_MAX, ax_idx=2),
            Constraint(value=0.05 * params_R3.V_GAS_STORAGE_MAX, ax_idx=2, color="grey"),
            Constraint(value=0.95 * params_R3.V_GAS_STORAGE_MAX, ax_idx=2, color="grey"),
            # adapt ylim of plots by adding invisible horizontal lines:
            Constraint(value=0, ax_idx=3, color="white"),  # gas production lower bound
            Constraint(value=7.6, ax_idx=4, color="white"),  # pH upper bound
            Constraint(value=7.2, ax_idx=4, color="white"),  # pH lower bound
            #Constraint(value=0, ax_idx=5, color="white"),  # inhibition lower bound
        ],
        dpi=dpi,
        show_plot=show_plot,
        height_ratios=[2, 1, 3, 3, 1],
        color_background_indices=(0,2,3),
        plot_olr=True,
    )
controller_plotting_2a(scenario_names)
