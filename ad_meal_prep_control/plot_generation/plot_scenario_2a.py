from ad_meal_prep_control.postprocessing import (
    PlotVarProperty,
    MPLProperties,
    PostProcessing,
    Constraint,
)
from ad_meal_prep_control import params_R3
import sys
from ad_meal_prep_control.plot_generation.Controller_output_plotting.Output_Scenario_2a import controller_plotting_2a

scenario_names = ["Scenario_2a_dynamic_nominal_feedback_mismatch_3std_8tsap",
                  'Scenario_2a_dynamic_nominal_ideal_feedback_8tsap',
                  'Scenario_2a_dynamic_robust_feedback_mismatch_3std_8tsap',
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
        dpi = 600
        show_plot = False

    if not 'nominal_feedback' in scenario_name:
        post_processing.plot(
            [
                #(
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
                #),
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
                            mpl_properties=MPLProperties(linewidth=1, linestyle="-",
                                                     color="black"),
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
            plot_save_name=plot_save_name,
            constraints=[
                #Constraint(value=0.0, ax_idx=1),
                Constraint(value=0.0, ax_idx=1),
                Constraint(value=params_R3.V_GAS_STORAGE_MAX, ax_idx=1),
                Constraint(value=0.05 * params_R3.V_GAS_STORAGE_MAX, ax_idx=1, color="grey"),
                Constraint(value=0.95 * params_R3.V_GAS_STORAGE_MAX, ax_idx=1, color="grey"),
                Constraint(value=0, ax_idx=4, color="white"),
            ],
            dpi=dpi,
            show_plot=show_plot,
            height_ratios=[2, 2, 2, 1, 1],
        )
    else:
        post_processing.plot(
            [
                #(
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
                #),
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
                            mpl_properties=MPLProperties(linewidth=1, linestyle="-",
                                                     color="black"),
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
            plot_save_name=plot_save_name,
            constraints=[
                #Constraint(value=0.0, ax_idx=1),
                Constraint(value=0.0, ax_idx=1),
                Constraint(value=params_R3.V_GAS_STORAGE_MAX, ax_idx=1),
                Constraint(value=0.05 * params_R3.V_GAS_STORAGE_MAX, ax_idx=1, color="grey"),
                Constraint(value=0.95 * params_R3.V_GAS_STORAGE_MAX, ax_idx=1, color="grey"),
            ],
            other_inset_axes=[
                #{
                #    "plot_idx": 1,
                #    "days": (18, 23),
                #    "ylimit": (-10, 120),
                #    "inset_axis_specs": (0.4, 0.2, 0.35, 0.35)
                #},
                {
                    "plot_idx": 1,
                    "days": (18, 23),
                    "ylimit": (20, 300),
                    "inset_axis_specs": (0.4, 0.2, 0.35, 0.35)
                },
                #{
                #    "plot_idx": 3,
                #    "days": (3, 5),
                #    "ylimit": (100, 700),
                #    "inset_axis_specs": (0.4, 0.4, 0.25, 0.25)
                #},
                # {
                #    "plot_idx": 3,
                #    "days": (17.1, 17.4),
                #    "ylimit": (350, 500),
                #    "inset_axis_specs": (0.2, 0.4, 0.25, 0.25)
                # },
            ],
            dpi=dpi,
            show_plot=show_plot,
            height_ratios=[2, 2, 2, 1, 1],
        )
controller_plotting_2a()