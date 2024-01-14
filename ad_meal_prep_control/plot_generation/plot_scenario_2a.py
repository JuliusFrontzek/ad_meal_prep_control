from ad_meal_prep_control.postprocessing import (
    PlotVarProperty,
    MPLProperties,
    PostProcessing,
)

post_processing = PostProcessing(
    result_directory="/home/julius/Projects/ad_meal_prep_control/results",
    scenario_name="Scenario_2a",
)

post_processing.plot(
    [
        (
            "Gas storage fill volume",
            {
                "v_gas_storage": PlotVarProperty(
                    mpl_properties=MPLProperties(linewidth=1.5, linestyle="-")
                )
            },
        ),
        (
            r"$\dot V_{CH_4}$" + "\n" + r"$[m^3/d]$",
            {
                "v_ch4_dot_tank_in": PlotVarProperty(
                    mpl_properties=MPLProperties(linewidth=1.5, linestyle="-"),
                    label="",
                ),
            },
        ),
        (
            r"$\dot V_g$" + "\n" + r"$[m^3/d]$",
            {
                "y_1": PlotVarProperty(
                    mpl_properties=MPLProperties(linewidth=1.5, linestyle="-"), label=""
                )
            },
        ),
        (
            r"$pH$" + "\n" + r"$[1]$",
            {
                "y_4": PlotVarProperty(
                    mpl_properties=MPLProperties(linewidth=1.5, linestyle="-")
                )
            },
        ),
    ],
    plot_save_name="scenario_2a",
)
