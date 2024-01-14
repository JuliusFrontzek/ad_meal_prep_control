from ad_meal_prep_control.postprocessing import (
    PlotVarProperty,
    MPLProperties,
    PostProcessing,
)

post_processing = PostProcessing(
    result_directory="/home/julius/Projects/ad_meal_prep_control/results",
    scenario_name="Scenario_2b",
)

post_processing.plot(
    [
        (
            "Forced substrates",
            {
                "dictated_sub_feed_1": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="black", linewidth=1.5, linestyle="-"
                    ),
                    label="Test",
                ),
            },
        ),
        (
            "Gas storage fill volume",
            {
                "v_gas_storage": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="red", linewidth=1.5, linestyle="-"
                    )
                )
            },
        ),
        (
            "Volume flow ch4",
            {
                "v_ch4_dot_tank_in": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="red", linewidth=1.5, linestyle="-"
                    )
                ),
            },
        ),
        (
            "V'g",
            {
                "y_1": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="red", linewidth=1.5, linestyle="-"
                    )
                )
            },
        ),
        (
            "pH",
            {
                "y_4": PlotVarProperty(
                    mpl_properties=MPLProperties(
                        color="red", linewidth=1.5, linestyle="-"
                    )
                )
            },
        ),
    ],
    plot_save_name="scenario_2b",
)
