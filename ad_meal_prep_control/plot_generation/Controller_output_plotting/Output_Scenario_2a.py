import numpy as np
from do_mpc.data import save_results, load_results
import do_mpc
import matplotlib.pyplot as plt
import pickle as pkl
from ad_meal_prep_control.params_R3 import *
from ad_meal_prep_control.utils import remove_duplicate_labels
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl

mpl.rcParams['figure.dpi'] = 1000


def controller_plotting_2a(scenario_names=None):
    if scenario_names is None:
        scenario_names = ["Scenario_2a_dynamic_nominal_feedback_mismatch_3std_8tsap",
                          # 'Scenario_2a_dynamic_nominal_ideal_feedback_8tsap',
                          # 'Scenario_2a_dynamic_robust_feedback_mismatch_3std_8tsap'
                          ]
    for scenario in scenario_names:
        mpc = load_results(f'../scenarios/results/{scenario}_mpc_results.pkl')
        metadata = load_results(f'../scenarios/results/{scenario}_scenario_meta_data.pkl')
        with open(f'../scenarios/results/plots/{scenario}.pkl', 'rb') as file:
            fig = pkl.load(file)

        # modify Graph
        if metadata['feedback']:
            try:
                # shift the data by creating a NaN array
                predicted_data = np.genfromtxt(f'../scenarios/results/Predicted Data {scenario}.csv', delimiter=' ')
                data_shift = np.zeros((predicted_data.shape[0], metadata['t_stp_ahead_pred'] - 1))
                predicted_data = np.hstack((data_shift, predicted_data))
                predicted_data = predicted_data[:, 0:predicted_data.shape[1] - metadata['t_stp_ahead_pred']]
                # We remove the value, so they don't appear on the plot
                predicted_data[:, 0:metadata['t_stp_ahead_pred'] + 1] = np.NaN

                # __SH: create time vector:
                n_steps = predicted_data.shape[1]
                time_step = metadata['t_step']
                start_time = 0
                end_time = n_steps * time_step
                simulated_time = np.arange(start_time, end_time, time_step)
                predicted_time = simulated_time + metadata['t_stp_ahead_pred'] * time_step

                if mpc['mpc'].meta_data['n_robust'] == 0:  # __SH nominal MPC
                    # fig.axes[1].plot(simulated_time, predicted_data[0, :]*100, 'black',
                    #                 linestyle = 'dotted', label=r"$V_{CH_4,tank, MPC}$", linewidth=1)
                    # fig.axes[1].plot(simulated_time, predicted_data[1, :]*100, 'blue',
                    #                 linestyle = 'dotted', label=r"$V_{CO_2,tank, MPC}$", linewidth=1)
                    fig.axes[2].plot(simulated_time,
                                     predicted_data[2, :],
                                     'black', linestyle='dotted', label=r"$V_{GS,MPC}$",
                                     linewidth=1)  # controller predictions
                    # fig.axes[1].set_ylim(bottom=0, top=500)  # ylim of gas storage filling level
                    fig.axes[3].plot(simulated_time,
                                     predicted_data[3, :],
                                     'blue', linestyle='dotted', label=r"$\dot V_{CH_4,MPC}$", linewidth=1)
                    fig.axes[3].plot(simulated_time,
                                     predicted_data[4, :],
                                     'black', linestyle='dotted', label=r"$\dot V_{g,MPC}$", linewidth=1)
                    fig.axes[3].set_ylim(bottom=-50, top=1000)  # ylim of gas production
                    fig.axes[4].plot(simulated_time,
                                     predicted_data[5, :], 'black',
                                     linestyle='dotted', label=r'$pH_{MPC}$', linewidth=1)
                    fig.axes[4].set_ylim(bottom=4, top=12)  # ylim of pH
                    fig.axes[4].set_yticks(range(4, 12, 2))  # Set specific tick positions for pH

                    # plot insets
                    axins_1 = fig.axes[2].inset_axes([0.1, 0.4, 0.3, 0.4], xlim=(0, 3.7), ylim=(-10, 350))
                    axins_1.plot(simulated_time,
                                 predicted_data[2, :],
                                 'black', linestyle='dotted', linewidth=1)
                    axins_1.plot(simulated_time,
                                 mpc['mpc']['_aux', "v_gas_storage"],
                                 'black', linewidth=1)
                    # hard and soft constraints:
                    axins_1.plot(simulated_time,
                                 np.zeros((predicted_data.shape[1], 1)),
                                 'black', linestyle='dashed')  # __SH lower bound GS
                    axins_1.plot(simulated_time,
                                 np.ones((predicted_data.shape[1], 1)) * V_GAS_STORAGE_MAX,
                                 'black', linestyle='dashed')  # __SH upper bound GS
                    axins_1.plot(simulated_time,
                                 np.ones((predicted_data.shape[1], 1)) * V_GAS_STORAGE_MAX * (
                                     metadata['controller_params']['gas_storage_bound_fraction']),
                                 'grey', linestyle='dashed')  # __SH lower soft constraint GS
                    axins_1.plot(simulated_time,
                                 np.ones((predicted_data.shape[1], 1)) * V_GAS_STORAGE_MAX * (
                                         1 - metadata['controller_params']['gas_storage_bound_fraction']),
                                 'grey', linestyle='dashed')  # __SH lower soft constraint GS
                    axins_1.tick_params(axis='x', which='both', labelbottom=True)  # Ensure labels are shown
                    axins_1.grid(True, linestyle='dashed')  # show grid
                    fig.axes[2].indicate_inset_zoom(
                        axins_1,
                        edgecolor="black",
                        linewidth=1.0,
                    )

                    # axins_2 = fig.axes[2].inset_axes([0.55, 0.5, 0.24, 0.4], xlim=(22, 26), ylim=(50, 380))
                    # axins_2.plot(simulated_time,
                    #              predicted_data[2, :],
                    #              'black', linewidth=1, linestyle='dotted')
                    # # SH: this is false, cause these are the mpc predictions (only without time step ahead predictions):
                    # # axins_2.plot(simulated_time, mpc['mpc']['_aux', "v_gas_storage"],
                    # #             'black', linewidth=1)
                    # # SH: this is correct, cause this is the plant performance:
                    # axins_2.plot(simulated_time,
                    #              mpc['simulator']._aux[:,26],
                    #              'black', linewidth=1)
                    # # hard and soft constraints:
                    # axins_2.plot(simulated_time,
                    #              np.ones((predicted_data.shape[1], 1)) * V_GAS_STORAGE_MAX,
                    #              'black', linestyle='dashed')
                    # axins_2.plot(simulated_time,
                    #              np.ones((predicted_data.shape[1], 1)) * V_GAS_STORAGE_MAX * (
                    #                      1 - metadata['controller_params']['gas_storage_bound_fraction']),
                    #              'grey', linestyle='dashed')
                    # axins_2.tick_params(axis='x', which='both', labelbottom=True)  # Ensure labels are shown
                    # axins_2.grid(True, linestyle='dashed')  # show grid
                    # fig.axes[2].indicate_inset_zoom(
                    #     axins_2,
                    #     edgecolor="black",
                    #     linewidth=1.0,
                    # )

                    # axins_3.plot(simulated_time, predicted_data[3, :], 'rebeccapurple',
                    #                 linestyle = 'dotted', linewidth=1)
                    # axins_3.plot(simulated_time, predicted_data[4, :], 'blue',
                    #                 linestyle = 'dotted', linewidth=1)

                if mpc['mpc'].meta_data['n_robust'] > 0:
                    # fig.axes[1].plot(simulated_time, predicted_data[[0,1], :].transpose()*100, 'black',
                    #                 linestyle = 'dotted', label=r"$V_{CH_4,tank, MPC}$", linewidth=1)
                    # fig.axes[1].plot(simulated_time, predicted_data[[2,3], :].transpose()*100, 'blue',
                    #                 linestyle = 'dotted', label=r"$V_{CO_2,tank, MPC}$", linewidth=1)
                    fig.axes[2].plot(simulated_time,
                                     predicted_data[[4, 5], :].transpose(), 'black',
                                     linestyle='dotted', label=r"$V_{GS,MPC}$", linewidth=1)
                    fig.axes[3].plot(simulated_time,
                                     predicted_data[[6, 7], :].transpose(), 'blue',
                                     linestyle='dotted', label=r"$\dot V_{CH_4,MPC}$", linewidth=1)
                    fig.axes[3].plot(simulated_time,
                                     predicted_data[[8, 9], :].transpose(), 'black',
                                     linestyle='dotted', label=r"$\dot V_{g,MPC}$", linewidth=1)
                    fig.axes[4].plot(simulated_time,
                                     predicted_data[[10, 11], :].transpose(), 'black',
                                     linestyle='dotted', label=r'$pH_{MPC}$', linewidth=1)

            except:
                pass

        if not metadata['feedback']:
            plant_output = np.genfromtxt(f'../scenarios/results/Plant Output {scenario}.csv', delimiter=' ')
            plant_output = plant_output[:-1, :]
            plant_output[0, :] = np.NaN

            fig.axes[1].plot(np.linspace(0, metadata['n_days_mpc'], num=plant_output.shape[0]),
                             plant_output[:, 0] * 100, 'limegreen',
                             linestyle='dotted', label=r"$V_{CH_4,tank, MPC}$", linewidth=1)
            fig.axes[1].plot(np.linspace(0, metadata['n_days_mpc'], num=plant_output.shape[0]),
                             plant_output[:, 1] * 100, 'darkgreen',
                             linestyle='dotted', label=r"$V_{CO_2,tank, MPC}$", linewidth=1)
            fig.axes[2].plot(np.linspace(0, metadata['n_days_mpc'], num=plant_output.shape[0]), plant_output[:, 2],
                             'rebeccapurple',
                             linestyle='dotted', label=r"$V_{g, tank, MPC}$", linewidth=1)
            fig.axes[3].plot(np.linspace(0, metadata['n_days_mpc'], num=plant_output.shape[0]), plant_output[:, 3],
                             'rebeccapurple',
                             linestyle='dotted', label=r"$\dot V_{CH_4,AD, MPC}$", linewidth=1)
            fig.axes[3].plot(np.linspace(0, metadata['n_days_mpc'], num=plant_output.shape[0]), plant_output[:, 4],
                             'blue',
                             linestyle='dotted', label=r"$\dot V_{g, AD, MPC}$", linewidth=1)
            fig.axes[4].plot(np.linspace(0, metadata['n_days_mpc'], num=plant_output.shape[0]), plant_output[:, 5],
                             'rebeccapurple',
                             linestyle='dotted', label=r'$pH_{model}$', linewidth=1)
            fig.axes[0].legend()
            fig.axes[1].legend(bbox_to_anchor=(0.95, 1))
            fig.axes[2].legend()
            fig.axes[3].legend(bbox_to_anchor=(0.95, 1))
            fig.axes[4].legend()
            plt.savefig(f'../scenarios/results/plots/Plant Output {scenario}.png')

        fig.axes[2].legend()
        fig.axes[3].legend()
        fig.axes[4].legend()

        #remove_duplicate_labels(fig, 1, legend_location='upper left', bbox_to_anchor=(0, 1))
        remove_duplicate_labels(fig, 2, legend_location='upper right')#, bbox_to_anchor=(0, 0.6))
        remove_duplicate_labels(fig, 3, legend_location='center right')#, bbox_to_anchor=(0, 0.6))
        remove_duplicate_labels(fig, 4, legend_location='center right')#, bbox_to_anchor=(0, 0.6))
        #remove_duplicate_labels(fig, 5, legend_location='center left', bbox_to_anchor=(0, 0.6))
        plt.savefig(f'../scenarios/results/plots/Scenario 2/{scenario}.png')
    return
