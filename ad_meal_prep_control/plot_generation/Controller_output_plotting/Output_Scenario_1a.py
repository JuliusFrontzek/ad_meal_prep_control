import numpy as np
from do_mpc.data import save_results, load_results
import do_mpc
import matplotlib.pyplot as plt
import pickle as pkl
from ad_meal_prep_control.utils import remove_duplicate_labels, nRMSE


def controller_plotting_1a(scenario_names=None):
    if scenario_names is None:
        scenario_names = ['Scenario_1a_quadratic_nominal_feedback_mismatch_3std_3tsap',
                          'Scenario_1a_quadratic_nominal_ideal_feedback_3tsap',
                          'Scenario_1a_quadratic_robust_feedback_mismatch_3std_3tsap'
                          ]
    elif scenario_names in 'Sensitivity':
        scenario_names = ['Scenario_1a_quadratic_no_feedback_mismatch_1std_ch',
                          'Scenario_1a_quadratic_no_feedback_mismatch_1std_pr',
                          'Scenario_1a_quadratic_no_feedback_mismatch_1std_li']

    for scenario in scenario_names:
        with open(f'../scenarios/results/plots/{scenario}_substrate_costs.pkl', 'rb') as file:
            fig = pkl.load(file)
        mpc = load_results(f'../scenarios/results/{scenario}_mpc_results.pkl')
        metadata = load_results(f'../scenarios/results/{scenario}_scenario_meta_data.pkl')
        # Graph
        if not metadata['feedback']:
            plant_output = np.genfromtxt(f'../scenarios/results/Plant Output {scenario}.csv', delimiter=' ')

            fig.axes[2].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 1], 'black',
                             linestyle='dotted', label=r'$pH_{controller}$')
            fig.axes[1].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 0], 'black',
                             linestyle='dotted', label=r"$\dot V_{g, controller}$")
            fig.axes[1].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 2], 'blue',
                             linestyle='dotted', label=r"$\dot V_{CH_4, controller}$")

            error_ph = (f'{scenario}_NRMSE_pH = ', nRMSE(x_est=mpc['mpc']['_aux', 'y_4'], x_true=plant_output[:, [1]]))
            error_gas = (
                f'{scenario}_NRMSE_gas = ', nRMSE(x_est=mpc['mpc']['_aux', 'y_1'], x_true=plant_output[:, [0]]))
            error_ch4 = (f'{scenario}_NRMSE_ch4 = ',
                         nRMSE(x_est=mpc['mpc']['_aux', 'v_ch4_dot_tank_in'], x_true=plant_output[:, [2]]))

            results = [
                [f"{scenario}_NRMSE_pH", error_ph],
                [f"{scenario}_NRMSE_gas", error_gas],
                [f"{scenario}_NRMSE_ch4", error_ch4]
            ]

            output_file = f'../scenarios/results/plots/Sensitivity/{scenario}_NRMSE.txt'
            with open(output_file, 'w') as f:
                f.write("Metric, Value\n")  # Header for clarity
                for metric, value in results:
                    f.write(f"{metric}, {value}\n")

            remove_duplicate_labels(fig, 0, legend_location='center right', bbox_to_anchor=(1, 1))
            remove_duplicate_labels(fig, 1, legend_location='center right', bbox_to_anchor=(1, 0.45))
            remove_duplicate_labels(fig, 2, legend_location='center right', bbox_to_anchor=(1, 0.7))
            # remove_duplicate_labels(fig, 3, legend_location='center right', bbox_to_anchor=(1,0.55))
            fig.axes[0].legend(loc='lower right', bbox_to_anchor=(0.67, 0.6))
            plt.savefig(f'../scenarios/results/plots/Sensitivity/{scenario}_substrate_costs.png')

        if metadata['feedback']:
            predicted_data = np.genfromtxt(f'../scenarios/results/Predicted Data {scenario}.csv', delimiter=' ')
            # We create a NaN array to shift our data
            data_shift = np.zeros((predicted_data.shape[0], metadata['t_stp_ahead_pred'] - 1))
            predicted_data = np.hstack((data_shift, predicted_data))
            predicted_data = predicted_data[:, 0:predicted_data.shape[1] - metadata['t_stp_ahead_pred']]
            # We remove the values, so they don't appear on the plot
            predicted_data[:, 0:metadata['t_stp_ahead_pred'] + 1] = np.NaN

            if mpc['mpc'].meta_data['n_robust'] == 0:
                fig.axes[2].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[1, :], 'black',
                                 linestyle='dotted', label=r'$pH_{controller}$')
                fig.axes[1].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[0, :], 'black',
                                 linestyle='dotted', label=r"$\dot V_{g, controller}$")
                fig.axes[1].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[2, :], 'blue',
                                 linestyle='dotted', label=r"$\dot V_{CH_4, controller}$")

            if mpc['mpc'].meta_data['n_robust'] > 0:
                fig.axes[2].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[2, 3], :].transpose(),
                                 'black',
                                 linestyle='dotted', label=r'$pH_{controller}$')
                fig.axes[1].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[0, 1], :].transpose(),
                                 'black',
                                 linestyle='dotted', label=r"$\dot V_{g, controller}$")
                fig.axes[1].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[4, 5], :].transpose(),
                                 'blue',
                                 linestyle='dotted', label=r"$\dot V_{CH_4, controller}$")

            remove_duplicate_labels(fig, 0, legend_location='center right', bbox_to_anchor=(1, 1))
            remove_duplicate_labels(fig, 2, legend_location='center right', bbox_to_anchor=(1, 0.7))
            remove_duplicate_labels(fig, 3, legend_location='center right', bbox_to_anchor=(1, 0.6))
            '''
                    axins_1 = fig.axes[1].inset_axes([0.4, 0.435, 0.25, 0.25],
                        xlim=(5, 10), ylim=(500, 900))
    
                    if 'nominal' in scenario:
    
                        axins_1.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[0, :], 'blue',
                                     linestyle = 'dotted', linewidth=1)
    
                    if 'robust' in scenario:
    
                        axins_1.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[0, 1], :].transpose(),
                                         'blue',
                                         linestyle='dotted', linewidth=1)
            '''
            # axins_1.grid(True, linestyle='dashed')
            fig.axes[0].legend(loc='lower right', bbox_to_anchor=(0.67, 0.45))
            remove_duplicate_labels(fig, 1, legend_location='center right')
            plt.savefig(f'../scenarios/results/plots/Scenario 1/{scenario}_substrate_costs.png')
