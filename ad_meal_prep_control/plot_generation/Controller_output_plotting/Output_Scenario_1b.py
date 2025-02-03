import numpy as np
from do_mpc.data import save_results, load_results
import do_mpc
import matplotlib.pyplot as plt
import pickle as pkl
from ad_meal_prep_control.utils import remove_duplicate_labels


def controller_plotting_1b(scenario_names=None):
    if scenario_names is None:
        scenario_names = ["Scenario_1b_quadratic_nominal_feedback_mismatch_2std_3tsap"]  # default name

    for scenario in scenario_names:
        with open(f'../scenarios/results/plots/{scenario}.pkl', 'rb') as file:
            fig = pkl.load(file)
            fig.set_dpi(1000)
        mpc = load_results(f'../scenarios/results/{scenario}_mpc_results.pkl')
        metadata = load_results(f'../scenarios/results/{scenario}_scenario_meta_data.pkl')
        # Graph
        if not metadata['feedback']:
            plant_output = np.genfromtxt(f'../scenarios/results/{scenario}.csv', delimiter=' ')

            fig.axes[3].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 1], 'rebeccapurple',
                             linestyle='dotted', label=r'$pH_{controller}$')
            fig.axes[2].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 0], 'blue',
                             linestyle='dotted', label=r"$\dot V_{g, controller}$")
            fig.axes[2].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 2], 'rebeccapurple',
                             linestyle='dotted', label=r"$\dot V_{CH_4, controller}$")

            remove_duplicate_labels(fig, 0)
            remove_duplicate_labels(fig, 2)
            remove_duplicate_labels(fig, 3)
            remove_duplicate_labels(fig, 1, bbox_to_anchor=(1, 1))
            fig.axes[0].legend(loc="lower right")
            plt.savefig(f'../scenarios/results/plots/Plant Output {scenario}.png')
            plt.show()

        if metadata['feedback']:
            predicted_data = np.genfromtxt(f'../scenarios/results/Predicted Data {scenario}.csv', delimiter=' ')
            # We create a NaN array to shift our data
            data_shift = np.zeros((predicted_data.shape[0], metadata['t_stp_ahead_pred'] - 1))
            predicted_data = np.hstack((data_shift, predicted_data))
            predicted_data = predicted_data[:, 0:predicted_data.shape[1] - metadata['t_stp_ahead_pred']]
            # We remove the values, so they don't appear on the plot
            predicted_data[:, 0:metadata['t_stp_ahead_pred'] + 1] = np.NaN

            if mpc['mpc'].meta_data['n_robust'] == 0:
                fig.axes[3].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[1, :], 'black',
                                 linestyle='dotted', label=r'$pH_{controller}$')
                fig.axes[2].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[0, :], 'black',
                                 linestyle='dotted', label=r"$\dot V_{g, controller}$")
                fig.axes[2].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[2, :], 'blue',
                                 linestyle='dotted', label=r"$\dot V_{CH_4, controller}$")

            if mpc['mpc'].meta_data['n_robust'] > 0:
                fig.axes[3].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[2, 3], :].transpose(),
                                 'black',
                                 linestyle='dotted', label=r'$pH_{controller}$')
                fig.axes[2].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[0, 1], :].transpose(),
                                 'black',
                                 linestyle='dotted', label=r"$\dot V_{g, controller}$")
                fig.axes[2].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[4, 5], :].transpose(),
                                 'blue',
                                 linestyle='dotted', label=r"$\dot V_{CH_4, controller}$")

            remove_duplicate_labels(fig, 0, legend_location='upper left', bbox_to_anchor=(0.1, 1))
            remove_duplicate_labels(fig, 1, legend_location='center right')
            remove_duplicate_labels(fig, 2, legend_location='upper right')
            remove_duplicate_labels(fig, 3, legend_location='lower right')
            remove_duplicate_labels(fig, 4, legend_location='center right', bbox_to_anchor=(1, 0.6))
            # fig.axes[0].legend(loc='lower right', bbox_to_anchor=(0.62, 0.45))
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
            plt.savefig(f'../scenarios/results/plots/Scenario 1/{scenario}.png')
    return
