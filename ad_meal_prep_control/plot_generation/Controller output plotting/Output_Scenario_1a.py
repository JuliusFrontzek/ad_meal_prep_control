import numpy as np
from do_mpc.data import save_results, load_results
import do_mpc
import matplotlib.pyplot as plt
import pickle as pkl
from ad_meal_prep_control.utils import remove_duplicate_labels

scenario_names = ["Scenario_1a_quadratic_no_feedback_mismatch_1std_ch",
                  'Scenario_1a_quadratic_no_feedback_mismatch_1std_pr',
                  'Scenario_1a_quadratic_no_feedback_mismatch_1std_li']

#scenario_names = ["Scenario_1a_quadratic_feedback_mismatch_5std_3tsap",
#                  'Scenario_1a_quadratic_nominal_ideal_feedback_3tsap',
#                  'Scenario_1a_quadratic_robust_feedback_mismatch_5std_3tsap']

for scenario in scenario_names:
    with open(f'./results/plots/{scenario}_substrate_costs.pkl', 'rb') as file:
        fig = pkl.load(file)
    mpc = load_results(f'./results/{scenario}_mpc_results.pkl')
    metadata = load_results(f'./results/{scenario}_scenario_meta_data.pkl')
    # Graph
    if not metadata['feedback']:
        plant_output = np.genfromtxt(f'./results/Plant Output {scenario}.csv', delimiter=' ')

        fig.axes[1].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 1], 'rebeccapurple',
                         linestyle = 'dotted', label=r'$pH_{controller}$')
        fig.axes[2].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 0], 'blue',
                         linestyle = 'dotted', label=r"$\dot V_{g, controller}$")
        fig.axes[2].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 2], 'rebeccapurple',
                         linestyle = 'dotted', label=r"$\dot V_{CH_4, controller}$")

        remove_duplicate_labels(fig, 0)
        remove_duplicate_labels(fig, 1)
        remove_duplicate_labels(fig, 2, bbox_to_anchor=(1, 1))
        fig.axes[0].legend(loc="lower right")
        plt.savefig(f'./results/plots/Plant Output {scenario}.png')
        plt.show()

    if metadata['feedback']:
        predicted_data = np.genfromtxt(f'./results/Predicted Data {scenario}.csv', delimiter=' ')
        # We create a NaN array to shift our data
        data_shift = np.zeros((predicted_data.shape[0], metadata['t_stp_ahead_pred']-1))
        predicted_data = np.hstack((data_shift, predicted_data))
        predicted_data = predicted_data[:, 0:predicted_data.shape[1] - metadata['t_stp_ahead_pred']]
        # We remove the values, so they don't appear on the plot
        predicted_data[:, 0:metadata['t_stp_ahead_pred']+1] = np.NaN

        if mpc['mpc'].meta_data['n_robust'] == 0:
            fig.axes[1].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[1, :], 'rebeccapurple',
                         linestyle = 'dotted', label=r'$pH_{controller}$')
            fig.axes[2].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[0, :], 'blue',
                         linestyle = 'dotted', label=r"$\dot V_{g, controller}$")
            fig.axes[2].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[2, :], 'rebeccapurple',
                         linestyle = 'dotted', label=r"$\dot V_{CH_4, controller}$")

        if mpc['mpc'].meta_data['n_robust'] > 0:
            fig.axes[1].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[2, 3], :].transpose(),
                             'rebeccapurple',
                             linestyle='dotted', label=r'$pH_{controller}$')
            fig.axes[2].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[0, 1], :].transpose(),
                             'blue',
                             linestyle='dotted', label=r"$\dot V_{g, controller}$")
            fig.axes[2].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[4, 5], :].transpose(),
                             'rebeccapurple',
                             linestyle='dotted', label=r"$\dot V_{CH_4, controller}$")

        remove_duplicate_labels(fig, 0)
        remove_duplicate_labels(fig, 1)
        remove_duplicate_labels(fig, 2)

        axins_1 = fig.axes[2].inset_axes([0.4, 0.435, 0.25, 0.25],
            xlim=(5.5, 8), ylim=(1000, 1300))

        if 'nominal' in scenario:

            axins_1.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[0, :], 'blue',
                         linestyle = 'dotted', linewidth=1)

        if 'robust' in scenario:

            axins_1.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[0, 1], :].transpose(),
                             'blue',
                             linestyle='dotted', linewidth=1)

        fig.axes[0].legend(loc = 'lower right')
        axins_1.grid(True, linestyle='dashed')
        remove_duplicate_labels(fig, 2, bbox_to_anchor=(1, 1))
        plt.savefig(f'./results/plots/predicted_data {scenario}.png')
        plt.show()