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

scenario_names = ["Scenario_2a_dynamic_nominal_feedback_mismatch_3std_8tsap",
                  'Scenario_2a_dynamic_nominal_ideal_feedback_8tsap',
                  'Scenario_2a_dynamic_robust_feedback_mismatch_3std_8tsap']

for scenario in scenario_names:
    mpc = load_results(f'./results/{scenario}_mpc_results.pkl')
    metadata = load_results(f'./results/{scenario}_scenario_meta_data.pkl')
    with open(f'./results/plots/{scenario}.pkl', 'rb') as file:
        fig = pkl.load(file)

    # Graph

    if not metadata['feedback']:
        plant_output = np.genfromtxt(f'./results/Plant Output {scenario}.csv', delimiter=' ')
        plant_output = plant_output[:-1,:]
        plant_output[0, :] = np.NaN

        fig.axes[1].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 0]*100, 'limegreen',
                         linestyle = 'dotted', label=r"$V_{CH_4,tank, controller}$", linewidth=1)
        fig.axes[1].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 1]*100, 'darkgreen',
                         linestyle = 'dotted', label=r"$V_{CO_2,tank, controller}$", linewidth=1)
        fig.axes[2].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 2], 'rebeccapurple',
                         linestyle = 'dotted', label=r"$V_{g, tank, controller}$", linewidth=1)
        fig.axes[3].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 3], 'rebeccapurple',
                         linestyle = 'dotted', label=r"$\dot V_{CH_4,AD, controller}$", linewidth=1)
        fig.axes[3].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 4], 'blue',
                         linestyle = 'dotted', label=r"$\dot V_{g, AD, controller}$", linewidth=1)
        fig.axes[4].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 5], 'rebeccapurple',
                         linestyle = 'dotted', label=r'$pH_{model}$', linewidth=1)
        fig.axes[0].legend()
        fig.axes[1].legend(bbox_to_anchor=(0.95, 1))
        fig.axes[2].legend()
        fig.axes[3].legend(bbox_to_anchor=(0.95, 1))
        fig.axes[4].legend()
        plt.savefig(f'./results/plots/Plant Output {scenario}.png')
    if metadata['feedback']:
        predicted_data = np.genfromtxt(f'./results/Predicted Data {scenario}.csv', delimiter=' ')
        # We create a NaN array to shift our data
        data_shift = np.zeros((predicted_data.shape[0], metadata['t_stp_ahead_pred']-1))
        predicted_data = np.hstack((data_shift, predicted_data))
        predicted_data = predicted_data[:, 0:predicted_data.shape[1] - metadata['t_stp_ahead_pred']]
        # We remove the value, so they don't appear on the plot
        predicted_data[:, 0:metadata['t_stp_ahead_pred'] + 1] = np.NaN

        if mpc['mpc'].meta_data['n_robust'] == 0:
            fig.axes[1].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[0, :]*100, 'limegreen',
                             linestyle = 'dotted', label=r"$V_{CH_4,tank, controller}$", linewidth=1)
            fig.axes[1].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[1, :]*100, 'darkgreen',
                             linestyle = 'dotted', label=r"$V_{CO_2,tank, controller}$", linewidth=1)
            fig.axes[2].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[2, :], 'rebeccapurple',
                             linestyle = 'dotted', label=r"$V_{g, tank, controller}$", linewidth=1)
            fig.axes[3].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[3, :], 'rebeccapurple',
                             linestyle = 'dotted', label=r"$\dot V_{CH_4,AD, controller}$", linewidth=1)
            fig.axes[3].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[4, :], 'blue',
                             linestyle = 'dotted', label=r"$\dot V_{g, AD, controller}$", linewidth=1)
            fig.axes[4].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[5, :], 'rebeccapurple',
                             linestyle = 'dotted', label=r'$pH_{controller}$', linewidth=1)

        if mpc['mpc'].meta_data['n_robust'] > 0:
            fig.axes[1].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[0,1], :].transpose()*100, 'limegreen',
                             linestyle = 'dotted', label=r"$V_{CH_4,tank, controller}$", linewidth=1)
            fig.axes[1].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[2,3], :].transpose()*100, 'darkgreen',
                             linestyle = 'dotted', label=r"$V_{CO_2,tank, controller}$", linewidth=1)
            fig.axes[2].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[4,5], :].transpose(), 'rebeccapurple',
                             linestyle = 'dotted', label=r"$V_{g, tank, controller}$", linewidth=1)
            fig.axes[3].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[6,7], :].transpose(), 'rebeccapurple',
                             linestyle = 'dotted', label=r"$\dot V_{CH_4,AD, controller}$", linewidth=1)
            fig.axes[3].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[8,9], :].transpose(), 'blue',
                             linestyle = 'dotted', label=r"$\dot V_{g, AD, controller}$", linewidth=1)
            fig.axes[4].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[10,11], :].transpose(), 'rebeccapurple',
                             linestyle = 'dotted', label=r'$pH_{controller}$', linewidth=1)

        remove_duplicate_labels(fig, 0, bbox_to_anchor=(0.35, 0.35))
        remove_duplicate_labels(fig, 1, bbox_to_anchor=(0.945, 1))
        remove_duplicate_labels(fig, 2, bbox_to_anchor=(0.95, 1))
        remove_duplicate_labels(fig, 3, bbox_to_anchor=(0.95, 1))
        remove_duplicate_labels(fig, 4)

        axins_1 = fig.axes[1].inset_axes([0.2, 0.15, 0.25, 0.25],
            xlim=(17, 17.5), ylim=(15, 35))
        axins_2 = fig.axes[2].inset_axes([0.2, 0.15, 0.25, 0.25],
                                         xlim=(17.15, 17.45), ylim=(180, 250))
        axins_3 = fig.axes[3].inset_axes([0.2, 0.15, 0.25, 0.25],
                                         xlim=(21.5, 22.5), ylim=(300, 400))
        #axins_3 = fig.axes[3].inset_axes([0.2, 0.15, 0.25, 0.25],
        #                                 xlim=(17.1, 17.4), ylim=(350, 500))
        if 'nominal' in scenario:

            axins_1.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[0, :]*100, 'limegreen',
                             linestyle = 'dotted', linewidth=1)
            axins_1.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[1, :]*100, 'darkgreen',
                             linestyle = 'dotted', linewidth=1)
            axins_2.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[2, :], 'rebeccapurple',
                             linestyle = 'dotted', linewidth=1)
            axins_3.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[3, :], 'rebeccapurple',
                             linestyle = 'dotted', linewidth=1)
            axins_3.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[4, :], 'blue',
                             linestyle = 'dotted', linewidth=1)

        if 'robust' in scenario:

            axins_1.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[0,1], :].transpose()*100, 'limegreen',
                             linestyle = 'dotted', linewidth=1)
            axins_1.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[2,3], :].transpose()*100, 'darkgreen',
                             linestyle = 'dotted', linewidth=1)
            axins_2.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[4,5], :].transpose(), 'rebeccapurple',
                             linestyle = 'dotted', linewidth=1)
            axins_3.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[6,7], :].transpose(), 'rebeccapurple',
                             linestyle = 'dotted', linewidth=1)
            axins_3.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[8,9], :].transpose(), 'blue',
                             linestyle = 'dotted', linewidth=1)

        axins_2.plot(np.linspace(0, 30, num=predicted_data.shape[1]),
                     np.ones((predicted_data.shape[1],1))*V_GAS_STORAGE_MAX, 'red',
                     linestyle='dashed')
        axins_2.plot(np.linspace(0, 30, num=predicted_data.shape[1]),
                     np.ones((predicted_data.shape[1],1))*V_GAS_STORAGE_MAX*(1-metadata['controller_params']['gas_storage_bound_fraction']), 'blue',
                     linestyle='dashed')
        axins_1.grid(True, linestyle='dashed')
        axins_2.grid(True, linestyle='dashed')
        axins_3.grid(True, linestyle='dashed')

    plt.savefig(f'./results/plots/predicted_data scenario_2a_dynamic_{scenario}.png')
    plt.show()