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


def controller_plotting_2c(scenario_names=None):
    if scenario_names is None:
        scenario_names = ["Scenario_2c_dynamic"]
    for scenario in scenario_names:
        mpc = load_results(f'../scenarios/results/{scenario}_mpc_results.pkl')
        metadata = load_results(f'../scenarios/results/{scenario}_scenario_meta_data.pkl')
        with open(f'../scenarios/results/plots/{scenario}.pkl', 'rb') as file:
            fig = pkl.load(file)

        # Graph

        if not metadata['feedback']:
            plant_output = np.genfromtxt(f'../scenarios/results/Plant Output {scenario}.csv', delimiter=' ')
            plant_output = plant_output[:-1, :]
            plant_output[0, :] = np.NaN

            fig.axes[1].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 0] * 100, 'limegreen',
                             linestyle='dotted', label=r"$V_{CH_4,tank, controller}$", linewidth=1)
            fig.axes[1].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 1] * 100, 'darkgreen',
                             linestyle='dotted', label=r"$V_{CO_2,tank, controller}$", linewidth=1)
            fig.axes[2].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 2], 'rebeccapurple',
                             linestyle='dotted', label=r"$V_{g, tank, controller}$", linewidth=1)
            fig.axes[3].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 3], 'rebeccapurple',
                             linestyle='dotted', label=r"$\dot V_{CH_4,AD, controller}$", linewidth=1)
            fig.axes[3].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 4], 'blue',
                             linestyle='dotted', label=r"$\dot V_{g, AD, controller}$", linewidth=1)
            fig.axes[4].plot(np.linspace(0, 30, num=plant_output.shape[0]), plant_output[:, 5], 'rebeccapurple',
                             linestyle='dotted', label=r'$pH_{model}$', linewidth=1)
            fig.axes[0].legend()
            fig.axes[1].legend()
            fig.axes[2].legend()
            fig.axes[3].legend()
            fig.axes[4].legend()
            plt.savefig(f'../scenarios/results/plots/Plant Output {scenario}.png')
        else:
            pass
        '''
            # We create a NaN array to shift our data
            try:
                predicted_data = np.genfromtxt(f'../scenarios/results/Predicted Data {scenario}.csv', delimiter=' ')
                data_shift = np.zeros((predicted_data.shape[0], metadata['t_stp_ahead_pred']-1))
                predicted_data = np.hstack((data_shift, predicted_data))
                predicted_data = predicted_data[:, 0:predicted_data.shape[1] - metadata['t_stp_ahead_pred']]
                # We remove the value, so they don't appear on the plot
                predicted_data[:, 0:metadata['t_stp_ahead_pred'] + 1] = np.NaN

                if mpc['mpc'].meta_data['n_robust'] == 0:
                    fig.axes[1].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[0, :]*100, 'black',
                                     linestyle = 'dotted', label=r"$V_{CH_4,tank, controller}$", linewidth=1)
                    fig.axes[1].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[1, :]*100, 'blue',
                                     linestyle = 'dotted', label=r"$V_{CO_2,tank, controller}$", linewidth=1)
                    fig.axes[2].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[2, :], 'black',
                                     linestyle = 'dotted', label=r"$V_{g, tank, controller}$", linewidth=1)
                    fig.axes[3].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[3, :], 'blue',
                                     linestyle = 'dotted', label=r"$\dot V_{CH_4,AD, controller}$", linewidth=1)
                    fig.axes[3].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[4, :], 'black',
                                     linestyle = 'dotted', label=r"$\dot V_{g, AD, controller}$", linewidth=1)
                    fig.axes[4].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[5, :], 'black',
                                     linestyle = 'dotted', label=r'$pH_{controller}$', linewidth=1)

                if mpc['mpc'].meta_data['n_robust'] > 0:
                    fig.axes[1].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[0,1], :].transpose()*100, 'black',
                                     linestyle = 'dotted', label=r"$V_{CH_4,tank, controller}$", linewidth=1)
                    fig.axes[1].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[2,3], :].transpose()*100, 'blue',
                                     linestyle = 'dotted', label=r"$V_{CO_2,tank, controller}$", linewidth=1)
                    fig.axes[2].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[4,5], :].transpose(), 'black',
                                     linestyle = 'dotted', label=r"$V_{g, tank, controller}$", linewidth=1)
                    fig.axes[3].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[6,7], :].transpose(), 'blue',
                                     linestyle = 'dotted', label=r"$\dot V_{CH_4,AD, controller}$", linewidth=1)
                    fig.axes[3].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[8,9], :].transpose(), 'black',
                                     linestyle = 'dotted', label=r"$\dot V_{g, AD, controller}$", linewidth=1)
                    fig.axes[4].plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[10,11], :].transpose(), 'black',
                                     linestyle = 'dotted', label=r'$pH_{controller}$', linewidth=1)

                    #axins_1 = fig.axes[1].inset_axes([0.4, 0.15, 0.25, 0.25],
                    #                                 xlim=(3, 5), ylim=(10, 50))
                    #axins_2 = fig.axes[2].inset_axes([0.4, 0.15, 0.25, 0.25],
                    #                                 xlim=(3, 5), ylim=(150, 300))
                    #axins_3 = fig.axes[3].inset_axes([0.4, 0.15, 0.25, 0.25],
                    #                                 xlim=(3, 5), ylim=(100, 700))

                if 'nominal_feedback' in scenario:
                    axins_1 = fig.axes[1].inset_axes([0.4, 0.6, 0.35, 0.35],
                                                     xlim=(18, 23), ylim=(-10, 120))
                    axins_2 = fig.axes[2].inset_axes([0.4, 0.6, 0.35, 0.35],
                                                     xlim=(18, 23), ylim=(20, 300))

                    axins_1.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[0, :]*100, 'black',
                                     linestyle = 'dotted', linewidth=1)
                    axins_1.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[1, :]*100, 'blue',
                                     linestyle = 'dotted', linewidth=1)
                    axins_2.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[2, :], 'black',
                                     linestyle = 'dotted', linewidth=1)
                    #axins_3.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[3, :], 'rebeccapurple',
                    #                 linestyle = 'dotted', linewidth=1)
                    #axins_3.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[4, :], 'blue',
                    #                 linestyle = 'dotted', linewidth=1)
                if 'robust' in scenario:
                    pass
                    #axins_1.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[0,1], :].transpose()*100, 'limegreen',
                    #                 linestyle = 'dotted', linewidth=1)
                    #axins_1.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[2,3], :].transpose()*100, 'darkgreen',
                    #                 linestyle = 'dotted', linewidth=1)
                    #axins_2.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[4,5], :].transpose(), 'rebeccapurple',
                    #                 linestyle = 'dotted', linewidth=1)
                    #axins_3.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[6,7], :].transpose(), 'rebeccapurple',
                    #                 linestyle = 'dotted', linewidth=1)
                    #axins_3.plot(np.linspace(0, 30, num=predicted_data.shape[1]), predicted_data[[8,9], :].transpose(), 'blue',
                    #                 linestyle = 'dotted', linewidth=1)

                axins_1.plot(np.linspace(0, 30, num=predicted_data.shape[1]),
                             np.ones((predicted_data.shape[1],1))*0, 'black',
                             linestyle='dashed')
                axins_2.plot(np.linspace(0, 30, num=predicted_data.shape[1]),
                             np.ones((predicted_data.shape[1],1))*V_GAS_STORAGE_MAX, 'black',
                             linestyle='dashed')
                axins_2.plot(np.linspace(0, 30, num=predicted_data.shape[1]),
                             np.ones((predicted_data.shape[1],1))*V_GAS_STORAGE_MAX*(1-metadata['controller_params']['gas_storage_bound_fraction']), 'grey',
                             linestyle='dashed')
                axins_1.grid(True, linestyle='dashed')
                axins_2.grid(True, linestyle='dashed')
                axins_2.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
                axins_1.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
                #axins_3.grid(True, linestyle='dashed')
        
            except:
                pass
        '''
        # remove_duplicate_labels(fig, 0, legend_location='center', bbox_to_anchor=(0.5, 0.5))
        remove_duplicate_labels(fig, 1, legend_location='center right')  # , bbox_to_anchor=(1, 1))
        remove_duplicate_labels(fig, 2, legend_location='center right')  # , bbox_to_anchor=(0, 0.6))
        remove_duplicate_labels(fig, 3, legend_location='center right')  # , bbox_to_anchor=(0, 0.6))
        remove_duplicate_labels(fig, 4, legend_location='center right')  # , bbox_to_anchor=(0, 0.6))
        remove_duplicate_labels(fig, 5, legend_location='center right', bbox_to_anchor=(1, 0.6))
        # fig.axes[0].legend(loc='lower right', bbox_to_anchor=(0.62, 0.45))
        # remove_duplicate_labels(fig, 6, legend_location='center left', bbox_to_anchor=(0, 0.6))
        plt.savefig(f'../scenarios/results/plots/Scenario 2/{scenario}.png')
    return
