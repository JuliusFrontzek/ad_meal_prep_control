import do_mpc
import matplotlib.pyplot as plt

plt.rcParams["axes.grid"] = True

data = do_mpc.data.load_results("./results/methanation_mpc_results.pkl")

fig, ax, graphics = do_mpc.graphics.default_plot(
    data=data["mpc"], states_list=[], aux_list=["v_ch4_dot_tank_in", "y_4"]
)

ax[0].legend(labels=[i for i in range(4)])

ax[-1].set_xlabel("Time [d]")

plt.show()
