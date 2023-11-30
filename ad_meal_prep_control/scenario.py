from dataclasses import dataclass
import pygame
import ad_meal_prep_control.visualization as vis
from mpc import mpc_setup
import numpy as np
import do_mpc
import substrates
from uncertainties import ufloat
from models import adm1_r3_frac_norm
from state_estimator import StateFeedback, mhe_setup
from simulator import simulator_setup
import copy
import matplotlib.pyplot as plt
import params_R3
from utils import StateObserver, ScenarioData
import os
from pathlib import Path

np.random.seed(seed=42)


@dataclass(kw_only=True)
class Scenario:
    scenario_data: ScenarioData

    def __post_init__(self):
        self._n_steps_steady_state = round(
            self.scenario_data.n_days_steady_state / self.scenario_data.t_step
        )
        self._n_steps_mpc = round(
            self.scenario_data.n_days_mpc / self.scenario_data.t_step
        )

        self._t_mpc = np.linspace(
            0,
            self.scenario_data.n_days_mpc,
            num=round(self.scenario_data.n_days_mpc / self.scenario_data.t_step),
            endpoint=True,
        )

        # Take ownership of Tx and x0 because these values are subject to change within the instances of this class
        self.Tx = np.copy(self.scenario_data.Tx)
        self.x0_norm_true = np.copy(self.scenario_data.x0_true) / self.Tx

        # Compute estimated state vector
        self.x0_norm_estimated = np.copy(self.x0_norm_true) * (
            1.0 + 0.1 * np.random.randn(self.x0_norm_true.shape[0])
        )

        self._v_ch4_norm_true_0 = self.x0_norm_true[18]
        self._v_co2_norm_true_0 = self.x0_norm_true[19]
        self._v_ch4_norm_estimated_0 = self.x0_norm_estimated[18]
        self._v_co2_norm_estimated_0 = self.x0_norm_estimated[19]

        # Look for HSL solver
        self._hsllib = None
        for x in os.listdir("./"):
            if x.startswith("coinhsl") and os.path.isdir(x):
                path_to_check = Path(x, "builddir", "libcoinhsl.so")
                if path_to_check.exists():
                    self._hsllib = path_to_check
                    break

    @property
    def x0_norm_true(self) -> np.ndarray:
        return self._x0_norm_true

    @x0_norm_true.setter
    def x0_norm_true(self, new_x0_norm_true: np.ndarray):
        self._x0_norm_true = new_x0_norm_true

    @property
    def x0_norm_estimated(self) -> np.ndarray:
        return self._x0_norm_estimated

    @x0_norm_estimated.setter
    def x0_norm_estimated(self, new_x0_norm_estimated: np.ndarray):
        self._x0_norm_estimated = new_x0_norm_estimated

    @property
    def x0_denorm_true(self) -> np.ndarray:
        return (self._x0_norm_true.T * self.Tx).T

    @property
    def x0_denorm_estimated(self) -> np.ndarray:
        return (self._x0_norm_estimated.T * self.Tx).T

    @property
    def Tx(self) -> np.ndarray:
        return self._Tx

    @Tx.setter
    def Tx(self, new_Tx: np.ndarray):
        self._Tx = new_Tx

    def setup(self):
        self._substrate_setup()
        self._setup_state_and_normalization_vectors()
        self._model_setup()

        if self.scenario_data.pygame_vis:
            self._pygame_setup()

        if self.scenario_data.simulate_steady_state:
            self._sim_setup(ch4_outflow_rate=np.zeros(self._n_steps_steady_state))
        else:
            self._sim_setup(self.scenario_data.ch4_outflow_rate)

        # Estimator setup
        self._estimator_setup()

    def run(self):
        if self.scenario_data.simulate_steady_state:
            self._run_steady_state_sim()

            self._set_new_Tx_x0()
            self._substrate_setup()
            self._model_setup()

            if self.scenario_data.external_gas_storage_model:
                self.x0_norm_true[18] = self._v_ch4_norm_true_0
                self.x0_norm_true[19] = self._v_co2_norm_true_0
                self.x0_norm_estimated[18] = self._v_ch4_norm_estimated_0
                self.x0_norm_estimated[19] = self._v_co2_norm_estimated_0
            self._sim_setup(self.scenario_data.ch4_outflow_rate)
            self._estimator_setup()

        if self.scenario_data.simulate_mpc:
            self._mpc = mpc_setup(
                model=self.model,
                t_step=self.scenario_data.t_step,
                n_horizon=self.scenario_data.controller_params.mpc_n_horizon,
                n_robust=self.scenario_data.controller_params.mpc_n_robust,
                xi_ch_norm=self._xi_ch_mpc_mhe_norm,
                xi_pr_norm=self._xi_pr_mpc_mhe_norm,
                xi_li_norm=self._xi_li_mpc_mhe_norm,
                compile_nlp=self.scenario_data.compile_nlp,
                ch4_outflow_rate=self.scenario_data.ch4_outflow_rate,
                cost_func=self.scenario_data.controller_params.cost_func,
                substrate_costs=[sub.cost for sub in self._subs],
                consider_substrate_costs=self.scenario_data.controller_params.consider_substrate_costs,
                store_full_solution=self.scenario_data.mpc_live_vis,
                bounds=self.scenario_data.controller_params.bounds,
                nl_cons=self.scenario_data.controller_params.nl_cons,
                rterm=self.scenario_data.controller_params.rterm,
                hsllib=self._hsllib,
            )

            if self.scenario_data.mpc_live_vis:
                self._graphics_setup()

            self._run_mpc()

    def _pygame_setup(self):
        pygame.init()
        screen_size = (1280, 720)
        self._screen = pygame.display.set_mode(screen_size)
        pygame.time.Clock()

        self._bga = vis.BioGasPlantVis(params_R3.V_GAS_STORAGE_MAX, self._screen)
        self._data = vis.DataVis(self._screen)

    def _substrate_setup(self):
        # Get the substrate objects
        self._subs = []
        for sub_name in self.scenario_data.sub_names:
            self._subs.append(getattr(substrates, sub_name))

        # Set substrate feeding limits
        self._limited_subs_indices = []
        if isinstance(self.scenario_data.limited_substrates, list):
            for lim_sub in self.scenario_data.limited_substrates:
                for idx, sub in enumerate(self._subs):
                    if sub.name == lim_sub.name:
                        sub.set_limit(
                            lim_sub.amount_remaining,
                            lim_sub.days_remaining,
                            self.scenario_data.t_step,
                        )
                        self._limited_subs_indices.append(idx)

        xi = [sub.xi for sub in self._subs]

        # Normalize xi values (also includes ch, pr and li but these are ignored in the model and later replaced with the uncertain ones)
        self._xi_norm = list(np.array([val / self.Tx[:18] for val in xi]).T)

        # Compute the uncertain xis (ch, pr and li)
        uncertain_xis = [sub.get_uncertain_xi_ch_pr_li() for sub in self._subs]

        # Get nominal values and standard deviations for uncertain xis of all substrates
        xi_ch_nom = np.array([un_xi[0].nominal_value for un_xi in uncertain_xis])
        xi_pr_nom = np.array([un_xi[1].nominal_value for un_xi in uncertain_xis])
        xi_li_nom = np.array([un_xi[2].nominal_value for un_xi in uncertain_xis])

        xi_ch_std_dev = np.array([un_xi[0].std_dev for un_xi in uncertain_xis])
        xi_pr_std_dev = np.array([un_xi[1].std_dev for un_xi in uncertain_xis])
        xi_li_std_dev = np.array([un_xi[2].std_dev for un_xi in uncertain_xis])

        if (
            self.scenario_data.controller_params.num_std_devs == 0.0
            or self.scenario_data.controller_params.mpc_n_robust == 0
        ):
            self._xi_ch_mpc_mhe = np.array([xi_ch_nom])
            self._xi_pr_mpc_mhe = np.array([xi_pr_nom])
            self._xi_li_mpc_mhe = np.array([xi_li_nom])
        else:
            self._xi_ch_mpc_mhe = np.array(
                [
                    xi_ch_nom
                    - self.scenario_data.controller_params.num_std_devs * xi_ch_std_dev,
                    xi_ch_nom
                    + self.scenario_data.controller_params.num_std_devs * xi_ch_std_dev,
                ]
            )
            self._xi_pr_mpc_mhe = np.array(
                [
                    xi_pr_nom
                    - self.scenario_data.controller_params.num_std_devs * xi_pr_std_dev,
                    xi_pr_nom
                    + self.scenario_data.controller_params.num_std_devs * xi_pr_std_dev,
                ]
            )
            self._xi_li_mpc_mhe = np.array(
                [
                    xi_li_nom
                    - self.scenario_data.controller_params.num_std_devs * xi_li_std_dev,
                    xi_li_nom
                    + self.scenario_data.controller_params.num_std_devs * xi_li_std_dev,
                ]
            )

        if self.scenario_data.num_std_devs_sim == 0.0:
            self._x_ch_in_sim = np.array([xi_ch_nom])
            self._x_pr_in_sim = np.array([xi_pr_nom])
            self._x_li_in_sim = np.array([xi_li_nom])
        else:
            self._x_ch_in_sim = np.array(
                [
                    xi_ch_nom
                    + np.random.choice([-1, 1])
                    * self.scenario_data.num_std_devs_sim
                    * xi_ch_std_dev
                ]
            )
            self._x_pr_in_sim = np.array(
                [
                    xi_pr_nom
                    + np.random.choice([-1, 1])
                    * self.scenario_data.num_std_devs_sim
                    * xi_pr_std_dev
                ]
            )
            self._x_li_in_sim = np.array(
                [
                    xi_li_nom
                    + np.random.choice([-1, 1])
                    * self.scenario_data.num_std_devs_sim
                    * xi_li_std_dev
                ]
            )

        # Normalize uncertain xi's
        self._xi_ch_mpc_mhe_norm = self._xi_ch_mpc_mhe / self.Tx[5]
        self._xi_pr_mpc_mhe_norm = self._xi_pr_mpc_mhe / self.Tx[7]
        self._xi_li_mpc_mhe_norm = self._xi_li_mpc_mhe / self.Tx[8]

        self._x_ch_in_sim_norm = self._x_ch_in_sim / self.Tx[5]
        self._x_pr_in_sim_norm = self._x_pr_in_sim / self.Tx[7]
        self._x_li_in_sim_norm = self._x_li_in_sim / self.Tx[8]

    def _estimator_setup(self):
        if self.scenario_data.state_observer == StateObserver.MHE:
            self._estimator = mhe_setup(
                model=self.model,
                t_step=self.scenario_data.t_step,
                n_horizon=self.scenario_data.mhe_n_horizon,
                xi_ch_norm=self._xi_ch_mpc_mhe_norm,
                xi_pr_norm=self._xi_pr_mpc_mhe_norm,
                xi_li_norm=self._xi_li_mpc_mhe_norm,
                P_x=np.diag((self.x0_norm_estimated - self.x0_norm_true) ** 2),
                P_v=0.0001 * np.ones((8, 8)),
                ch4_outflow_rate=self.scenario_data.ch4_outflow_rate,
                store_full_solution=self.scenario_data.mpc_live_vis,
                hsllib=self._hsllib,
            )

            self._estimator.x0 = np.copy(self.x0_norm_estimated)
            self._estimator.set_initial_guess()
        elif self.scenario_data.state_observer == StateObserver.STATEFEEDBACK:
            self._estimator = StateFeedback(self.model, self._simulator)

    def _setup_state_and_normalization_vectors(self):
        """
        Shortens the state and normalization vectors if they've been handed too long if the gas storage had been considered in a previous simulation but not anymore.
        Also normalizes the state vector and sets up the normalization vector for the inputs.
        """
        if not self.scenario_data.external_gas_storage_model:
            self.x0_norm_true = self.x0_norm_true[:18]
            self.x0_norm_estimated = self.x0_norm_estimated[:18]
            self.Tx = self.Tx[:18]

        self.Tu = np.array([self.scenario_data.u_max[sub.state] for sub in self._subs])

    def _model_setup(self):
        # Model
        self.model = adm1_r3_frac_norm(
            xi_norm=self._xi_norm,
            Tu=self.Tu,
            Tx=self.Tx,
            Ty=self.scenario_data.Ty,
            external_gas_storage_model=self.scenario_data.external_gas_storage_model,
            limited_subs_indices=self._limited_subs_indices,
        )

    def _sim_setup(self, ch4_outflow_rate: np.ndarray):
        self._simulator = simulator_setup(
            model=self.model,
            t_step=self.scenario_data.t_step,
            xi_ch_norm=self._x_ch_in_sim_norm,
            xi_pr_norm=self._x_pr_in_sim_norm,
            xi_li_norm=self._x_li_in_sim_norm,
            ch4_outflow_rate=ch4_outflow_rate,
        )

        # Set normalized x0
        self._simulator.x0 = np.copy(self.x0_norm_true)
        self._simulator.set_initial_guess()

    def _graphics_setup(self):
        # Initialize graphic:
        self._graphics = {}
        self._graphics["mpc_graphics"] = do_mpc.graphics.Graphics(self._mpc.data)
        self._graphics["sim_graphics"] = do_mpc.graphics.Graphics(self._simulator.data)
        if self.scenario_data.state_observer == StateObserver.MHE:
            self._graphics["mhe_graphics"] = do_mpc.graphics.Graphics(
                self._estimator.data
            )

        # Configure plot:
        plt.rcParams["axes.grid"] = True
        num_plots = (
            len(self.scenario_data.plot_vars) + 1
            if self.scenario_data.external_gas_storage_model
            else len(self.scenario_data.plot_vars)
        )
        self._fig, self._ax = plt.subplots(num_plots, sharex=True)

        for idx, var in enumerate(self.scenario_data.plot_vars):
            if var[0] == "u":
                self._graphics["mpc_graphics"].add_line(
                    var_type=f"_{var[0]}", var_name=var, axis=self._ax[idx]
                )
                self._ax[idx].legend(
                    labels=[
                        sub.lower().replace("_", " ")
                        for sub in self.scenario_data.sub_names
                    ],
                    title="Substrates",
                    loc="center right",
                )
                self._ax[idx].set_ylabel(var)
            elif var[0] == "x":
                self._graphics["sim_graphics"].add_line(
                    var_type=f"_{var[0]}", var_name=var, axis=self._ax[idx]
                )
                self._ax[idx].legend(
                    self._graphics["sim_graphics"].result_lines[var],
                    labels=["Actual state"],
                    loc="center right",
                )

                self._graphics["mpc_graphics"].add_line(
                    var_type=f"_{var[0]}", var_name=var, axis=self._ax[idx]
                )

                if self.scenario_data.state_observer == StateObserver.MHE:
                    self._graphics["mhe_graphics"].add_line(
                        var_type=f"_{var[0]}", var_name=var, axis=self._ax[idx]
                    )

                    self._ax[idx].legend(
                        labels=[
                            "Actual state",
                            "MPC",
                            "MPC prediction",
                            "MHE estimation",
                        ],
                        loc="center right",
                    )
                else:
                    self._ax[idx].legend(
                        labels=[
                            "Actual state",
                            "MPC",
                            "MPC prediction",
                        ],
                        loc="center right",
                    )
                self._ax[idx].set_ylabel(
                    self.scenario_data._state_names[int(var.split("_")[-1]) - 1]
                )

        self._ax[-1].set_xlabel("Time [d]")

        if self.scenario_data.state_observer == StateObserver.MHE:
            for line_i in self._graphics["mhe_graphics"].result_lines.full:
                line_i.set_alpha(0.4)
                line_i.set_linewidth(6)
        # Update properties for all prediction lines:
        for line_i in self._graphics["mpc_graphics"].pred_lines.full:
            line_i.set_linewidth(1)

        self._fig.align_ylabels()
        # self._fig.tight_layout()
        plt.ion()

    def _run_steady_state_sim(self):
        # Run steady state simulation
        for _ in range(self._n_steps_steady_state):
            self._screen.fill("white")

            u_norm_steady_state = np.array([[0.1 for _ in self._subs]]).T

            y_next = self._simulator.make_step(u_norm_steady_state)
            self.x0_norm_true = np.array(self._simulator.x0.master)
            self.x0_norm_estimated = self._estimator.make_step(y_next)

        # Save normalized x values
        np.savetxt(
            f"./results/{self.scenario_data.name}_steady_state_x.csv",
            self._simulator.data._x,
            delimiter=",",
        )

    def _run_mpc(self):
        # Initialize simulator and controller
        self._mpc.x0 = np.copy(self.x0_norm_estimated)
        self._mpc.u0 = np.ones(len(self._subs)) * 0.5
        self._mpc.set_initial_guess()

        u_norm_computed = None

        # MPC
        for k in range(self._n_steps_mpc):
            # fill the screen with a color to wipe away anything from last frame
            self._screen.fill("white")

            u_norm_computed_old = copy.copy(u_norm_computed)
            u_norm_computed = self._mpc.make_step(self.x0_norm_estimated)

            u_norm_actual = np.copy(u_norm_computed)

            # Manipulate the actual feed to the biogas plant 'u_norm_actual'
            # based on the set disturbances
            if not self.scenario_data.disturbances.feed_computation_stuck is None:
                stuck_start_idx = (
                    self.scenario_data.disturbances.feed_computation_stuck[0]
                )
                stuck_end_idx = (
                    stuck_start_idx
                    + self.scenario_data.disturbances.feed_computation_stuck[1]
                )
                if k in range(stuck_start_idx, stuck_end_idx):
                    u_norm_actual = copy.copy(u_norm_computed_old)
                    u_norm_computed = copy.copy(u_norm_computed_old)

            if not self.scenario_data.disturbances.clogged_feeding is None:
                for (
                    sub_idx,
                    val,
                ) in self.scenario_data.disturbances.clogged_feeding.items():
                    if k in range(val[0], val[0] + val[1]):
                        u_norm_actual[sub_idx, 0] = 0.0

            if not self.scenario_data.disturbances.max_feeding_error is None:
                u_norm_actual *= (
                    1.0
                    + np.array(
                        [
                            np.random.uniform(
                                low=-1.0, high=1.0, size=u_norm_actual.shape[0]
                            )
                            * self.scenario_data.disturbances.max_feeding_error
                        ]
                    ).T
                )

            y_next = self._simulator.make_step(u_norm_actual)
            self.x0_norm_true = np.array(self._simulator.x0.master)
            self.x0_norm_estimated = self._estimator.make_step(y_next)

            if not self.scenario_data.disturbances.state_jumps is None:
                for x_idx, val in self.scenario_data.disturbances.state_jumps.items():
                    if k == val[0]:
                        self.x0_norm_estimated[x_idx] += val[1]

            if self.scenario_data.mpc_live_vis:
                for g in self._graphics.values():
                    g.plot_results(t_ind=k)
                self._graphics["mpc_graphics"].plot_predictions(t_ind=k)
                self._graphics["mpc_graphics"].reset_axes()
                self._graphics["sim_graphics"].reset_axes()
                for idx, var in enumerate(self.scenario_data.plot_vars):
                    if var[0] == "y":
                        y_num = int(var.split("_")[-1])

                        self._ax[idx].scatter(
                            self._t_mpc[: k + 1],
                            self._simulator.data._y[:, self.model.u.size + y_num - 1],
                            color="red",
                        )
                        self._ax[idx].set_ylabel(
                            self.scenario_data._meas_names[int(var.split("_")[-1]) - 1]
                        )

                if self.scenario_data.external_gas_storage_model:
                    self._ax[-1].scatter(
                        self._t_mpc[: k + 1],
                        self.scenario_data.ch4_outflow_rate[: k + 1],
                        color="red",
                    )
                    self._ax[-1].set_ylabel("CH4 outflow\nvolume flow")

                plt.show()
                plt.pause(0.01)

            if self.scenario_data.pygame_vis:
                vis.visualize(
                    self.scenario_data,
                    self._bga,
                    self._data,
                    self.x0_norm_true,
                    self.Tx,
                    self._simulator,
                    u_norm_actual,
                )

        if self.scenario_data.save_results:
            self._save_results()

    def _save_results(self):
        do_mpc.data.save_results(
            save_list=[self._mpc, self._simulator],
            result_name=f"{self.scenario_data.name}_mpc_results",
            result_path="./results/",
            overwrite=True,
        )

    def _set_new_Tx_x0(self):
        """
        Sets the new Tx normalization vector based on the simulation results from the initial steady state simulation as well as the new x0 value.
        """
        # Get denormalized steady state state-vector
        x_steady_state_true = np.copy(self._simulator.x0.master * self.Tx)
        x_steady_state_estimated = np.copy(self.x0_denorm_estimated)

        # Set new Tx based on max absolute values of states during steady state simulation
        self.Tx *= np.max(np.abs(self._simulator.data._x), axis=0)

        # Set new normalized initial state (for MPC)
        self.x0_norm_true = (x_steady_state_true.T / self.Tx).T
        self.x0_norm_estimated = (x_steady_state_estimated.T / self.Tx).T
