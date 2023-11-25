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
from utils import ScenarioType, StateObserver, ScenarioData


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

        self._sim_setup()

        # Estimator setup
        self._estimator_setup()

    def run(self):
        if self.scenario_data.simulate_steady_state:
            self._run_steady_state_sim()

            self._set_new_Tx_x0()
            self._substrate_setup()
            self._model_setup()
            self._sim_setup()
            self._estimator_setup()

        if self.scenario_data.simulate_mpc:
            self._mpc = mpc_setup(
                model=self.model,
                t_step=self.scenario_data.t_step,
                n_horizon=self.scenario_data.mpc_n_horizon,
                n_robust=self.scenario_data.mpc_n_robust,
                x_ch_in=self._x_ch_in,
                x_pr_in=self._x_pr_in,
                x_li_in=self._x_li_in,
                compile_nlp=self.scenario_data.compile_nlp,
                vol_flow_rate=self.scenario_data.vol_flow_rate,
                cost_func=self.scenario_data.cost_func,
                substrate_costs=[sub.cost for sub in self._subs],
                consider_substrate_costs=self.scenario_data.consider_substrate_costs,
                bounds=self.scenario_data.bounds,
                nl_cons=self.scenario_data.nl_cons,
                rterm=self.scenario_data.rterm,
            )

            if self.scenario_data.mpc_live_vis:
                self._graphics_setup()

            self._run_mpc()

    def _pygame_setup(self):
        pygame.init()
        screen_size = (1280, 720)
        self._screen = pygame.display.set_mode(screen_size)
        pygame.time.Clock()

        self._bga = vis.BioGasPlantVis(150.0, self._screen)
        self._data = vis.DataVis(self._screen)

    def _substrate_setup(self):
        self._subs = []
        for sub_name in self.scenario_data.sub_names:
            self._subs.append(getattr(substrates, sub_name))
        xi = [sub.xi for sub in self._subs]

        self._xi_norm = list(np.array([val / self.Tx[:18] for val in xi]).T)

        uncertain_xis = [sub.get_uncertain_xi_ch_pr_li() for sub in self._subs]

        if not self.scenario_data.consider_uncertainty:
            uncertain_xis[0] = (
                ufloat(xi[0][5], 0.0),
                ufloat(xi[0][7], 0.0),
                ufloat(xi[0][8], 0.0),
            )

        x_ch_nom = np.array([un_xi[0].nominal_value for un_xi in uncertain_xis])
        x_pr_nom = np.array([un_xi[1].nominal_value for un_xi in uncertain_xis])
        x_li_nom = np.array([un_xi[2].nominal_value for un_xi in uncertain_xis])

        x_ch_std_dev = np.array([un_xi[0].std_dev for un_xi in uncertain_xis])
        x_pr_std_dev = np.array([un_xi[1].std_dev for un_xi in uncertain_xis])
        x_li_std_dev = np.array([un_xi[2].std_dev for un_xi in uncertain_xis])

        if self.scenario_data.consider_uncertainty:
            self._x_ch_in = np.array(
                [
                    x_ch_nom - self.scenario_data.num_std_devs * x_ch_std_dev,
                    x_ch_nom + self.scenario_data.num_std_devs * x_ch_std_dev,
                ]
            )
            self._x_pr_in = np.array(
                [
                    x_pr_nom - self.scenario_data.num_std_devs * x_pr_std_dev,
                    x_pr_nom + self.scenario_data.num_std_devs * x_pr_std_dev,
                ]
            )
            self._x_li_in = np.array(
                [
                    x_li_nom - self.scenario_data.num_std_devs * x_li_std_dev,
                    x_li_nom + self.scenario_data.num_std_devs * x_li_std_dev,
                ]
            )
        else:
            self._x_ch_in = np.array([x_ch_nom])
            self._x_pr_in = np.array([x_pr_nom])
            self._x_li_in = np.array([x_li_nom])

        # Normalize uncertain xi's
        self._x_ch_in /= self.Tx[5]
        self._x_pr_in /= self.Tx[7]
        self._x_li_in /= self.Tx[8]

    def _estimator_setup(self):
        if self.scenario_data.state_observer == StateObserver.MHE:
            self._estimator = mhe_setup(
                model=self.model,
                t_step=self.scenario_data.t_step,
                n_horizon=self.scenario_data.mhe_n_horizon,
                x_ch_in=self._x_ch_in,
                x_pr_in=self._x_pr_in,
                x_li_in=self._x_li_in,
                P_x=np.diag((self.x0_norm_estimated - self.x0_norm_true) ** 2),
                P_v=0.0001 * np.ones((8, 8)),
                vol_flow_rate=self.scenario_data.vol_flow_rate,
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
        if self.scenario_data.scenario_type == ScenarioType.METHANATION:
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
            scenario_type=self.scenario_data.scenario_type,
        )

    def _sim_setup(self):
        self._simulator = simulator_setup(
            model=self.model,
            t_step=self.scenario_data.t_step,
            x_ch_in=self._x_ch_in,
            x_pr_in=self._x_pr_in,
            x_li_in=self._x_li_in,
            vol_flow_rate=self.scenario_data.vol_flow_rate,
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
        self._fig, self._ax = plt.subplots(
            len(self.scenario_data.plot_vars), sharex=True
        )

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
            # np.array(
            #     [[8.0 / self.Tu[0] for _ in range(len(self._subs))]]
            # ).T  # 1.0 / 1e4

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
                    params_R3.V_GAS_STORAGE_MAX,
                )

        if self.scenario_data.store_results:
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
