from dataclasses import dataclass
import pygame
import ad_meal_prep_control.visualization as vis
from mpc import mpc_setup
import numpy as np
import do_mpc
import substrates
from uncertainties import ufloat
import ad_meal_prep_control.utils as utils
from models import adm1_r3_frac_norm
import state_estimator
from simulator import simulator_setup
from typing import Union
import copy
import matplotlib.pyplot as plt
import params_R3
from utils import ScenarioType


@dataclass(kw_only=True)
class Scenario:
    name: str
    scenario_type: ScenarioType
    n_horizon: int
    n_robust: int
    t_step: float  # Time in days
    n_days_steady_state: float
    n_days_mpc: float
    sub_names: list[str]
    disturbances: utils.Disturbances
    x0: np.ndarray
    Tx: np.ndarray
    Ty: np.ndarray
    u_max: dict[str, float]
    num_std_devs: float  # Lower and upper bound of uncertainties is determined by the number of standard deviations that we consider
    plot_vars: list[str]
    mterm: Union[str, None] = None
    lterm: Union[str, None] = None
    consider_uncertainty: bool = True
    simulate_steady_state: bool = True
    simulate_mpc: bool = True
    mpc_live_vis: bool = True
    pygame_vis: bool = False
    store_results: bool = True
    compile_nlp: bool = False
    vol_flow_rate: Union[np.ndarray, None] = None

    _state_names = [
        "S_ac",
        "S_ch4",
        "S_IC",
        "S_IN",
        "S_h2o",
        "X_ch_f",
        "X_ch_s",
        "X_pr",
        "X_li",
        "X_bac",
        "X_ac",
        "X_ash",
        "S_ion",
        "S_ac−",
        "S_hco3−",
        "S_nh3",
        "S_ch4_gas",
        "S_co2_gas",
        "V_CH4",
        "V_CO2",
    ]

    _meas_names = ["V´_g", "p_CH4", "p_CO2", "pH", "S_IN", "TS", "VS", "S_ac"]

    def setup(self):
        self._substrate_setup()
        self._setup_state_and_normalization_vectors()
        self._model_setup()

        if self.pygame_vis:
            self._pygame_setup()

        # Estimator setup
        self._estimator = state_estimator.StateEstimator(self.model)

        if self.simulate_steady_state:
            self._n_steps_steady_state = round(self.n_days_steady_state / self.t_step)

        if self.simulate_mpc:
            cost_func = utils.CostFunction(
                mterm=eval(self.mterm), lterm=eval(self.mterm)
            )

            self.mpc = mpc_setup(
                model=self.model,
                t_step=self.t_step,
                n_horizon=self.n_horizon,
                n_robust=self.n_robust,
                x_ch_in=self._x_ch_in,
                x_pr_in=self._x_pr_in,
                x_li_in=self._x_li_in,
                compile_nlp=self.compile_nlp,
                vol_flow_rate=self.vol_flow_rate,
                cost_func=cost_func,
            )

            self._n_steps_mpc = round(self.n_days_mpc / self.t_step)

        if self.mpc_live_vis:
            # Initialize graphic:
            self._graphics = do_mpc.graphics.Graphics(self.mpc.data)

            # Configure plot:
            fig, ax = plt.subplots(len(self.plot_vars), sharex=True)
            for idx, var in enumerate(self.plot_vars):
                self._graphics.add_line(
                    var_type=f"_{var[0]}", var_name=var, axis=ax[idx]
                )
                ax[idx].set_ylabel(var)

            # Update properties for all prediction lines:
            for line_i in self._graphics.pred_lines.full:
                line_i.set_linewidth(1)

            fig.align_ylabels()
            fig.tight_layout()
            plt.ion()

    def run(self):
        if self.simulate_steady_state:
            self._run_steady_state_sim()

        if self.simulate_mpc:
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
        for sub_name in self.sub_names:
            self._subs.append(getattr(substrates, sub_name))
        xi = [sub.xi for sub in self._subs]

        self._xi_norm = list(np.array([val / self.Tx[:-2] for val in xi]).T)

        uncertain_xis = [sub.get_uncertain_xi_ch_pr_li() for sub in self._subs]

        if not self.consider_uncertainty:
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

        if self.consider_uncertainty:
            self._x_ch_in = np.array(
                [
                    x_ch_nom - self.num_std_devs * x_ch_std_dev,
                    x_ch_nom + self.num_std_devs * x_ch_std_dev,
                ]
            )
            self._x_pr_in = np.array(
                [
                    x_pr_nom - self.num_std_devs * x_pr_std_dev,
                    x_pr_nom + self.num_std_devs * x_pr_std_dev,
                ]
            )
            self._x_li_in = np.array(
                [
                    x_li_nom - self.num_std_devs * x_li_std_dev,
                    x_li_nom + self.num_std_devs * x_li_std_dev,
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

    def _setup_state_and_normalization_vectors(self):
        """
        Shortens the state and normalization vectors if they've been handed too long if the gas storage had been considered in a previous simulation but not anymore.
        Also normalizes the state vector and sets up the normalization vector for the inputs.
        """
        if self.scenario_type == ScenarioType.METHANATION:
            self.x0 = self.x0[:18]
            self.Tx = self.Tx[:18]

        self.x0_norm = np.copy(self.x0)
        self.x0_norm /= self.Tx
        self.Tu = np.array([self.u_max[sub.state] for sub in self._subs])

    def _model_setup(self):
        # Model
        self.model = adm1_r3_frac_norm(
            self._xi_norm, self.Tu, self.Tx, self.Ty, self.scenario_type
        )

    def _run_steady_state_sim(self):
        # Setup the simulator
        steady_state_simulator = simulator_setup(
            model=self.model,
            t_step=self.t_step,
            x_ch_in=self._x_ch_in,
            x_pr_in=self._x_pr_in,
            x_li_in=self._x_li_in,
            vol_flow_rate=np.zeros(self._n_steps_steady_state),
        )

        # Set normalized x0
        steady_state_simulator.x0 = self.x0_norm

        # Run steady state simulation
        for _ in range(self._n_steps_steady_state):
            self._screen.fill("white")

            u_norm_steady_state = np.array([[0.1 for _ in self._subs]]).T
            # np.array(
            #     [[8.0 / self.Tu[0] for _ in range(len(self._subs))]]
            # ).T  # 1.0 / 1e4

            y_next = steady_state_simulator.make_step(u_norm_steady_state)
            self.x0_norm = self._estimator.estimate_x(y_next)

        # Save normalized x values
        np.savetxt(
            f"{self.name}_steady_state_x.csv",
            steady_state_simulator.data._x,
            delimiter=",",
        )

    def _run_mpc(self):
        mpc_simulator = simulator_setup(
            model=self.model,
            t_step=self.t_step,
            x_ch_in=self._x_ch_in,
            x_pr_in=self._x_pr_in,
            x_li_in=self._x_li_in,
            vol_flow_rate=self.vol_flow_rate,
        )

        # Initialize simulator and controller
        mpc_simulator.x0 = self.x0_norm
        self.mpc.x0 = self.x0_norm
        self.mpc.u0 = np.ones(len(self._subs)) * 0.5
        self.mpc.set_initial_guess()

        u_norm_computed = None

        # MPC
        for k in range(self._n_steps_mpc):
            # fill the screen with a color to wipe away anything from last frame
            self._screen.fill("white")

            u_norm_computed_old = copy.copy(u_norm_computed)
            u_norm_computed = self.mpc.make_step(self.x0_norm)

            u_norm_actual = np.copy(u_norm_computed)

            # Manipulate the actual feed to the biogas plant 'u_norm_actual'
            # based on the set disturbances
            if not self.disturbances.feed_computation_stuck is None:
                stuck_start_idx = self.disturbances.feed_computation_stuck[0]
                stuck_end_idx = (
                    stuck_start_idx + self.disturbances.feed_computation_stuck[1]
                )
                if k in range(stuck_start_idx, stuck_end_idx):
                    u_norm_actual = copy.copy(u_norm_computed_old)
                    u_norm_computed = copy.copy(u_norm_computed_old)

            if not self.disturbances.clogged_feeding is None:
                for sub_idx, val in self.disturbances.clogged_feeding.items():
                    if k in range(val[0], val[0] + val[1]):
                        u_norm_actual[sub_idx, 0] = 0.0

            if not self.disturbances.max_feeding_error is None:
                u_norm_actual *= (
                    1.0
                    + np.array(
                        [
                            np.random.uniform(
                                low=-1.0, high=1.0, size=u_norm_actual.shape[0]
                            )
                            * self.disturbances.max_feeding_error
                        ]
                    ).T
                )

            y_next = mpc_simulator.make_step(u_norm_actual)
            self.x0_norm = self._estimator.estimate_x(y_next)

            if not self.disturbances.state_jumps is None:
                for x_idx, val in self.disturbances.state_jumps.items():
                    if k == val[0]:
                        self.x0_norm[x_idx] += val[1]

            if self.mpc_live_vis:
                self._graphics.plot_results(t_ind=k)
                self._graphics.plot_predictions(t_ind=k)
                self._graphics.reset_axes()
                plt.show()
                plt.pause(0.01)

            if self.pygame_vis:
                vis.visualize(
                    self._bga,
                    self._data,
                    self._state_names,
                    self._meas_names,
                    self.Tx,
                    self.Ty,
                    self.x0_norm,
                    mpc_simulator,
                    u_norm_actual,
                    params_R3.V_GAS_STORAGE_MAX,
                    scenario_type=self.scenario_type,
                )

        if self.store_results:
            do_mpc.data.save_results(
                save_list=[self.mpc, mpc_simulator],
                result_name=f"{self.name}_mpc_results",
                result_path="./results/",
                overwrite=True,
            )
