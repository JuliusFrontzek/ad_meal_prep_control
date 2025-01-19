from dataclasses import dataclass
import pygame
import ad_meal_prep_control.visualization as vis
from ad_meal_prep_control.mpc import mpc_setup
import numpy as np
import do_mpc
from ad_meal_prep_control import substrates
from ad_meal_prep_control.models import adm1_r3_frac_norm
from ad_meal_prep_control.state_estimator import StateFeedback, mhe_setup
from ad_meal_prep_control.simulator import simulator_setup
from ad_meal_prep_control.simulator_plant import simulator_plant_setup
import copy
import matplotlib.pyplot as plt
from ad_meal_prep_control import params_R3
from ad_meal_prep_control.utils import (
    StateObserver,
    Scenario,
    typical_ch4_vol_flow_rate,
    Disturbances,
)
import os
from pathlib import Path
from tqdm import tqdm
import pickle
import sys
import time

np.random.seed(seed=42)


def theta_dummy_tvp_fun(theta: np.ndarray):
    theta[:, 0] = [
        params_R3.kchF,
        params_R3.kchS,
        params_R3.kpr,
        params_R3.kli,
        params_R3.kdec,
        params_R3.mu_m_ac,
        params_R3.K_S_ac,
        params_R3.K_I_nh3,
        params_R3.fracChFast,
    ]


@dataclass(kw_only=True)
class Simulation:
    scenario: Scenario

    def __post_init__(self):
        self._n_steps_steady_state = round(
            self.scenario.n_days_steady_state / self.scenario.t_step
        )
        self._n_steps_mpc = round(self.scenario.n_days_mpc / self.scenario.t_step)

        self._t_mpc = np.linspace(
            0,
            self.scenario.n_days_mpc,
            num=round(self.scenario.n_days_mpc / self.scenario.t_step),
            endpoint=True,
        )

        # Take ownership of Tx and x0 because these values are subject to change within the instances of this class
        self.Tx = np.copy(self.scenario.Tx)
        self.x0_norm_true = np.copy(self.scenario.x0_true) / self.Tx

        # Compute estimated state vector
        self.x0_norm_estimated = np.copy(self.x0_norm_true) * (
            1.0 + 0.1 * np.random.randn(self.x0_norm_true.shape[0])
        )

        if self.scenario.external_gas_storage_model:
            self._v_ch4_norm_true_0 = self.x0_norm_true[18]
            self._v_co2_norm_true_0 = self.x0_norm_true[19]
            self._v_ch4_norm_estimated_0 = self.x0_norm_estimated[18]
            self._v_co2_norm_estimated_0 = self.x0_norm_estimated[19]
            self._Tx_gas_storage_states = np.copy(self.Tx[18:])

        if self.scenario.disturbances.dictated_feeding is not None:
            self._num_dictated_subs = len(self.scenario.disturbances.dictated_feeding)
        else:
            self._num_dictated_subs = 0

        # Check if plotting setup is valid
        if len(self.scenario.plot_vars) == 0 and self.scenario.mpc_live_vis:
            raise ValueError("No variables provided for plotting!")

        for plot_var in self.scenario.plot_vars:
            if plot_var.startswith("dictated_sub_feed"):
                feed_num = int(plot_var.split("_")[-1])
                assert (
                    feed_num >= 0 and feed_num <= self._num_dictated_subs
                ), f"Plotting variable {plot_var} is incompatible with the simulation setup."

        # Set up ch4 outflow
        if self.scenario.P_el_chp is None:
            self._ch4_outflow_rate = None
        else:
            self._ch4_outflow_rate = typical_ch4_vol_flow_rate(
                max_power=self.scenario.P_el_chp,
                n_steps=self._n_steps_mpc
                + self.scenario.controller_params.mpc_n_horizon,
                t_step=self.scenario.t_step,
            )

        if not self.scenario.mpc_live_vis and not self.scenario.pygame_vis:
            self._suppress_ipopt_output = False
        else:
            self._suppress_ipopt_output = False

        self._theta = np.zeros(shape=(9, 1))

        self._mpc_computation_times_mikro_secs = np.zeros(self._n_steps_mpc)

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

        if self.scenario.pygame_vis:
            self._pygame_setup()

        if self.scenario.simulate_steady_state:
            self._sim_setup(
                ch4_outflow_rate=np.zeros(self._n_steps_steady_state),
                disturbances=Disturbances(),
            )
        else:
            self._sim_setup(
                self._ch4_outflow_rate, disturbances=self.scenario.disturbances
            )

        # Estimator setup
        self._estimator_setup()

    def run(self):
        if self.scenario.simulate_steady_state:
            self._run_steady_state_sim()

            self._set_new_Tx_x0()
            self._substrate_setup()
            self._model_setup()

            if self.scenario.external_gas_storage_model:
                self.x0_norm_true[18] = self._v_ch4_norm_true_0
                self.x0_norm_true[19] = self._v_co2_norm_true_0
                self.x0_norm_estimated[18] = self._v_ch4_norm_estimated_0
                self.x0_norm_estimated[19] = self._v_co2_norm_estimated_0
            self._sim_setup(
                self._ch4_outflow_rate, disturbances=self.scenario.disturbances
            )
            self._estimator_setup()

        if self.scenario.simulate_mpc:
            self._mpc = mpc_setup(
                model=self.model,
                t_step=self.scenario.t_step,
                n_horizon=self.scenario.controller_params.mpc_n_horizon,
                n_robust=self.scenario.controller_params.mpc_n_robust,
                xi_ch_norm=self._xi_ch_mpc_mhe_norm,
                xi_pr_norm=self._xi_pr_mpc_mhe_norm,
                xi_li_norm=self._xi_li_mpc_mhe_norm,
                compile_nlp=self.scenario.compile_nlp,
                ch4_outflow_rate=self._ch4_outflow_rate,
                cost_func=self.scenario.controller_params.cost_func,
                substrate_costs=[sub.cost for sub in self._subs_controlled],
                substrate_cost_formulation=self.scenario.controller_params.substrate_cost_formulation,
                store_full_solution=True,
                disturbances=self.scenario.disturbances,
                gas_storage_bound_fraction=self.scenario.controller_params.gas_storage_bound_fraction,
                bounds=self.scenario.controller_params.bounds,
                nl_cons=self.scenario.controller_params.nl_cons,
                rterm=self.scenario.controller_params.rterm,
                ch4_set_point_function=self.scenario.controller_params.ch4_set_point_function,
                theta=self._theta,
                suppress_ipopt_output=self._suppress_ipopt_output,
            )

            if self.scenario.mpc_live_vis:
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
        self._subs_controlled = []
        for sub_name in self.scenario.sub_names:
            self._subs_controlled.append(getattr(substrates, sub_name))

        self._subs_all = copy.copy(self._subs_controlled)

        if self._num_dictated_subs > 0:
            for sub_name in self.scenario.disturbances.dictated_feeding.keys():
                self._subs_all.append(getattr(substrates, sub_name))

        # Set substrate feeding limits
        self._limited_subs_indices = []
        if isinstance(self.scenario.limited_substrates, list):
            for lim_sub in self.scenario.limited_substrates:
                for idx, sub in enumerate(self._subs_controlled):
                    if sub.name == lim_sub.name:
                        sub.set_limit(
                            lim_sub.amount_remaining,
                            lim_sub.days_remaining,
                            self.scenario.t_step,
                        )
                        self._limited_subs_indices.append(idx)

        xi = [sub.xi for sub in self._subs_all]

        # Normalize xi values (also includes ch, pr and li but these are ignored in the model and later replaced with the uncertain ones)
        self._xi_norm = list(np.array([val / self.Tx[:18] for val in xi]).T)

        # Compute the uncertain xis (ch, pr and li)
        uncertain_xis = [sub.get_uncertain_xi_ch_pr_li() for sub in self._subs_all]

        # Get nominal values and standard deviations for uncertain xis of all substrates
        xi_ch_nom = np.array([un_xi[0].nominal_value for un_xi in uncertain_xis])
        xi_pr_nom = np.array([un_xi[1].nominal_value for un_xi in uncertain_xis])
        xi_li_nom = np.array([un_xi[2].nominal_value for un_xi in uncertain_xis])

        xi_ch_std_dev = np.array([un_xi[0].std_dev for un_xi in uncertain_xis])
        xi_pr_std_dev = np.array([un_xi[1].std_dev for un_xi in uncertain_xis])
        xi_li_std_dev = np.array([un_xi[2].std_dev for un_xi in uncertain_xis])

        if (
            self.scenario.controller_params.num_std_devs == 0.0
            or self.scenario.controller_params.mpc_n_robust == 0
            or not self.scenario.mismatch
        ):
            self._xi_ch_mpc_mhe = np.array([xi_ch_nom])
            self._xi_pr_mpc_mhe = np.array([xi_pr_nom])
            self._xi_li_mpc_mhe = np.array([xi_li_nom])
        else:
            self._xi_ch_mpc_mhe = np.array(
                [
                    xi_ch_nom
                    - self.scenario.controller_params.num_std_devs * xi_ch_std_dev,
                    xi_ch_nom
                    + self.scenario.controller_params.num_std_devs * xi_ch_std_dev,
                ]
            )
            self._xi_pr_mpc_mhe = np.array(
                [
                    xi_pr_nom
                    #- self.scenario.controller_params.num_std_devs * xi_pr_std_dev,
                    #xi_pr_nom
                    #+ self.scenario.controller_params.num_std_devs * xi_pr_std_dev,
                ]
            )
            self._xi_li_mpc_mhe = np.array(
                [
                    xi_li_nom
                    #- self.scenario.controller_params.num_std_devs * xi_li_std_dev,
                    #xi_li_nom
                    #+ self.scenario.controller_params.num_std_devs * xi_li_std_dev,
                ]
            )

        if (
            self.scenario.num_std_devs_sim == 0.0
            or not self.scenario.mismatch
        ):
            self._xi_ch_sim = np.array([xi_ch_nom])
            self._xi_pr_sim = np.array([xi_pr_nom])
            self._xi_li_sim = np.array([xi_li_nom])
        else:
            if (
                self.scenario.feedback          # for regular simulations do-MPC paper
                or '_ch' in self.scenario.name  # for sensitivity analysis
            ):
                self._xi_ch_sim = np.array(
                    [
                        xi_ch_nom
                        # + np.random.choice([-1, 1], size=xi_ch_nom.size)
                        + self.scenario.num_std_devs_sim
                        * xi_ch_std_dev
                    ]
                )
                self._xi_pr_sim = np.array(
                    [
                        xi_pr_nom
                    ]
                )
                self._xi_li_sim = np.array(
                    [
                        xi_li_nom
                    ]
                )
            elif '_pr' in self.scenario.name:
                self._xi_ch_sim = np.array(
                    [
                        xi_ch_nom
                    ]
                )
                self._xi_pr_sim = np.array(
                    [
                        xi_pr_nom
                        #+ np.random.choice([-1, 1], size=xi_pr_nom.size)
                        + self.scenario.num_std_devs_sim
                        * xi_pr_std_dev
                    ]
                )
                self._xi_li_sim = np.array(
                    [
                        xi_li_nom
                    ]
                )
            elif '_li' in self.scenario.name:
                self._xi_ch_sim = np.array(
                    [
                        xi_ch_nom
                    ]
                )
                self._xi_pr_sim = np.array(
                    [
                        xi_pr_nom
                    ]
                )
                self._xi_li_sim = np.array(
                    [
                        xi_li_nom
                        #+ np.random.choice([-1, 1], size=xi_li_nom.size)
                        + self.scenario.num_std_devs_sim
                        * xi_li_std_dev
                    ]
                )

        # Normalize uncertain xi's
        self._xi_ch_mpc_mhe_norm = self._xi_ch_mpc_mhe / self.Tx[5]
        self._xi_pr_mpc_mhe_norm = self._xi_pr_mpc_mhe / self.Tx[7]
        self._xi_li_mpc_mhe_norm = self._xi_li_mpc_mhe / self.Tx[8]

        self._xi_ch_sim_norm = self._xi_ch_sim / self.Tx[5]
        self._xi_pr_sim_norm = self._xi_pr_sim / self.Tx[7]
        self._xi_li_sim_norm = self._xi_li_sim / self.Tx[8]

        assert np.all(
            self._xi_ch_mpc_mhe_norm > 0.0
        ), f"Negative xi values for carbohydrates in MPC/MHE encountered."
        assert np.all(
            self._xi_pr_mpc_mhe_norm > 0.0
        ), f"Negative xi values for proteins in MPC/MHE encountered."
        assert np.all(
            self._xi_li_mpc_mhe_norm > 0.0
        ), f"Negative xi values for lipids in MPC/MHE encountered."
        assert np.all(
            self._xi_ch_sim_norm > 0.0
        ), f"Negative xi values for carbohydrates in simulation encountered."
        assert np.all(
            self._xi_pr_sim_norm > 0.0
        ), f"Negative xi values for proteins in simulation encountered."
        assert np.all(
            self._xi_li_sim_norm > 0.0
        ), f"Negative xi values for lipids in simulation encountered."

    def _estimator_setup(self):
        if self.scenario.state_observer == StateObserver.MHE:
            self._estimator = mhe_setup(
                model=self.model,
                t_step=self.scenario.t_step,
                n_horizon=self.scenario.mhe_n_horizon,
                xi_ch_norm=self._xi_ch_mpc_mhe_norm,
                xi_pr_norm=self._xi_pr_mpc_mhe_norm,
                xi_li_norm=self._xi_li_mpc_mhe_norm,
                P_x=np.diag((self.x0_norm_estimated - self.x0_norm_true) ** 2),
                P_v=0.0001 * np.ones((8, 8)),
                ch4_outflow_rate=self._ch4_outflow_rate,
                store_full_solution=True,
                suppress_ipopt_output=self._suppress_ipopt_output,
            )

            self._estimator.x0 = np.copy(self.x0_norm_estimated)
            self._estimator.set_initial_guess()
        elif self.scenario.state_observer == StateObserver.STATEFEEDBACK:
            self._estimator = StateFeedback(self.model, self._simulator)

    def _setup_state_and_normalization_vectors(self):
        """
        Shortens the state and normalization vectors if they've been handed too long if the gas storage had been considered in a previous simulation but not anymore.
        Also normalizes the state vector and sets up the normalization vector for the inputs.
        """
        if not self.scenario.external_gas_storage_model:
            self.x0_norm_true = self.x0_norm_true[:18]
            self.x0_norm_estimated = self.x0_norm_estimated[:18]
            self.Tx = self.Tx[:18]

        self.Tu = np.array([self.scenario.u_max[sub.state] for sub in self._subs_all])
        self.scenario.Tu = self.Tu

    def _model_setup(self):
        # Model
        self.model = adm1_r3_frac_norm(
            xi_norm=self._xi_norm,
            Tu=self.Tu,
            Tx=self.Tx,
            Ty=self.scenario.Ty,
            external_gas_storage_model=self.scenario.external_gas_storage_model,
            num_dictated_subs=self._num_dictated_subs,
            limited_subs_indices=self._limited_subs_indices,
        )

    def _sim_setup(
        self, ch4_outflow_rate: np.ndarray, disturbances: Disturbances = None
    ):
        self._simulator = simulator_setup(
            model=self.model,
            t_step=self.scenario.t_step,
            xi_ch_norm=self._xi_ch_sim_norm,
            xi_pr_norm=self._xi_pr_sim_norm,
            xi_li_norm=self._xi_li_sim_norm,
            ch4_outflow_rate=ch4_outflow_rate,
            disturbances=disturbances,
            theta=self._theta,
        )

        # Set normalized x0
        self._simulator.x0 = np.copy(self.x0_norm_true)
        self._simulator.set_initial_guess()

        if not self.scenario.feedback:
            #Plant
            # Compute the uncertain xis (ch, pr and li)
            uncertain_xis = [sub.get_uncertain_xi_ch_pr_li() for sub in self._subs_all]

            # Get nominal values and standard deviations for uncertain xis of all substrates
            xi_ch_nom = np.array([un_xi[0].nominal_value for un_xi in uncertain_xis])
            xi_pr_nom = np.array([un_xi[1].nominal_value for un_xi in uncertain_xis])
            xi_li_nom = np.array([un_xi[2].nominal_value for un_xi in uncertain_xis])

            # Here we compute the output of the simulator for the nominal value, so what the plant output
            # would be if there were no mismatch

            self._simulator_plant = simulator_plant_setup(
                model=self.model,
                t_step=self.scenario.t_step,
                xi_ch_norm=xi_ch_nom/self.Tx[5],
                xi_pr_norm=xi_pr_nom/self.Tx[7],
                xi_li_norm=xi_li_nom/self.Tx[8],
                ch4_outflow_rate=ch4_outflow_rate,
                disturbances=disturbances,
                theta=self._theta,
            )

            # Set normalized x0
            self._simulator_plant.x0 = np.copy(self.x0_norm_true)
            self._simulator_plant.set_initial_guess()

    def _graphics_setup(self):
        # Initialize graphic:
        self._graphics = {}
        self._graphics["mpc_graphics"] = do_mpc.graphics.Graphics(self._mpc.data)
        self._graphics["sim_graphics"] = do_mpc.graphics.Graphics(self._simulator.data)
        if self.scenario.state_observer == StateObserver.MHE:
            self._graphics["mhe_graphics"] = do_mpc.graphics.Graphics(
                self._estimator.data
            )

        # Configure plot:
        plt.rcParams["axes.grid"] = True
        num_plots = (
            len(self.scenario.plot_vars) + 1
            if self.scenario.external_gas_storage_model
            else len(self.scenario.plot_vars)
        )
        self._fig, self._ax = plt.subplots(num_plots, sharex=True)
        if num_plots == 1:
            self._ax = [self._ax]

        for idx, var in enumerate(self.scenario.plot_vars):
            if var[0] == "u":
                self._graphics["mpc_graphics"].add_line(
                    var_type=f"_{var[0]}", var_name=var, axis=self._ax[idx]
                )
                self._ax[idx].legend(
                    labels=[
                        sub.lower().replace("_", " ") for sub in self.scenario.sub_names
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

                if self.scenario.state_observer == StateObserver.MHE:
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
                    self.scenario._state_names[int(var.split("_")[-1]) - 1]
                )

        self._ax[-1].set_xlabel("Time [d]")

        if self.scenario.state_observer == StateObserver.MHE:
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
            theta_dummy_tvp_fun(self._theta)

            if self.scenario.pygame_vis:
                self._screen.fill("white")

            u_norm_steady_state = np.array(
                [
                    [
                        1.0 / _tu / len(self._subs_controlled)
                        if sub.state == "solid"
                        else 3.0 / _tu / len(self._subs_controlled)
                        for _tu, sub in zip(self.Tu, self._subs_controlled)
                    ]
                ]
            ).T

            y_next = self._simulator.make_step(u_norm_steady_state)
            self.x0_norm_true = np.array(self._simulator.x0.master)
            self.x0_norm_estimated = self._estimator.make_step(y_next)

        # Save normalized x values
        np.savetxt(
            f"./results/{self.scenario.name}_steady_state_x.csv",
            self._simulator.data._x,
            delimiter=",",
        )

    def _run_mpc(self):
        # Initialize variables
        vg_ch4_model = self._mpc.data.prediction(('_aux', 'v_ch4_dot_tank_in'))[0][0]
        vg_model = self._mpc.data.prediction(('_aux', 'y_1'))[0][0]
        ph_model = self._mpc.data.prediction(('_aux', 'y_4'))[0][0]

        if 'Scenario_2' in self.scenario.name:
            vg_ch4_tank_model = self._mpc.data.prediction(('_x', 'x_19'))[0][0]
            vg_co2_tank_model = self._mpc.data.prediction(('_x', 'x_20'))[0][0]
            vg_tank_model = self._mpc.data.prediction(('_aux', 'v_gas_storage'))[0][0]

        # Initialize simulator and controller
        self._mpc.x0 = np.copy(self.x0_norm_estimated)
        self._mpc.u0 = np.ones(len(self._subs_controlled)) * 0.5
        self._mpc.set_initial_guess()

        u_norm_computed = None

        # MPC
        try:
            for k in tqdm(range(self._n_steps_mpc)):
                time_start_mikro_secs = time.time_ns() / 1000
                theta_dummy_tvp_fun(self._theta)
                # fill the screen with a color to wipe away anything from last frame
                if self.scenario.pygame_vis:
                    self._screen.fill("white")

                u_norm_computed_old = copy.copy(u_norm_computed)
                u_norm_computed = self._mpc.make_step(self.x0_norm_estimated)

                u_norm_actual = np.copy(u_norm_computed)

                # Manipulate the actual feed to the biogas plant 'u_norm_actual'
                # based on the set disturbances
                if not self.scenario.disturbances.feed_computation_stuck is None:
                    stuck_start_idx = self.scenario.disturbances.feed_computation_stuck[
                        0
                    ]
                    stuck_end_idx = (
                        stuck_start_idx
                        + self.scenario.disturbances.feed_computation_stuck[1]
                    )
                    if k in range(stuck_start_idx, stuck_end_idx):
                        u_norm_actual = copy.copy(u_norm_computed_old)
                        u_norm_computed = copy.copy(u_norm_computed_old)

                if not self.scenario.disturbances.clogged_feeding is None:
                    for (
                        sub_idx,
                        val,
                    ) in self.scenario.disturbances.clogged_feeding.items():
                        if k in range(val[0], val[0] + val[1]):
                            u_norm_actual[sub_idx, 0] = 0.0

                if not self.scenario.disturbances.max_feeding_error is None:
                    u_norm_actual *= (
                        1.0
                        + np.array(
                            [
                                np.random.uniform(
                                    low=-1.0, high=1.0, size=u_norm_actual.shape[0]
                                )
                                * self.scenario.disturbances.max_feeding_error
                            ]
                        ).T
                    )
                if self.scenario.feedback:
                    # Prediction of MPC using model
                    vg_ch4_model = np.column_stack((vg_ch4_model, self._mpc.data.prediction(
                        ('_aux', 'v_ch4_dot_tank_in'))[0][self.scenario.t_stp_ahead_pred]))
                    vg_model = np.column_stack((vg_model, self._mpc.data.prediction(
                        ('_aux', 'y_1'))[0][self.scenario.t_stp_ahead_pred]))
                    ph_model = np.column_stack((ph_model, self._mpc.data.prediction(
                        ('_aux', 'y_4'))[0][self.scenario.t_stp_ahead_pred]))

                    if 'Scenario_2' in self.scenario.name:
                        # Prediction of MPC using model
                        vg_ch4_tank_model = np.column_stack((vg_ch4_tank_model, self._mpc.data.prediction(
                            ('_x', 'x_19'))[0][self.scenario.t_stp_ahead_pred]))
                        vg_co2_tank_model = np.column_stack((vg_co2_tank_model, self._mpc.data.prediction(
                            ('_x', 'x_20'))[0][self.scenario.t_stp_ahead_pred]))
                        vg_tank_model = np.column_stack((vg_tank_model, self._mpc.data.prediction(
                            ('_aux', 'v_gas_storage'))[0][self.scenario.t_stp_ahead_pred]))

                if not self.scenario.feedback:
                    # Plant results
                    # y_next order: 4 substrates (CORN_SILAGE, GRASS_SILAGE, CATTLE_MANURE,
                    # SUGAR_BEET_SILAGE) + 6 measurements
                    y_next_plant = self._simulator_plant.make_step(u_norm_actual) # Plant output

                y_next = self._simulator.make_step(u_norm_actual)
                self.x0_norm_true = np.array(self._simulator.x0.master)
                self.x0_norm_estimated = self._estimator.make_step(y_next)

                if not self.scenario.disturbances.state_jumps is None:
                    for (
                        x_idx,
                        val_list,
                    ) in self.scenario.disturbances.state_jumps.items():
                        for val in val_list:
                            if k == val[0]:
                                self.x0_norm_estimated[x_idx] += val[1]

                if self.scenario.mpc_live_vis:
                    for g in self._graphics.values():
                        g.plot_results(t_ind=k)
                    self._graphics["mpc_graphics"].plot_predictions(t_ind=k)
                    self._graphics["mpc_graphics"].reset_axes()
                    self._graphics["sim_graphics"].reset_axes()
                    for idx, var in enumerate(self.scenario.plot_vars):
                        if var[0] == "y":
                            y_num = int(var.split("_")[-1])

                            self._ax[idx].scatter(
                                self._t_mpc[k],
                                self._simulator.data._y[
                                    -1,
                                    self.model.u.size
                                    + y_num
                                    + self._num_dictated_subs
                                    - 1,
                                ],
                                color="red",
                            )
                            self._ax[idx].set_ylabel(
                                self.scenario._meas_names[int(var.split("_")[-1]) - 1]
                            )
                        elif var.startswith("dictated_sub_feed"):
                            feed_num = int(var.split("_")[-1])

                            self._ax[idx].scatter(
                                self._t_mpc[k],
                                self._simulator.data._y[
                                    -1, self.model.u.size + feed_num - 1
                                ],
                                color="red",
                            )
                            self._ax[idx].set_ylabel(
                                f"Dictated substrate num. {feed_num}"
                            )
                        elif var[0] != "u" and var[0] != "x":
                            try:
                                aux_expression_idx = np.where(
                                    np.array(self.model.aux.keys()) == var
                                )[0][0]
                                self._ax[idx].scatter(
                                    self._t_mpc[k],
                                    self._simulator.data._aux[-1, aux_expression_idx],
                                    color="red",
                                )
                                self._ax[idx].set_ylabel(var)
                            except IndexError:
                                raise ValueError(
                                    f"'{var}' was detected as an aux expression. However, it was either wrongly identified or is not defined as an aux expression in the do-mpc model."
                                )

                    if self.scenario.external_gas_storage_model:
                        self._ax[-1].scatter(
                            self._t_mpc[k],
                            self._ch4_outflow_rate[k],
                            color="red",
                        )
                        self._ax[-1].set_ylabel(
                            f"CH4 outflow\nvolume flow {r'$[m^3/d]$'}"
                        )

                    plt.show()
                    plt.pause(0.01)

                time_stop_mikro_secs = time.time_ns() / 1000

                self._mpc_computation_times_mikro_secs[k] = (
                    time_stop_mikro_secs - time_start_mikro_secs
                )

                if self.scenario.pygame_vis:
                    vis.visualize(
                        self.scenario,
                        self._bga,
                        self._data,
                        self.x0_norm_true,
                        self.Tx,
                        self._simulator,
                        u_norm_actual,
                        model=self.model,
                    )

            if self.scenario.feedback:
                #save predicted data
                if 'Scenario_1' in self.scenario.name:
                    predicted_data = np.concatenate((vg_model, ph_model, vg_ch4_model), axis=0)
                    np.savetxt(f'./results/Predicted Data {self.scenario.name}.csv', predicted_data)
                if 'Scenario_2' in self.scenario.name:
                    predicted_data = np.concatenate((vg_ch4_tank_model, vg_co2_tank_model, vg_tank_model,
                                                   vg_ch4_model, vg_model, ph_model), axis=0)
                    np.savetxt(f'./results/Predicted Data {self.scenario.name}.csv', predicted_data)

            if not self.scenario.feedback:
                vg_ch4_model = self._simulator_plant.data['_aux', 'v_ch4_dot_tank_in']
                vg_model = self._simulator_plant.data['_aux', 'y_1']
                ph_model = self._simulator_plant.data['_aux', 'y_4']

                if 'Scenario_2' in self.scenario.name:
                    vg_ch4_tank_model = self._simulator_plant.data['_x', 'x_19']
                    vg_co2_tank_model = self._simulator_plant.data['_x', 'x_20']
                    vg_tank_model = self._simulator_plant.data['_aux', 'v_gas_storage']

                #Save plant output
                if 'Scenario_1' in self.scenario.name:
                    plant_output = np.concatenate((vg_model, ph_model, vg_ch4_model), axis=1)
                    np.savetxt(f'./results/Plant Output {self.scenario.name}.csv', plant_output)

                if 'Scenario_2' in self.scenario.name:
                    plant_output = np.concatenate((vg_ch4_tank_model, vg_co2_tank_model, vg_tank_model,
                                                   vg_ch4_model, vg_model, ph_model), axis=1)

                    np.savetxt(f'./results/Plant Output {self.scenario.name}.csv', plant_output)

        except (KeyboardInterrupt, SystemError):
            pass
        finally:
            if self.scenario.save_results:
                self._save_results()


    def _save_results(self):
        # Save do_mpc data
        do_mpc.data.save_results(
            save_list=[self._mpc, self._simulator],
            result_name=f"{self.scenario.name}_mpc_results",
            result_path="./results/",
            overwrite=True,
        )

        # Save scenario meta data
        scenario_dict = self.scenario.to_dict()
        scenario_dict["aux_var_names"] = self.model.aux.keys()

        with open(f"./results/{self.scenario.name}_scenario_meta_data.pkl", "wb") as fp:
            pickle.dump(scenario_dict, fp)

        np.savetxt(
            fname=f"./results/{self.scenario.name}_mpc_computation_times_mikro_secs.txt",
            X=self._mpc_computation_times_mikro_secs,
        )

    def _set_new_Tx_x0(self):
        """
        Sets the new Tx normalization vector based on the simulation results from the initial
        steady state simulation as well as the new x0 value.
        """
        # Get denormalized steady state state-vector
        x_steady_state_true = np.copy(self._simulator.x0.master * self.Tx)
        x_steady_state_estimated = np.copy(self.x0_denorm_estimated)

        # Set new Tx based on max absolute values of states during steady state simulation
        self.Tx = np.max(np.abs(self._simulator.data._x), axis=0)

        if self.scenario.external_gas_storage_model:
            self.Tx[18:] = self._Tx_gas_storage_states

        # Set new normalized initial state (for MPC)
        self.x0_norm_true = (x_steady_state_true.T / self.Tx).T
        self.x0_norm_estimated = (x_steady_state_estimated.T / self.Tx).T
