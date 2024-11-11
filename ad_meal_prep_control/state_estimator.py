import do_mpc
import numpy as np
from pathlib import Path
from ad_meal_prep_control.utils import Disturbances


class StateFeedback(do_mpc.estimator.Estimator):
    """
    For now this functions as a state-feedback but will facilitate the EKF
    as well as add some noise to the measured variables.
    """

    def __init__(self, model, simulator):
        super().__init__(model)
        self._simulator = simulator

    def make_step(self, y) -> np.ndarray:
        """Return the actual state of the system. The parameter y is passed as a dummy but not actually used to match the signature of the MHE 'make_step' method."""
        return np.array(self._simulator.x0.master)


def mhe_setup(
    model: do_mpc.model.Model,
    t_step: float,
    n_horizon: int,
    xi_ch_norm: np.ndarray,
    xi_pr_norm: np.ndarray,
    xi_li_norm: np.ndarray,
    P_x: np.ndarray,
    P_v: np.ndarray,
    ch4_outflow_rate: np.ndarray,
    store_full_solution: bool,
    suppress_ipopt_output: bool = False,
) -> do_mpc.estimator.MHE:
    """
    MHE setup
    """
    raise NotImplementedError("Not properly implemented at the moment. Must be fixed.")
    num_states = model._x.size
    mhe = do_mpc.estimator.MHE(
        model, p_est_list=[f"x_{i+1}" for i in range(num_states)]
    )
    setup_mhe = {
        "t_step": t_step,
        "n_horizon": n_horizon,
        "store_full_solution": store_full_solution,
        "meas_from_data": True,
        'nlpsol_opts': {"ipopt.linear_solver": "ma27", "ipopt.hsllib": "/usr/local/lib/libcoinhsl.so.2.2.5"}
    }

    if suppress_ipopt_output:
        setup_mhe["nlpsol_opts"]["ipopt.print_level"] = 0

    mhe.set_param(**setup_mhe)

    mhe.set_default_objective(P_x=P_x, P_v=P_v)

    tvp_template = mhe.get_tvp_template()

    if num_states == 20:

        def tvp_fun(t_now):
            t_now_idx = int(np.round(t_now / t_step))
            n_steps = min(t_now_idx, n_horizon)

            for k in range(n_steps):
                tvp_template["_tvp", k, "v_ch4_dot_tank_out"] = ch4_outflow_rate[
                    t_now_idx + k
                ]

            return tvp_template

    mhe.set_tvp_fun(tvp_fun)

    p_num = mhe.get_p_template()
    p_num["x_ch_in"] = np.array(
        [np.mean(xi_ch_norm[:, i]) for i in range(xi_ch_norm.shape[1])]
    )
    p_num["x_pr_in"] = np.array(
        [np.mean(xi_pr_norm[:, i]) for i in range(xi_pr_norm.shape[1])]
    )
    p_num["x_li_in"] = np.array(
        [np.mean(xi_li_norm[:, i]) for i in range(xi_li_norm.shape[1])]
    )

    def p_fun(t_now):
        return p_num

    mhe.set_p_fun(p_fun)

    # Every state but the 13th state must not be negative
    for i in range(1, num_states + 1):
        if i != 13:
            mhe.bounds["lower", "_x", f"x_{i}"] = 0.0

    mhe.setup()

    return mhe
